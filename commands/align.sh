#!/usr/bin/env bash
set -euo pipefail

# =========================
# help message
# =========================
usage() {
  cat <<EOF
Usage: multiEpiPrep align -i FASTQ_DIR -o OUT_DIR -g REF_GENOME

Required:
  -i, --input      Directory containing merged demultiplexed FASTQ files
                   Expected naming: <comb>.R1.fastq.gz and <comb>.R2.fastq.gz
  -o, --output     Output directory for BAM files (will be created; cleaned by default)
  -g, --genome     Reference genome name (hg38 or mm10)

Optional:
  --min-len        Min fragment length for bowtie2 -I (default: 10)
  --max-len        Max fragment length for bowtie2 -X (default: 1300)
  -p, --threads    Threads for bowtie2 (default: 8)
  --java-mem       Java memory for Picard (default: 16g)
  --picard-jar     Path to picard.jar (default: auto-detect via \$ROOT/bin/ref_prep.py get_picard_jar)
  -h, --help       Show this help message and exit

Output:
  OUT_DIR/<comb>.bam
  OUT_DIR/<comb>.bam.bai
  OUT_DIR/<comb>_picard_metrics.txt

Description:
  Align merged demultiplexed paired-end FASTQ files (per barcode combination) to a reference
  using bowtie2, then post-process:
    1) bowtie2 -> SAM
    2) SAM -> BAM
    3) Filter by MAPQ >= 10
    4) Remove chrM / random / chrUn
    5) Sort (coordinate)
    6) Picard MarkDuplicates (REMOVE_DUPLICATES=true)
    7) Sort (final) and index BAM

  Combination discovery:
    - scans FASTQ_DIR for *.R1.fastq.gz and requires matching *.R2.fastq.gz

Example:
  multiEpiPrep align -i ./merged_fastq -o ./bam_out -g hg38 -p 16 --java-mem 24g
EOF
}

# =========================
# argument setting
# =========================
FASTQ_DIR=""
OUT_DIR=""
REF_GENOME=""

MIN_LEN="10"
MAX_LEN="1300"
THREADS="8"
JAVA_MEM="16g"
PICARD_JAR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)
      FASTQ_DIR="$2"
      shift 2
      ;;
    -o|--output)
      OUT_DIR="$2"
      shift 2
      ;;
    -g|--genome)
      REF_GENOME="$2"
      shift 2
      ;;
    --min-len)
      MIN_LEN="$2"
      shift 2
      ;;
    --max-len)
      MAX_LEN="$2"
      shift 2
      ;;
    -p|--threads)
      THREADS="$2"
      shift 2
      ;;
    --java-mem)
      JAVA_MEM="$2"
      shift 2
      ;;
    --picard-jar)
      PICARD_JAR="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "[ERROR] Unknown option: $1"
      usage
      exit 1
      ;;
  esac
done

# =========================
# argument check
# =========================
if [[ -z "${FASTQ_DIR}" ]]; then
  echo "[ERROR] --input is required"
  usage
  exit 1
fi
if [[ -z "${OUT_DIR}" ]]; then
  echo "[ERROR] --output is required"
  usage
  exit 1
fi
if [[ -z "${REF_GENOME}" ]]; then
  echo "[ERROR] --genome is required"
  usage
  exit 1
fi

if [[ ! -d "${FASTQ_DIR}" ]]; then
  echo "[ERROR] FASTQ directory not found: ${FASTQ_DIR}"
  exit 1
fi

if [[ "${REF_GENOME}" != "hg38" && "${REF_GENOME}" != "mm10" ]]; then
  echo "[ERROR] --genome must be hg38 or mm10"
  exit 1
fi

if ! [[ "${MIN_LEN}" =~ ^[0-9]+$ ]] || [[ "${MIN_LEN}" -lt 0 ]]; then
  echo "[ERROR] --min-len must be a non-negative integer"
  exit 1
fi
if ! [[ "${MAX_LEN}" =~ ^[0-9]+$ ]] || [[ "${MAX_LEN}" -lt 0 ]]; then
  echo "[ERROR] --max-len must be a non-negative integer"
  exit 1
fi
if [[ "${MAX_LEN}" -lt "${MIN_LEN}" ]]; then
  echo "[ERROR] --max-len must be >= --min-len"
  exit 1
fi
if ! [[ "${THREADS}" =~ ^[0-9]+$ ]] || [[ "${THREADS}" -lt 1 ]]; then
  echo "[ERROR] --threads must be a positive integer"
  exit 1
fi

# =========================
# dependency check
# =========================
command -v bowtie2 >/dev/null 2>&1 || { echo "[ERROR] bowtie2 not found in PATH"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "[ERROR] samtools not found in PATH"; exit 1; }
command -v java >/dev/null 2>&1 || { echo "[ERROR] java not found in PATH"; exit 1; }
command -v sed >/dev/null 2>&1 || { echo "[ERROR] sed not found in PATH"; exit 1; }
command -v gzip >/dev/null 2>&1 || { echo "[ERROR] gzip not found in PATH"; exit 1; }

# =========================
# output dir
# =========================
mkdir -p "${OUT_DIR}"

# =========================
# bowtie2 index prefix (via ref_prep.py)
# =========================
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

INDEX_PREFIX="$(
  "$ROOT/bin/ref_prep.py" get_bowtie2_index --genome "$REF_GENOME" \
  | tail -n 1 \
  | sed -E 's/^"(.*)"$/\1/'
)"

# Bowtie2 index consists of multiple files like *.1.bt2 / *.rev.1.bt2 etc.
# We validate by checking one representative file.
if [[ ! -f "${INDEX_PREFIX}.1.bt2" && ! -f "${INDEX_PREFIX}.1.bt2l" ]]; then
  echo "[ERROR] bowtie2 index not found for prefix: ${INDEX_PREFIX}"
  echo "        Expected ${INDEX_PREFIX}.1.bt2 (or .1.bt2l)"
  exit 1
fi

# =========================
# picard jar (optional auto-detect)
# =========================
if [[ -z "${PICARD_JAR}" ]]; then
  PICARD_JAR="$(
    "$ROOT/bin/ref_prep.py" get_picard_jar \
    | tail -n 1 \
    | sed -E 's/^"(.*)"$/\1/'
  )"
fi

if [[ -z "${PICARD_JAR}" || ! -f "${PICARD_JAR}" ]]; then
  echo "[ERROR] picard jar not found: ${PICARD_JAR}"
  exit 1
fi

# =========================
# helpers
# =========================
align_one_comb() {
  local r1="$1"
  local r2="$2"
  local comb="$3"
  local out_prefix="$4"

  local sam="${out_prefix}.sam"
  local bam="${out_prefix}.bam"
  local filtered="${out_prefix}.filtered.bam"
  local clean="${out_prefix}.filtered.clean.bam"
  local temp_sam="${out_prefix}.temp.sam"
  local sorted_bam="${out_prefix}.sort.bam"
  local dup_bam="${out_prefix}.filtered.clean.dup.bam"
  local metrics="${out_prefix}_picard_metrics.txt"
  local dup_sort="${out_prefix}.filtered.clean.dup.sort.bam"
  local final_bam="${out_prefix}.bam"

  echo "[align] Processing comb: ${comb}"

  # 1) Bowtie2 -> SAM
  bowtie2 \
    -x "${INDEX_PREFIX}" \
    -I "${MIN_LEN}" \
    -X "${MAX_LEN}" \
    -5 45 \
    -p "${THREADS}" \
    -q \
    --local \
    --very-sensitive-local \
    --no-unal \
    --no-mixed \
    --no-discordant \
    --rg-id "${comb}" \
    --rg "SM:${comb}" \
    -1 "${r1}" \
    -2 "${r2}" \
    -S "${sam}"

  # 2) SAM -> BAM
  samtools view -bS -o "${bam}" "${sam}"
  rm -f "${sam}"

  # 3) Filter MAPQ >= 10
  samtools view -b -q 10 "${bam}" -o "${filtered}"
  rm -f "${bam}"

  # 4) Drop unwanted chromosomes (chrM, random, chrUn)
  samtools view -h "${filtered}" > "${temp_sam}"
  # delete alignments to unwanted references (keeps headers)
  sed '/chrM/d;/random/d;/chrUn/d' "${temp_sam}" | samtools view -Sb - > "${clean}"
  rm -f "${filtered}" "${temp_sam}"

  # 5) Sort by coordinate
  samtools sort "${clean}" -o "${sorted_bam}"
  rm -f "${clean}"

  # 6) Mark duplicates (remove duplicates)
  java "-Xmx${JAVA_MEM}" -jar "${PICARD_JAR}" MarkDuplicates \
    --INPUT "${sorted_bam}" \
    --OUTPUT "${dup_bam}" \
    --METRICS_FILE "${metrics}" \
    --REMOVE_DUPLICATES true
  rm -f "${sorted_bam}"

  # 7) Final sort (kept consistent with your python: samtools sort without -n)
  samtools sort "${dup_bam}" -o "${dup_sort}"
  rm -f "${dup_bam}"

  # 8) Rename to final BAM + index
  mv -f "${dup_sort}" "${final_bam}"
  samtools index "${final_bam}"

  echo "[align] Final BAM: ${final_bam}"
}

# =========================
# discover combinations + run
# =========================
shopt -s nullglob

r1_files=( "${FASTQ_DIR}"/*.R1.fastq.gz )
if [[ ${#r1_files[@]} -eq 0 ]]; then
  echo "[ERROR] No *.R1.fastq.gz found in: ${FASTQ_DIR}"
  exit 1
fi

missing=0
for r1 in "${r1_files[@]}"; do
  base="$(basename "${r1}")"
  comb="${base%.R1.fastq.gz}"
  r2="${FASTQ_DIR}/${comb}.R2.fastq.gz"

  if [[ ! -f "${r2}" ]]; then
    echo "[ERROR] Missing matching R2 for comb '${comb}': ${r2}"
    missing=1
    continue
  fi

  out_prefix="${OUT_DIR}/${comb}"
  align_one_comb "${r1}" "${r2}" "${comb}" "${out_prefix}"
done

if [[ "${missing}" -ne 0 ]]; then
  echo "[ERROR] Some combinations missing R2; see errors above."
  exit 1
fi

echo "[align] Input FASTQ dir : ${FASTQ_DIR}"
echo "[align] Genome          : ${REF_GENOME}"
echo "[align] Index prefix    : ${INDEX_PREFIX}"
echo "[align] Picard jar       : ${PICARD_JAR}"
echo "[align] Output dir       : ${OUT_DIR}"
echo "[align] Done"
