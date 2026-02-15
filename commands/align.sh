#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$ROOT/bin/utils.sh"

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
  --trim-5p      Number of nucleotides to trim from 5' end of each read (default: 45)
  --min-len        Min fragment length for bowtie2 -I (default: 10)
  --max-len        Max fragment length for bowtie2 -X (default: 1300)
  -p, --threads    Threads for bowtie2 (default: auto-detect)
  --java-mem       Java memory for Picard (default: auto-detect)
  --picard-jar     Path to picard.jar (default: auto-detect)
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
# help function
# =========================
get_picard_jar_path() {
  local picard_script current_dir link_target jar_path
  declare -g PICARD_JAR=""

  picard_script="$(command -v picard || true)"
  if [[ -z "$picard_script" || ! -e "$picard_script" ]]; then
    echo "ERROR: picard is not installed or not in PATH" >&2
    return 1
  fi

  while [[ -L "$picard_script" ]]; do
    current_dir="$(cd "$(dirname "$picard_script")" && pwd)"
    link_target="$(readlink "$picard_script")"
    if [[ "$link_target" != /* ]]; then
      picard_script="$current_dir/$link_target"
    else
      picard_script="$link_target"
    fi
  done

  current_dir="$(cd "$(dirname "$picard_script")" && pwd)"
  jar_path="$current_dir/picard.jar"

  if [[ -f "$jar_path" ]]; then
    PICARD_JAR="$jar_path"
    return 0
  else
    echo "ERROR: picard.jar not found at expected location: $jar_path" >&2
    return 1
  fi
}

# =========================
# argument setting
# =========================
FASTQ_DIR=""
OUT_DIR=""
REF_GENOME=""

TRIM_5P="45"
MIN_LEN="10"
MAX_LEN="1300"
THREADS=""
JAVA_MEM=""
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
    --trim-5p)
      TRIM_5P="$2"
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

if ! [[ "${TRIM_5P}" =~ ^[0-9]+$ ]] || [[ "${TRIM_5P}" -lt 0 ]]; then
  echo "[ERROR] --trim-5p must be a non-negative integer"
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

if [[ -z "$THREADS" ]] || ! [[ "${THREADS}" =~ ^[0-9]+$ ]] || [[ "${THREADS}" -lt 1 ]]; then
  echo "Automatically detecting available CPU cores"
  THREADS="$(get_cpu_cores)"
fi
if [[ -z "$JAVA_MEM" ]] || ! [[ "${JAVA_MEM}" =~ ^[0-9]+[mg]$ ]]; then
  echo "Automatically detecting available memory for Java"
  JAVA_MEM="$(get_memory)"
fi

if [[ -z "${PICARD_JAR:-}" ]] || [[ ! -f "${PICARD_JAR}" ]]; then
  echo "Automatically searching for picard.jar" >&2
  PICARD_JAR="$(
    "$ROOT/bin/ref_prep.py" get_picard_jar \
    | tail -n 1 \
    | sed -E 's/^"(.*)"$/\1/'
  )"
fi

[[ -n "${PICARD_JAR:-}" && -f "${PICARD_JAR}" ]] || {
  echo "[ERROR] picard jar not found: ${PICARD_JAR:-<empty>}" >&2
  exit 1
}

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
# helpers
# =========================
align_one_comb() {
  local r1="$1"
  local r2="$2"
  local comb="$3"
  local out_prefix="$4"

  local sam="${out_prefix}.sam"
  local clean_sorted="${out_prefix}.filtered.clean.sort.bam"
  local dup_bam="${out_prefix}.filtered.clean.dup.bam"
  local metrics="${out_prefix}_picard_metrics.txt"
  local final_bam="${out_prefix}.bam"

  echo "[align] Processing comb: ${comb}" >&2

  # 1) Bowtie2 -> SAM
  bowtie2 \
    -x "${INDEX_PREFIX}" \
    -I "${MIN_LEN}" \
    -X "${MAX_LEN}" \
    -5 "${TRIM_5P}" \
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

  # 2–5) SAM -> BAM | MAPQ>=10 | drop chrM/random/chrUn | sort
  samtools view -@ "${THREADS}" -h -q 10 "${sam}" \
    | awk 'BEGIN{OFS="\t"} /^@/ {print; next} ($3!="chrM" && $3 !~ /random/ && $3 !~ /^chrUn/) {print}' \
    | samtools view -@ "${THREADS}" -b - \
    | samtools sort -@ "${THREADS}" -o "${clean_sorted}" -

  rm -f "${sam}"

  # 6) Mark duplicates (REMOVE_DUPLICATES=true)
  java "-Xmx${JAVA_MEM}" -jar "${PICARD_JAR}" MarkDuplicates \
    --INPUT "${clean_sorted}" \
    --OUTPUT "${dup_bam}" \
    --METRICS_FILE "${metrics}" \
    --REMOVE_DUPLICATES true

  rm -f "${clean_sorted}"

  # 7) Rename to final BAM + build index
  if samtools index "${dup_bam}"; then
    mv -f "${dup_bam}" "${final_bam}"
  else
    samtools sort -@ "${THREADS}" -o "${final_bam}" "${dup_bam}"
    rm -f "${dup_bam}"
    samtools index "${final_bam}"
  fi

  echo "[align] Final BAM: ${final_bam}" >&2
}

# =========================
# discover combinations + run
# =========================
shopt -s nullglob
fastq_files=( "${FASTQ_DIR}"/*.fastq.gz )
if [[ ${#fastq_files[@]} -eq 0 ]]; then
  echo "[ERROR] No *.fastq.gz found in: ${FASTQ_DIR}" >&2
  exit 1
fi

declare -A FWD_BY_PREFIX
declare -A REV_BY_PREFIX
get_fastq_prefix fastq_files FWD_BY_PREFIX REV_BY_PREFIX ".fastq.gz" || exit $?

prefixes=( "${!FWD_BY_PREFIX[@]}" )

for comb in "${prefixes[@]}"; do
  r1="${FWD_BY_PREFIX[$comb]}"
  r2="${REV_BY_PREFIX[$comb]}"
  out_prefix="${OUT_DIR}/${comb}"
  align_one_comb "${r1}" "${r2}" "${comb}" "${out_prefix}"
done

echo "[align] Input FASTQ dir : ${FASTQ_DIR}"
echo "[align] Genome          : ${REF_GENOME}"
echo "[align] Index prefix    : ${INDEX_PREFIX}"
echo "[align] Picard jar       : ${PICARD_JAR}"
echo "[align] Output dir       : ${OUT_DIR}"
echo "[align] Trim 5' end      : ${TRIM_5P}"
echo "[align] Min fragment len: ${MIN_LEN}"
echo "[align] Max fragment len: ${MAX_LEN}"
echo "[align] Done"
