#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$ROOT/bin/utils.sh"

# =========================
# help message
# =========================
usage() {
  cat <<EOF
Usage: multiEpiPrep demux -1 R1.fastq.gz -2 R2.fastq.gz -f FWD_BARCODES.fa -r REV_BARCODES.fa -o OUT_DIR

Required:
  -1, --r1       Input R1 FASTQ (.fastq.gz)
  -2, --r2       Input R2 FASTQ (.fastq.gz)
  -f, --fwd      Forward barcode FASTA (absolute or relative path)
  -r, --rev      Reverse barcode FASTA (absolute or relative path)
  -o, --output   Output directory (will be created if missing)

Optional:
  -e, --error    Max barcode mismatches for cutadapt (default: 0)
  -j, --threads  Threads for cutadapt (default: auto-detect)
  -h, --help     Show this help message and exit

Output:
  OUT_DIR/{name1}-{name2}.R1.fastq.gz
  OUT_DIR/{name1}-{name2}.R2.fastq.gz

Description:
  Demultiplex paired-end FASTQ files by barcode combinations using cutadapt.
  - Uses linked adapter files:
      -g ^file:<FWD_BARCODES.fa>
      -G ^file:<REV_BARCODES.fa>
  - No indels allowed: --no-indels
  - No trimming performed: --action none
  - By default removes empty output FASTQ files (0 lines after gzip decompression)

Example:
  multiEpiPrep demux \\
    -1 ./raw_R1.fastq.gz -2 ./raw_R2.fastq.gz \\
    -f ./barcodes_fwd.fa -r ./barcodes_rev.fa \\
    -o ./demux_out \\
    -e 2 -j 16
EOF
}

# =========================
# argument setting
# =========================
R1_FASTQ=""
R2_FASTQ=""
FWD_FASTA=""
REV_FASTA=""
OUT_DIR=""
ERROR_RATE="0"
THREADS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -1|--r1)
      R1_FASTQ="$2"
      shift 2
      ;;
    -2|--r2)
      R2_FASTQ="$2"
      shift 2
      ;;
    -f|--fwd)
      FWD_FASTA="$2"
      shift 2
      ;;
    -r|--rev)
      REV_FASTA="$2"
      shift 2
      ;;
    -o|--output)
      OUT_DIR="$2"
      shift 2
      ;;
    -e|--error)
      ERROR_RATE="$2"
      shift 2
      ;;
    -j|--threads)
      THREADS="$2"
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
if [[ -z "${R1_FASTQ}" ]]; then
  echo "[ERROR] --r1 is required"
  usage
  exit 1
fi
if [[ -z "${R2_FASTQ}" ]]; then
  echo "[ERROR] --r2 is required"
  usage
  exit 1
fi
if [[ -z "${FWD_FASTA}" ]]; then
  echo "[ERROR] --fwd is required"
  usage
  exit 1
fi
if [[ -z "${REV_FASTA}" ]]; then
  echo "[ERROR] --rev is required"
  usage
  exit 1
fi
if [[ -z "${OUT_DIR}" ]]; then
  echo "[ERROR] --output is required"
  usage
  exit 1
fi

if [[ ! -f "${R1_FASTQ}" ]]; then
  echo "[ERROR] R1 FASTQ not found: ${R1_FASTQ}"
  exit 1
fi
if [[ ! -f "${R2_FASTQ}" ]]; then
  echo "[ERROR] R2 FASTQ not found: ${R2_FASTQ}"
  exit 1
fi
if [[ ! -f "${FWD_FASTA}" ]]; then
  echo "[ERROR] Forward barcode FASTA not found: ${FWD_FASTA}"
  exit 1
fi
if [[ ! -f "${REV_FASTA}" ]]; then
  echo "[ERROR] Reverse barcode FASTA not found: ${REV_FASTA}"
  exit 1
fi

if ! [[ "${ERROR_RATE}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
  echo "[ERROR] --error must be a number (e.g. 2)"
  exit 1
fi

if [[ -z "$THREADS" ]] || ! [[ "${THREADS}" =~ ^[0-9]+$ ]] || [[ "${THREADS}" -lt 1 ]]; then
  echo "Automatically detecting available CPU cores"
  THREADS="$(get_cpu_cores)"
fi

# =========================
# dependency check
# =========================
command -v cutadapt >/dev/null 2>&1 || { echo "[ERROR] cutadapt not found in PATH"; exit 1; }
command -v gzip >/dev/null 2>&1 || { echo "[ERROR] gzip not found in PATH"; exit 1; }
command -v find >/dev/null 2>&1 || { echo "[ERROR] find not found in PATH"; exit 1; }

# =========================
# output dir
# =========================
if [[ -e "$OUT_DIR" ]]; then
  if [[ -d "$OUT_DIR" && -z "$(ls -A "$OUT_DIR")" ]]; then
    # exists and empty → OK
    :
  else
    echo "[ERROR] Output directory exists and is not empty: $OUT_DIR" >&2
    echo "Please remove it or choose a new one." >&2
    exit 1
  fi
else
  mkdir -p "$OUT_DIR"
fi

# =========================
# main: cutadapt demultiplex
# =========================
# Output naming follows your python logic: {name1}-{name2}.R1/R2.fastq.gz
out_r1="${OUT_DIR}/{name1}-{name2}.R1.fastq.gz"
out_r2="${OUT_DIR}/{name1}-{name2}.R2.fastq.gz"

cutadapt \
  -e "${ERROR_RATE}" \
  -j "${THREADS}" \
  --no-indels \
  --action none \
  -g "^file:${FWD_FASTA}" \
  -G "^file:${REV_FASTA}" \
  -o "${out_r1}" \
  -p "${out_r2}" \
  "${R1_FASTQ}" \
  "${R2_FASTQ}"

echo "[demux] cutadapt execution successful"

# =========================
# remove empty fastq.gz
# =========================
shopt -s nullglob
fastqs=( "${OUT_DIR}"/*.fastq.gz )
drop_empty_files fastqs 1
echo "[demux] removed empty FASTQ outputs"

# =========================
# merge symmetric fastq.gz
# =========================
declare -A FWD_BY_PREFIX
declare -A REV_BY_PREFIX

get_fastq_prefix fastqs FWD_BY_PREFIX REV_BY_PREFIX ".fastq.gz" || exit $?
prefixes=( "${!FWD_BY_PREFIX[@]}" )

for prefix in "${prefixes[@]}"; do
  if [[ -z "${FWD_BY_PREFIX[$prefix]:-}" || -z "${REV_BY_PREFIX[$prefix]:-}" ]]; then
    echo "[WARN] skip incomplete prefix: $prefix" >&2
    continue
  fi

  target1="${prefix%%-*}"
  target2="${prefix#*-}"
  if [[ "$target1" == "$target2" ]]; then
    continue
  fi
  if [[ "$target1" < "$target2" ]]; then
    continue
  fi

  rev_prefix="${target2}-${target1}"
  rev_r1="${FWD_BY_PREFIX[$prefix]}"
  rev_r2="${REV_BY_PREFIX[$prefix]}"
  
  if [[ -n "${FWD_BY_PREFIX[$rev_prefix]:-}" && -n "${REV_BY_PREFIX[$rev_prefix]:-}" ]]; then
    r1="${FWD_BY_PREFIX[$rev_prefix]}"
    r2="${REV_BY_PREFIX[$rev_prefix]}"
    tmp_r1="$(mktemp "${r1}.tmp.XXXX")"
    tmp_r2="$(mktemp "${r2}.tmp.XXXX")"

    gzip -cd "$r1" "$rev_r1" | gzip > "$tmp_r1"
    gzip -cd "$r2" "$rev_r2" | gzip > "$tmp_r2"
    mv "$tmp_r1" "$r1"
    mv "$tmp_r2" "$r2"
  else
    new_r1="${OUT_DIR}/${rev_prefix}.R1.fastq.gz"
    new_r2="${OUT_DIR}/${rev_prefix}.R2.fastq.gz"
    mv -f -- "$rev_r1" "$new_r1"
    mv -f -- "$rev_r2" "$new_r2"
  fi
done

echo "[demux] R1 input  : ${R1_FASTQ}"
echo "[demux] R2 input  : ${R2_FASTQ}"
echo "[demux] FWD fasta : ${FWD_FASTA}"
echo "[demux] REV fasta : ${REV_FASTA}"
echo "[demux] Output dir: ${OUT_DIR}"
echo "[demux] Done"
