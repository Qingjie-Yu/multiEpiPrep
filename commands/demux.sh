#!/usr/bin/env bash
set -euo pipefail

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
  -e, --error    Max barcode mismatches for cutadapt (default: 2)
  -j, --threads  Threads for cutadapt (default: 16)
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
ERROR_RATE="2"
THREADS="16"

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

# basic numeric validation
if ! [[ "${ERROR_RATE}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
  echo "[ERROR] --error must be a number (e.g. 2)"
  exit 1
fi
if ! [[ "${THREADS}" =~ ^[0-9]+$ ]] || [[ "${THREADS}" -lt 1 ]]; then
  echo "[ERROR] --threads must be a positive integer"
  exit 1
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
mkdir -p "${OUT_DIR}"

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
# remove empty fastq.gz (default)
# =========================
shopt -s nullglob
fastqs=( "${OUT_DIR}"/*.fastq.gz )
for fq in "${fastqs[@]}"; do
  if ! gzip -cd -- "${fq}" 2>/dev/null | head -n 1 | grep -q .; then
    rm -f -- "${fq}"
  fi
done
echo "[demux] removed empty FASTQ outputs"


echo "[demux] R1 input  : ${R1_FASTQ}"
echo "[demux] R2 input  : ${R2_FASTQ}"
echo "[demux] FWD fasta : ${FWD_FASTA}"
echo "[demux] REV fasta : ${REV_FASTA}"
echo "[demux] Output dir: ${OUT_DIR}"
echo "[demux] Done"
