#!/usr/bin/env bash
set -euo pipefail

# =========================
# help message
# =========================
usage() {
  cat <<EOF
Usage: multiEpiPrep cov -i BED_DIR -o OUT_DIR -g REF_GENOME

Required:
  -i, --input    Directory containing BED files
  -o, --output   Base output directory
  -g, --genome   Reference genome name (hg38 or mm10)

Optional:
  -h, --help     Show this help message and exit

Output:
  OUT_DIR/*.bedGraph

Description:
  Convert BED files to bedGraph coverage tracks using bedtools genomecov:
    bedtools genomecov -bg -i <BED> -g <chrom.sizes>

Example:
  multiEpiPrep cov -i ./bed -o ./output -g hg38
EOF
}

# =========================
# argument setting
# =========================
BED_DIR=""
OUT_DIR=""
REF_GENOME=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)
      BED_DIR="$2"
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
if [[ -z "${BED_DIR}" ]]; then
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

if [[ ! -d "${BED_DIR}" ]]; then
  echo "[ERROR] BED directory not found: ${BED_DIR}"
  exit 1
fi

if [[ "${REF_GENOME}" != "hg38" && "${REF_GENOME}" != "mm10" ]]; then
  echo "[ERROR] --genome must be hg38 or mm10"
  exit 1
fi

# =========================
# dependency check
# =========================
command -v bedtools >/dev/null 2>&1 || { echo "[ERROR] bedtools not found in PATH"; exit 1; }

# =========================
# chrom sizes
# =========================
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CHROM_SIZES="$(
  "$ROOT/bin/ref_prep.py" get_chromsizes --genome "$REF_GENOME" \
  | tail -n 1 \
  | sed -E 's/^"(.*)"$/\1/'
)"

if [[ ! -f "${CHROM_SIZES}" ]]; then
  echo "[ERROR] chrom.sizes file not found: ${CHROM_SIZES}"
  exit 1
fi

# =========================
# output dir
# =========================
mkdir -p "${OUT_DIR}"

# =========================
# main loop
# =========================
shopt -s nullglob
bed_files=( "${BED_DIR}"/*.bed )
if [[ ${#bed_files[@]} -eq 0 ]]; then
  echo "[ERROR] No .bed files found in: ${BED_DIR}"
  exit 1
fi

for bed in "${bed_files[@]}"; do
  prefix=$(basename "$bed" .bed)
  bg_path="${OUT_DIR}/${prefix}.bedGraph"
  bedtools genomecov -bg -i "$bed" -g "${CHROM_SIZES}" > "${bg_path}"
done

echo "[cov] Input BED dir : ${BED_DIR}"
echo "[cov] Genome        : ${REF_GENOME}"
echo "[cov] Output dir    : ${OUT_DIR}"
echo "[cov] Done"
