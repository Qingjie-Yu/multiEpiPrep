#!/usr/bin/env bash
set -euo pipefail

# =========================
# help message
# =========================
usage() {
  cat <<EOF
Usage: multiEpiPrep track -i BED_DIR -o OUT_DIR -g REF_GENOME [options]

Required:
  -i, --input        Directory containing BED files
  -o, --output       Base output directory
  -g, --genome       Reference genome name (hg38 or mm10)

Optional:
  -n, --normalized   Apply RPM normalization (default: true)
      --no-normalized  Disable RPM normalization
  -h, --help         Show this help message and exit

Output:
  OUT_DIR/*.bw

Description:
  Convert BED files to bigWig tracks:
    1) count lines in BED for total reads/fragments
    2) bedtools genomecov -bg (optionally -scale 1e6/total_reads)
    3) sort bedGraph
    4) bedGraphToBigWig

Example:
  multiEpiPrep track -i ./bed -o ./output -g hg38
  multiEpiPrep track -i ./bed -o ./output -g mm10 --no-normalized
EOF
}

# =========================
# argument setting
# =========================
BED_DIR=""
OUT_DIR=""
REF_GENOME=""
NORMALIZED="true"

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
    -n|--normalized)
      NORMALIZED="true"
      shift 1
      ;;
    --no-normalized)
      NORMALIZED="false"
      shift 1
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

if [[ "${NORMALIZED}" != "true" && "${NORMALIZED}" != "false" ]]; then
  echo "[ERROR] --normalized must be true or false (or use --no-normalized)"
  exit 1
fi

# =========================
# dependency check
# =========================
command -v bedtools >/dev/null 2>&1 || { echo "[ERROR] bedtools not found in PATH"; exit 1; }
command -v bedGraphToBigWig >/dev/null 2>&1 || { echo "[ERROR] bedGraphToBigWig not found in PATH"; exit 1; }

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
# output dirs
# =========================
BW_DIR="${OUT_DIR}"
mkdir -p "${BW_DIR}"

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

  total_reads=$(awk 'NF>0{c++} END{print c+0}' "$bed")
  if [[ "${total_reads}" -eq 0 ]]; then
    echo "[ERROR] Empty BED file: $bed"
    exit 1
  fi

  scale_opt=()
  if [[ "${NORMALIZED}" == "true" ]]; then
    scale_factor=$(awk -v n="${total_reads}" 'BEGIN{printf "%.12g", 1e6/n}')
    scale_opt=( -scale "${scale_factor}" )
  fi

  bg_path="${BW_DIR}/${prefix}.bedGraph"
  sorted_bg_path="${BW_DIR}/${prefix}.sorted.bedGraph"
  bw_path="${BW_DIR}/${prefix}.bw"

  bedtools genomecov -bg "${scale_opt[@]}" -i "$bed" -g "${CHROM_SIZES}" > "${bg_path}"
  LC_ALL=C sort -k1,1 -k2,2n "${bg_path}" > "${sorted_bg_path}"
  bedGraphToBigWig "${sorted_bg_path}" "${CHROM_SIZES}" "${bw_path}"

  rm -f "${bg_path}" "${sorted_bg_path}"
done


echo "[track] Input BED dir : ${BED_DIR}"
echo "[track] Genome        : ${REF_GENOME}"
echo "[track] Normalized    : ${NORMALIZED}"
echo "[track] Output dir    : ${BW_DIR}"
echo "[track] Done"
