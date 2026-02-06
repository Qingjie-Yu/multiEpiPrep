#!/usr/bin/env bash
set -euo pipefail

# =========================
# help message
# =========================
usage() {
  cat <<EOF
Usage: multiEpiPrep frag -i BAM_DIR -o OUT_DIR

Required:
  -i, --input    Directory containing BAM files
  -o, --output   Base output directory

Optional:
  -h, --help     Show this help message and exit

Output:
  OUT_DIR/*.bed

Description:
  Convert BAM files into sorted BED files by:
    1) name-sorting BAM
    2) converting to BEDPE
    3) extracting fragment coordinates
    4) sorting BED by genomic position

Example:
  multiEpiPrep frag -i ./bam -o ./output
EOF
}

# =========================
# argument setting
# =========================
BAM_DIR=""
OUT_DIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)
      BAM_DIR="$2"
      shift 2
      ;;
    -o|--output)
      OUT_DIR="$2"
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
if [[ -z "${BAM_DIR}" ]]; then
  echo "[ERROR] --input is required"
  usage
  exit 1
fi

if [[ -z "${OUT_DIR}" ]]; then
  echo "[ERROR] --output is required"
  usage
  exit 1
fi

if [[ ! -d "${BAM_DIR}" ]]; then
  echo "[ERROR] BAM directory not found: ${BAM_DIR}"
  exit 1
fi

# =========================
# output dir
# =========================
BED_DIR="${OUT_DIR}"
rm -rf "${BED_DIR}"
mkdir -p "${BED_DIR}"

# =========================
# main loop
# =========================
for bam in "${BAM_DIR}"/*.bam; do
  [[ -e "$bam" ]] || continue

  if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then
    echo "[ERROR] Missing BAM index for: $bam"
    exit 1
  fi

  prefix=$(basename "$bam" .bam)

  SORT_BAM="${BED_DIR}/${prefix}.sort.bam"
  BEDPE="${BED_DIR}/${prefix}.bedpe"
  TMP_BED="${BED_DIR}/${prefix}.tmp.bed"
  FINAL_BED="${BED_DIR}/${prefix}.bed"

  # 1. name sort BAM
  samtools sort -n "$bam" -o "$SORT_BAM"

  # 2. BAM -> BEDPE
  bedtools bamtobed -bedpe -i "$SORT_BAM" > "$BEDPE"

  # 3. BEDPE -> BED (fragment)
  awk 'BEGIN{OFS="\t"} {print $1,$2,$6}' "$BEDPE" > "$TMP_BED"

  # 4. sort BED
  sort -k1,1 -k2,2n "$TMP_BED" -o "$FINAL_BED"

  # cleanup
  rm -f "$SORT_BAM" "$BEDPE" "$TMP_BED"
done

echo "[frag] Input BAM dir : ${BAM_DIR}"
echo "[frag] Output BED   : ${BED_DIR}"
echo "[frag] Done"
