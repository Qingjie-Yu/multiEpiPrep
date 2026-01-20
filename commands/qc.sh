#!/usr/bin/env bash
set -euo pipefail

# help message
usage() {
  cat <<EOF
Usage: mutiEpiPrep qc -i BAM_DIR -o OUT_DIR [options]

Required:
  -i, --input       Directory containing input BAM files
  -o, --output      Output directory for QC results

Optional:
  -p, --percentile  Percentile threshold for filtering (default: 0.25)
  -h, --help        Show this help message and exit

Output:
  all_read_count.tsv
  filtered_read_count.tsv

Example:
  qc_count.sh -i ./bam -o ./qc -p 0.25
EOF
}

# argument setting
BAM_DIR=""
OUT_DIR=""
PERCENTILE="0.25"

while [[ $# -gt 0 ]]; do
  case $1 in
    -i|--input)
      BAM_DIR="$2"
      shift 2
      ;;
    -o|--output)
      OUT_DIR="$2"
      shift 2
      ;;
    -p|--percentile)
      PERCENTILE="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      usage
      exit 1
      ;;
  esac
done

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

if ! awk "BEGIN {exit !(${PERCENTILE} >= 0 && ${PERCENTILE} <= 1)}"; then
  echo "[ERROR] percentile must be between 0 and 1"
  exit 1
fi

# create output directory
mkdir -p "${OUT_DIR}"
ALL_TSV="${OUT_DIR}/all_read_count.tsv"
FILTERED_TSV="${OUT_DIR}/filtered_read_count.tsv"

# calculate read counts
echo -e "pair\tread_count" > "${ALL_TSV}"
for bam in "${BAM_DIR}"/*.bam; do
  [[ -e "$bam" ]] || continue

  if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then
    echo "[ERROR] Missing BAM index for: $bam (run: samtools index)"
    exit 1
  fi

  pair=$(basename "$bam" .bam)
  count=$(samtools idxstats "$bam" | awk '{sum += $3} END {print sum+0}')
  echo -e "${pair}\t${count}" >> "${ALL_TSV}"
done

# filter by percentile
THRESHOLD=$(awk 'NR>1 {print $2}' "${ALL_TSV}" \
  | sort -n \
  | awk -v p="${PERCENTILE}" '
      {a[NR]=$1}
      END {
        if (NR==0) exit 1
        idx = int((NR-1)*p) + 1
        print a[idx]
      }'
)

awk -v T="${THRESHOLD}" 'NR==1 || $2 >= T' "${ALL_TSV}" > "${FILTERED_TSV}"

echo "[QC] Input BAM dir : ${BAM_DIR}"
echo "[QC] Output dir    : ${OUT_DIR}"
echo "[QC] Percentile    : ${PERCENTILE}"
echo "[QC] Threshold     : ${THRESHOLD}"
echo "[QC] Done"