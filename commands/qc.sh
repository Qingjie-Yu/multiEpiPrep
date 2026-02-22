#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$ROOT/bin/utils.sh"

die()  { echo "[ERROR] $*" >&2; exit 1; }

# help message
usage() {
  cat <<EOF
Usage: multiEpiPrep qc -i BAM_DIR -o OUT_DIR [options]

Required:
  -i, --input       Directory containing input BAM files
  -o, --output      Output directory for QC results

Optional:
  -p, --percentile  Percentile threshold for filtering (default: 0.25)
  -e, --exclude     Comma-separated list of CRF names to exclude (default: unknown,IgG_control)
                    Exclusion is performed by matching substrings in BAM filenames. If -e/--exclude is provided without a value, the exclude list is empty.
  -h, --help        Show this help message and exit

Output:
  all_read_count.tsv
  filtered_read_count.tsv

Example:
  qc_count.sh -i ./bam -o ./qc -p 0.25
EOF
}
die() { echo "$*" >&2; exit 1; }

# argument setting
BAM_DIR=""
OUT_DIR=""
PERCENTILE="0.25"
EXCLUDE_STR="unknown,IgG_control"
EXCLUDE=()

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
    -e|--exclude)
      if [[ -n "${2:-}" && "$2" != -* ]]; then
        EXCLUDE_STR="$2"
        shift 2
      else
        EXCLUDE_STR=""
        shift 1
      fi
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

awk -v p="$PERCENTILE" 'BEGIN{
  if (p ~ /^[0-9]*\.?[0-9]+$/ && p>=0 && p<=1) exit 0;
  exit 1
}' || die "[ERROR] percentile must be a number between 0 and 1"

if [[ -n "$EXCLUDE_STR" ]]; then
  IFS=',' read -r -a EXCLUDE <<< "$EXCLUDE_STR"
else
  EXCLUDE=()
fi

# create output directory
mkdir -p "${OUT_DIR}"
ALL_TSV="${OUT_DIR}/all_read_count.tsv"
FILTERED_TSV="${OUT_DIR}/filtered_read_count.tsv"

printf "pair\tread_count\n" > "$ALL_TSV"
printf "pair\tread_count\n" > "$FILTERED_TSV"

# helper function
get_read_count() {
  local bam="$1"
  local pair count
  [[ -e "$bam" ]] || return 0
  if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then
    echo "[ERROR] Missing BAM index for: $bam (run: samtools index)"
    exit 1
  fi
  pair=$(basename "$bam" .bam)
  count=$(samtools idxstats "$bam" | awk '{sum += $3} END {print sum+0}')
  printf "%s\t%s\n" "$pair" "$count" >> "$ALL_TSV"
}

in_array() {
  local val="$1"; shift
  local e
  for e in "$@"; do
    [[ "$e" == "$val" ]] && return 0
  done
  return 1
}

# scan naming mode and detect empty files
shopt -s nullglob
bams=( "$BAM_DIR"/*.bam )
shopt -u nullglob
(( ${#bams[@]} > 0 )) || die "[ERROR] No .bam files found in: $BAM_DIR"

drop_empty_files bams 0

has_dash=0
for bam in "${bams[@]}"; do
  base="$(basename "$bam")"
  [[ "$base" == *-* ]] && ((++has_dash))
done

mode="mixed"
if (( has_dash == ${#bams[@]} )); then
  mode="pair"
elif (( has_dash == 0 )); then
  mode="single"
fi

# count with exclude
if (( ${#EXCLUDE[@]} > 0 )); then
  if [[ "$mode" == "pair" ]]; then
    for bam in "${bams[@]}"; do
      b="$(basename "$bam" .bam)"
      t1="${b%%-*}"
      t2="${b#*-}"
      in_array "$t1" "${EXCLUDE[@]}" && continue
      in_array "$t2" "${EXCLUDE[@]}" && continue
      get_read_count "$bam"
    done
  elif [[ "$mode" == "single" ]]; then
    for bam in "${bams[@]}"; do
      b="$(basename "$bam" .bam)"
      skip=0
      for item in "${EXCLUDE[@]}"; do
        [[ "$b" == "$item" ]] && { skip=1; break; }
      done
      (( skip )) && continue
      get_read_count "$bam"
    done
  else
    echo "[WARN] Non-standard BAM naming detected (not all files match 'CRF1-CRF2.bam'). Substring matching can be overly broad (e.g., excluding 'YY' also matches 'YY1')." >&2
    for bam in "${bams[@]}"; do
      b="$(basename "$bam")"
      skip=0
      for item in "${EXCLUDE[@]}"; do
        [[ "$b" == *"$item"* ]] && { skip=1; break; }
      done
      (( skip )) && continue
      get_read_count "$bam"
    done
  fi
else
  for bam in "${bams[@]}"; do
    get_read_count "$bam"
  done
fi

# sort ALL_TSV by read_count in descending numeric order
{
  head -n 1 "$ALL_TSV"
  tail -n +2 "$ALL_TSV" | sort -k2,2nr
} > "${ALL_TSV}.tmp"

mv "${ALL_TSV}.tmp" "$ALL_TSV"

# filter by percentile
THRESHOLD="$(
  awk 'NR>1 {print $2}' "$ALL_TSV" \
    | sort -n \
    | awk -v p="$PERCENTILE" '
        {a[NR]=$1}
        END {
          if (NR==0) {print 0; exit 0}
          idx = int((NR-1)*p) + 1
          print a[idx]
        }'
)"

awk -v T="${THRESHOLD}" 'NR==1 || $2 >= T' "${ALL_TSV}" > "${FILTERED_TSV}"

PAIR_COUNT=$(awk 'NR>1 {c++} END {print c+0}' "${FILTERED_TSV}")

echo "[QC] Input BAM dir : ${BAM_DIR}"
echo "[QC] Output dir    : ${OUT_DIR}"
echo "[QC] Percentile    : ${PERCENTILE}"
echo "[QC] Exclude       : ${EXCLUDE_STR:-<empty>}"
echo "[QC] Threshold     : ${THRESHOLD}"
echo "[QC] Retained pairs : ${PAIR_COUNT}"
echo "[QC] Done"