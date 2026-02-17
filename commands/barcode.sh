#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# =========================
# help message
# =========================
usage() {
  cat <<EOF
Usage: multiEpiPrep barcode -i BARCODE_FILE -o OUT_DIR

Required:
  -i, --input      Barcode spreadsheet (.xlsx/.xls/.csv/.tsv)
  -o, --output     Output base directory (barcode folder will be created inside)

Output:
  OUT_DIR/barcode/barcode_fwd.fasta
  OUT_DIR/barcode/barcode_rev.fasta

Description:
  Create forward and reverse barcode FASTA files from a spreadsheet.
  - Automatically detects target column (target/crf)
  - Automatically detects barcode column (bc/barcode/sequence/index)
  - Cleans target names:
        spaces/dots/slashes → _
        hyphens removed
  - Forward and reverse FASTA files contain identical content
  - Existing barcode directory will be replaced

Example:
  multiEpiPrep barcode \
    -i ./barcodes.xlsx \
    -o ./output
EOF
}

# =========================
# argument parsing
# =========================
INPUT=""
OUT_DIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)
      INPUT="$2"
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
if [[ -z "$INPUT" ]]; then
  echo "[ERROR] --input is required"
  usage
  exit 1
fi

if [[ -z "$OUT_DIR" ]]; then
  echo "[ERROR] --output is required"
  usage
  exit 1
fi

if [[ ! -f "$INPUT" ]]; then
  echo "[ERROR] Input file not found: $INPUT"
  exit 1
fi

# =========================
# dependency check
# =========================
command -v python3 >/dev/null 2>&1 || { echo "[ERROR] python3 not found in PATH"; exit 1; }

# =========================
# output directory
# =========================
mkdir -p "$OUT_DIR"
FWD_FASTA="${OUT_DIR}/barcode_fwd.fasta"
REV_FASTA="${OUT_DIR}/barcode_rev.fasta"

# =========================
# main function
# =========================
"$ROOT/bin/create_barcode_fasta.py" -i "$INPUT" -o "$OUT_DIR"

if [[ ! -f "$FWD_FASTA" ]]; then
  echo "[ERROR] Forward FASTA not generated: $FWD_FASTA"
  exit 1
fi

if [[ ! -f "$REV_FASTA" ]]; then
  echo "[ERROR] Reverse FASTA not generated: $REV_FASTA"
  exit 1
fi

if [[ ! -s "$FWD_FASTA" ]]; then
  echo "[ERROR] Forward FASTA is empty: $FWD_FASTA"
  exit 1
fi

if [[ ! -s "$REV_FASTA" ]]; then
  echo "[ERROR] Reverse FASTA is empty: $REV_FASTA"
  exit 1
fi

echo "[barcode] Input file : $INPUT"
echo "[barcode] Output dir : $OUT_DIR"
echo "[barcode] Done"