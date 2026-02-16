#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_DIR="$ROOT/bin"
TARGET="$BIN_DIR/multiEpiPrep"
BASHRC="$HOME/.bashrc"

if [[ ! -f "$TARGET" ]]; then
  echo "Error: $TARGET not found"
  exit 1
fi

chmod +x "$TARGET" || true
chmod +x "$BIN_DIR/ref_prep.py" || true
chmod +x "$BIN_DIR/create_barcode_fasta.py" || true

EXPORT_LINE="export PATH=\"$BIN_DIR:\$PATH\""

touch "$BASHRC"

if grep -Fxq "$EXPORT_LINE" "$BASHRC"; then
  echo "multiEpiPrep already in PATH (~/.bashrc)"
else
  {
    echo ""
    echo "# multiEpiPrep"
    echo "$EXPORT_LINE"
  } >> "$BASHRC"
  echo "Added multiEpiPrep to PATH in ~/.bashrc"
fi

echo ""
echo "Installation complete."
echo "Run the following to activate:"
echo "  source ~/.bashrc"