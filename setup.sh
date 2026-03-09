#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TARGET="$ROOT/multiEpiPrep"
BASHRC="$HOME/.bashrc"

if [[ ! -f "$TARGET" ]]; then
  echo "Error: executable not found: $TARGET"
  exit 1
fi

chmod +x "$TARGET" || true

EXPORT_LINE="export PATH=\"$ROOT:\$PATH\""

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
echo ""
echo "Then test with:"
echo "  multiEpiPrep --help"