#!/usr/bin/env bash
# Run from a complete local copy of the repository.
set -euo pipefail

if [[ "$(uname -s)" != "Darwin" ]]; then
    echo "This launcher is for macOS; see docs/INSTALLATION.md." >&2
    exit 1
fi

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

echo "=== Genotype PCR-RFLP: local Mac setup ==="
echo "Folder: $ROOT"
echo "Input files are processed by this local app; it does not send files to GitHub."

if command -v python3.11 >/dev/null 2>&1; then
    PYTHON="$(command -v python3.11)"
elif command -v brew >/dev/null 2>&1; then
    echo "Installing Python 3.11 using Homebrew..."
    brew install python@3.11
    PYTHON="$(brew --prefix python@3.11)/bin/python3.11"
else
    echo "Please install Python 3.11 from https://www.python.org/downloads/macos/" >&2
    echo "Then rerun the same one-time quick-start command." >&2
    exit 1
fi

"$PYTHON" -c 'import sys; assert sys.version_info[:2] == (3, 11), "Expected Python 3.11"'

if [[ ! -d ".venv" ]]; then
    "$PYTHON" -m venv .venv
elif [[ ! -x ".venv/bin/python" ]]; then
    echo "The existing .venv is invalid; rename it and rerun, without deleting user data." >&2
    exit 1
elif ! .venv/bin/python -c 'import sys; assert sys.version_info[:2] == (3, 11)' >/dev/null 2>&1; then
    echo "The existing .venv is not Python 3.11. Rename it and rerun." >&2
    exit 1
fi

VENV_PY="$ROOT/.venv/bin/python"
echo "[1/3] Installing all packages listed in requirements.txt..."
"$VENV_PY" -m pip install --upgrade pip
"$VENV_PY" -m pip install -r requirements.txt

echo "[2/3] Running software tests..."
"$VENV_PY" -m pytest -q

echo "[3/3] Starting Streamlit (local access only)..."
echo "Open the Local URL printed by Streamlit. Use Try example or Upload file."
echo "Press Control+C to stop."
"$VENV_PY" -m streamlit run app.py --server.address 127.0.0.1
