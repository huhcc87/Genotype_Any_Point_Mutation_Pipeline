# Locally launch Genotype PCR-RFLP on Windows PowerShell.
# Does not change the system execution policy or upload sequence files to GitHub.
Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$root = Split-Path -Parent $PSScriptRoot
Set-Location $root

Write-Host "=== Genotype PCR-RFLP: local Windows setup ==="
Write-Host "Repository: $root"
Write-Host "Uploaded sequences are processed by the local Streamlit server."

# Python.org Windows installations usually provide the 'py' launcher.
if (-not (Get-Command py -ErrorAction SilentlyContinue)) {
    throw "Python launcher 'py' was not found. Install Python 3.11 (including the py launcher) from https://www.python.org/downloads/windows/ and rerun the setup."
}
& py -3.11 -c "import sys; assert sys.version_info[:2] == (3, 11), 'Python 3.11 required'"
if ($LASTEXITCODE -ne 0) {
    throw "Python 3.11 is not installed or not accessible via 'py -3.11'. Install it from https://www.python.org/downloads/windows/ and rerun."
}

$venvDir = Join-Path $root ".venv"
$venvPython = Join-Path $venvDir "Scripts\python.exe"
if (-not (Test-Path -LiteralPath $venvDir)) {
    Write-Host "Creating isolated Python 3.11 environment..."
    & py -3.11 -m venv $venvDir
    if ($LASTEXITCODE -ne 0) {
        throw "Could not create the Python virtual environment."
    }
}
if (-not (Test-Path -LiteralPath $venvPython)) {
    throw "An existing .venv is incomplete. Rename it (to preserve its contents), then rerun."
}
& $venvPython -c "import sys; assert sys.version_info[:2] == (3, 11), 'Expected Python 3.11 virtual environment'"
if ($LASTEXITCODE -ne 0) {
    throw "Existing .venv is not Python 3.11. Rename it, then rerun; it was not deleted."
}

Write-Host "[1/3] Installing required Python packages..."
& $venvPython -m pip install --upgrade pip
if ($LASTEXITCODE -ne 0) { throw "pip upgrade failed." }
& $venvPython -m pip install -r requirements.txt
if ($LASTEXITCODE -ne 0) { throw "Dependency installation failed." }

Write-Host "[2/3] Running automated software checks..."
& $venvPython -m pytest -q
if ($LASTEXITCODE -ne 0) { throw "Tests failed. The app was not started. Review the error above." }

Write-Host "[3/3] Starting Streamlit on this Windows PC only..."
Write-Host "Open the Local URL shown below (usually http://localhost:8501)."
Write-Host "Use Try example, Upload file, or Paste sequence."
Write-Host "Press Ctrl+C in this PowerShell window to stop the app."
& $venvPython -m streamlit run app.py --server.address 127.0.0.1
if ($LASTEXITCODE -ne 0) { throw "Streamlit exited with an error." }
