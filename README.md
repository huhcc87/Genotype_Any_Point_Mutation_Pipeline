# Genotype — PCR-RFLP assay designer

**Locally run, research-use Streamlit application** for sequence upload → primer design → ranked restriction enzymes → predicted genotype bands and gel → complete downloadable reports.

## Start on macOS: copy ONE code block into Terminal

**Requirements:** Mac with internet access, Git, and either Python 3.11 or Homebrew. This command creates an independent installation at `~/Genotype_Pipeline_App`; it **does not modify** an existing research folder on an external drive. It checks for local changes before updating. The launcher installs dependencies, runs tests, and starts Streamlit bound only to this Mac.

```bash
(
  set -eu
  APP_DIR="$HOME/Genotype_Pipeline_App"
  BRANCH="main"
  REPO="https://github.com/huhcc87/Genotype_Any_Point_Mutation_Pipeline.git"

  if ! command -v git >/dev/null 2>&1; then
    echo "Install Git first: https://git-scm.com/downloads/mac"
    exit 1
  fi

  if [ -d "$APP_DIR/.git" ]; then
    if [ -n "$(git -C "$APP_DIR" status --porcelain)" ]; then
      echo "Existing uncommitted changes: $APP_DIR. Preserve them before updating."
      exit 1
    fi
    git -C "$APP_DIR" fetch origin "$BRANCH"
    git -C "$APP_DIR" switch "$BRANCH"
    git -C "$APP_DIR" pull --ff-only origin "$BRANCH"
  elif [ -e "$APP_DIR" ]; then
    echo "Folder already exists and is not this Git checkout: $APP_DIR"
    exit 1
  else
    git clone --branch "$BRANCH" --single-branch "$REPO" "$APP_DIR"
  fi

  bash "$APP_DIR/scripts/start_local_mac.sh"
)
```

If both Python 3.11 and Homebrew are unavailable, the launcher explains how to install Python 3.11 from [python.org](https://www.python.org/downloads/macos/) and then rerun the same command. If Git is absent, install it from [git-scm.com](https://git-scm.com/downloads/mac). The browser's **Local URL** is printed by Streamlit, typically `http://localhost:8501` (or `8502`). Keep Terminal open; press Control+C to stop.

## Start on Windows: copy ONE code block into PowerShell

**Requirements:** Windows 10/11, [Git for Windows](https://git-scm.com/downloads/win), and [Python 3.11](https://www.python.org/downloads/windows/) with the Python `py` launcher installed. Install these prerequisites once if needed. Then open **Windows PowerShell** (not Command Prompt), copy the entire block, paste, and press Enter. The command creates a separate folder under your user profile, safely updates a clean checkout, installs all Python packages, runs tests, and launches Streamlit for local access only.

```powershell
$ErrorActionPreference = "Stop"
$appDir = Join-Path $HOME "Genotype_Pipeline_App"
$repo = "https://github.com/huhcc87/Genotype_Any_Point_Mutation_Pipeline.git"

if (-not (Get-Command git -ErrorAction SilentlyContinue)) {
    throw "Git is missing. Install it from https://git-scm.com/downloads/win, then rerun."
}
if (Test-Path -LiteralPath (Join-Path $appDir ".git")) {
    $changes = & git -C $appDir status --porcelain
    if ($LASTEXITCODE -ne 0) { throw "Could not check the existing repository." }
    if ($changes) { throw "There are local changes in $appDir. Preserve them before updating." }
    & git -C $appDir fetch origin main
    if ($LASTEXITCODE -ne 0) { throw "Git fetch failed." }
    & git -C $appDir switch main
    if ($LASTEXITCODE -ne 0) { throw "Could not switch to main." }
    & git -C $appDir pull --ff-only origin main
    if ($LASTEXITCODE -ne 0) { throw "Git fast-forward update failed." }
} elseif (Test-Path -LiteralPath $appDir) {
    throw "Folder already exists and is not a Git clone: $appDir. No files were changed."
} else {
    & git clone --branch main --single-branch $repo $appDir
    if ($LASTEXITCODE -ne 0) { throw "Could not clone the repository." }
}

& powershell.exe -NoProfile -ExecutionPolicy Bypass -File (Join-Path $appDir "scripts\start_local_windows.ps1")
if ($LASTEXITCODE -ne 0) { throw "Setup or application run failed; read the error above." }
```

The process-scoped execution-policy setting is used only to run the checked-out launcher; it does **not** permanently change Windows security settings. Review the public launcher source before executing it. After startup, open the **Local URL** printed by Streamlit (usually `http://localhost:8501`). Keep PowerShell open; press **Ctrl+C** to stop. Your input sequences remain in the local app process rather than being sent to GitHub.

[Detailed Windows/macOS setup and dependencies](docs/INSTALLATION.md) · [Windows launcher source](scripts/start_local_windows.ps1)

**Further help:** [Complete installation / all dependency versions](docs/INSTALLATION.md) · [Step-by-step user guide](docs/USER_GUIDE.md) · [Troubleshooting](docs/TROUBLESHOOTING.md).

## What you can upload, analyze, and download

| Stage | Supported functionality |
| --- | --- |
| **Upload** | FASTA (`.fa`, `.fasta`, `.fna`), `.txt`, Word `.docx`, text-based `.pdf`, or pasted DNA; 10 MB max upload |
| **Analyze** | Single-base substitution in 1-based **local sequence** coordinates (example `C101T`); reference-base validation |
| **Primers** | Forward and reverse primer candidates, Tm, product size, color-highlighted DNA binding regions |
| **Enzymes** | Predicted differential restriction fragments, A–D **heuristic design** grades, ranking |
| **Gel** | Simulated WT, heterozygous and mutant digests (fragment sizes in bp); Thermo, NEB and generic ladder presets |
| **Download** | Full ZIP package plus individual PDF, Word, HTML, PNG, CSV, JSON and FASTA outputs |

### Downloads in one ZIP

`assay_report.pdf`, `assay_report.docx`, `assay_report.html`, `sequence_highlighted.html`, `virtual_gel.png`, `primer_map.png`, `all_primers.csv`, `ranked_enzymes.csv`, `expected_fragments.csv`, `analysis.json`, `input_reference.fasta`, `README_RESULTS.txt`.

### Dependencies: all installed automatically

Python **3.11**, and Python packages `biopython`, `primer3-py`, `click`, `pyyaml`, `python-docx`, `matplotlib`, `pandas`, `streamlit`, `pytest`, `pypdf`, `reportlab`. Exact supported ranges are maintained in [requirements.txt](requirements.txt) and enumerated in the [installation guide](docs/INSTALLATION.md). There is **no need to install these Python packages manually**.

### One-click example

In Streamlit's sidebar, choose **Try example** → **Analyze example**. A TTN fragment and local `C101T` substitution are included. Review all five results tabs and select **Prepare all reports** to download the ZIP.

**Limitations:** In-silico research candidate design, **not** a clinical genotype caller or validated assay. This app does not automatically map HGVS transcript/genome coordinates, verify whole-genome primer specificity, or validate enzyme conditions. Grades reflect a computational heuristic, not measured restriction enzyme efficiency. No scanned-PDF OCR. Never upload patient identifiers to a public deployment.

## More options

- [Full app user guide](docs/USER_GUIDE.md)
- [Installation and dependency list](docs/INSTALLATION.md)
- [Troubleshooting](docs/TROUBLESHOOTING.md)
- [Automated software tests](https://github.com/huhcc87/Genotype_Any_Point_Mutation_Pipeline/actions)
- [Merged Streamlit development pull request](https://github.com/huhcc87/Genotype_Any_Point_Mutation_Pipeline/pull/1)

The interactive application is available from the default `main` branch.
