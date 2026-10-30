# Install and run locally — one copy-and-paste setup

[Back to repository](../README.md) · [Complete user guide](USER_GUIDE.md) · [Troubleshooting](TROUBLESHOOTING.md)

## macOS: one code block (recommended, including Apple Silicon Macs)

Open **Terminal**, copy the **entire block**, paste it, and press Return. It creates a fresh dedicated folder under your home directory, downloads the correct GitHub feature branch, installs **Python 3.11** (using Homebrew if it is available and Python is missing), creates an isolated environment, installs every package from `requirements.txt`, runs the automated tests, and starts the local Streamlit app.

**Required initially:** Git, Python 3.11 **or** Homebrew (to install Python 3.11). If both Python 3.11 and Homebrew are absent, install Python 3.11 from [python.org](https://www.python.org/downloads/macos/) and paste the same block again. If Git is missing, install it from [git-scm.com](https://git-scm.com/downloads/mac) or via macOS Command Line Tools.

```bash
(
  set -eu
  APP_DIR="$HOME/Genotype_Pipeline_App"
  BRANCH="feature/streamlit-validated-pcr-rflp"
  REPO="https://github.com/huhcc87/Genotype_Any_Point_Mutation_Pipeline.git"

  if ! command -v git >/dev/null 2>&1; then
    echo "Git is required: https://git-scm.com/downloads/mac"
    exit 1
  fi

  if [ -d "$APP_DIR/.git" ]; then
    if [ -n "$(git -C "$APP_DIR" status --porcelain)" ]; then
      echo "Local changes found in $APP_DIR; preserve them before updating."
      exit 1
    fi
    git -C "$APP_DIR" fetch origin "$BRANCH"
    git -C "$APP_DIR" switch "$BRANCH"
    git -C "$APP_DIR" pull --ff-only origin "$BRANCH"
  elif [ -e "$APP_DIR" ]; then
    echo "$APP_DIR already exists but is not this Git checkout; nothing was changed."
    exit 1
  else
    git clone --branch "$BRANCH" --single-branch "$REPO" "$APP_DIR"
  fi

  bash "$APP_DIR/scripts/start_local_mac.sh"
)
```

The URL printed by Streamlit (normally `http://localhost:8501`, sometimes `8502`) is the address to open in the browser. The service binds to `127.0.0.1`, restricting access to **your own computer**. Keep the Terminal window open while using it. Press **Control+C** in Terminal to stop.

**If software tests fail, the launcher stops before opening the app.** Check [Troubleshooting](TROUBLESHOOTING.md), or open an issue using only non-sensitive sample data. The launcher does not alter your existing project at `/Volumes/Extreme_SSD/Genotype_Pipeline` because it makes a separate installation at `~/Genotype_Pipeline_App`.

### What the code installs

| Requirement | Version requirement | Purpose | Automatically installed? |
| --- | --- | --- | --- |
| Git | macOS compatible | Download/update GitHub repository | Must be present initially |
| Python | 3.11 (supported CI: 3.10, 3.11, 3.12) | Run the application | Launcher checks and may install via Homebrew |
| `biopython` | `>=1.85,<2` | Read DNA, restriction enzyme search | Yes |
| `primer3-py` | `>=2,<3` | Design PCR primers | Yes |
| `click` | `>=8,<9` | Command-line interface | Yes |
| `pyyaml` | `>=6,<7` | Enzyme preference configuration | Yes |
| `python-docx` | `>=1,<2` | Word uploads and reports | Yes |
| `matplotlib` | `>=3.7,<4` | Gel and primer-map figures | Yes |
| `pandas` | `>=2,<3` | Result tables/CSV exports | Yes |
| `streamlit` | `>=1.38,<2` | Web interface | Yes |
| `pytest` | `>=8,<10` | Automated tests | Yes |
| `pypdf` | `>=5,<7` | Extract text from supported PDFs | Yes |
| `reportlab` | `>=4,<5` | Generate PDF reports | Yes |

**You do not need to download these Python packages separately.** `pip install -r requirements.txt` installs them and their necessary transitive dependencies.

## Windows/Linux or existing Mac installation

The one-block script above is for macOS. For Windows/Linux, clone the [feature branch](https://github.com/huhcc87/Genotype_Any_Point_Mutation_Pipeline/tree/feature/streamlit-validated-pcr-rflp), create a virtual environment using Python 3.10–3.12, install `requirements.txt`, run `python -m pytest -q`, and start `python -m streamlit run app.py`. The repository's `scripts/start_local_mac.sh` intentionally checks for macOS to avoid pretending it supports other platforms.

## Upload and analyze

In the app, choose **Try example**, **Upload file**, or **Paste sequence**, then enter a substitution relative to the selected DNA fragment (e.g. `C101T`), and click **Analyze example** or **Analyze sequence**. Review the output tabs, choose a restriction enzyme and ladder, and generate the report bundle. Detailed instructions and file lists: [USER_GUIDE.md](USER_GUIDE.md).

## Important safety limits

Local Streamlit is preferable for unpublished research because it processes the selected file on your own computer. It is **not a certified system for patient or protected health data**. Predictions are *in silico* research candidates only; no patient genotype is inferred from experimental data. Do not treat enzyme grades as measured activity.
