# Building the Windows `.exe`

The TADF Workbench is packaged for Windows as a **one-folder PyInstaller
build**: a `TADF-Workbench.exe` launcher plus a folder of bundled Qt / OpenGL /
numpy runtime. The target machine needs **no Python install**.

> **Why not just build it on a Mac?** PyInstaller does not cross-compile — a
> Windows executable must be produced *on* Windows. We do that on a GitHub
> Actions Windows runner (no local Windows machine required), but the same spec
> also builds locally on any Windows box.

## Option A — Build in the cloud with GitHub Actions (recommended)

The workflow at [`.github/workflows/build-windows.yml`](../.github/workflows/build-windows.yml)
runs on a `windows-latest` runner.

**To get an `.exe`:**

1. Push this repo to GitHub (the workflow runs automatically on every push to
   `main`), **or** open the repo's **Actions** tab → **Build Windows EXE** →
   **Run workflow**.
2. When the run finishes (~3–5 min), open it and download the
   **`TADF-Workbench-windows`** artifact (a `.zip`) from the **Artifacts**
   section at the bottom of the run summary.
3. Unzip it on Windows and run `TADF-Workbench/TADF-Workbench.exe`.

**To cut a downloadable release:** push a version tag and the build is attached
to a GitHub Release automatically:

```bash
git tag v0.1.0
git push origin v0.1.0
```

## Option B — Build locally on a Windows machine

On a Windows box with Python 3.12:

```bat
py -3.12 -m venv .venv
.venv\Scripts\activate
pip install -e .
pip install pyinstaller

:: optional: regenerate the demo scan so the app has something to open
python scripts\generate_demo_scan.py

pyinstaller tadf-workbench.spec --noconfirm
```

The build lands in `dist\TADF-Workbench\`. Copy `data\demo_scan` next to the
exe if you want the bundled demo:

```bat
xcopy /E /I data\demo_scan dist\TADF-Workbench\demo_scan
```

Then run `dist\TADF-Workbench\TADF-Workbench.exe`.

## Using the packaged app

Double-clicking the exe opens an empty workbench. Load data via
**File → Open Scan Folder…** and point it at the bundled `demo_scan` folder (or
any folder of Gaussian `.log` files following the
[file-naming convention](../README.md#file-naming-convention-for-scan-folders)).

## Troubleshooting

- **Want to see a crash traceback?** Edit [`tadf-workbench.spec`](../tadf-workbench.spec),
  set `console=False` → `console=True`, and rebuild. A console window will then
  show Python errors at startup.
- **3D viewer fails to open / OpenGL errors.** The spec already force-bundles all
  of PyOpenGL via `collect_submodules("OpenGL")`. On very old or headless
  Windows VMs without a GPU/driver, install a software OpenGL fallback (e.g.
  Mesa `opengl32.dll`) next to the exe.
- **App picks the wrong Qt binding.** The spec excludes `PyQt6` / `PySide2` /
  `PySide6` so pyqtgraph always resolves to PyQt5. Don't add those to the build
  environment.
- **Antivirus / SmartScreen warning.** Unsigned PyInstaller exes commonly trip
  Windows SmartScreen ("Windows protected your PC"). Click *More info → Run
  anyway*, or code-sign the exe for distribution.
