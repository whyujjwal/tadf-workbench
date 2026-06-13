# TADF Workbench

Interactive analysis workbench for **donor–acceptor dihedral-angle scans** in
Thermally Activated Delayed Fluorescence (TADF) molecule design.

Given a folder of Gaussian TDDFT log files — one per scanned dihedral angle θ
between a donor and acceptor fragment — the workbench finds the angle that
minimises the singlet–triplet energy gap **ΔE(S₁ − T₁)**, the key figure of
merit for efficient reverse intersystem crossing in TADF emitters.

## At a glance

- **3D molecule viewer** with the donor and acceptor fragments visually
  distinguished, the rotatable bond highlighted, both fragment best-fit planes
  rendered as translucent quads, and the IUPAC dihedral angle drawn as a live
  arc + degree readout.
- **ΔE(S₁ − T₁) vs θ** hero plot with a configurable threshold line (default
  0.2 eV) and TADF-candidate points highlighted.
- Companion plots for the **S₁ transition dipole moment magnitude** and the
  **dominant-CI-coefficient probability weight 2·P²** vs θ.
- **Per-θ data table** with TADF rows highlighted and a one-click CSV export.
- **Smart fragment dialog**: pick the rotatable bond — the molecule
  auto-splits into donor + acceptor via graph BFS and reports the formula on
  each side. Ring bonds are rejected with a clear error.
- **Live cross-wiring**: drag the slider, click any plot point, or pick a
  table row — every other view updates in lockstep.

## Quick start

```bash
# Set up a virtualenv and install
uv venv
source .venv/bin/activate
uv pip install -e ".[dev]"

# Generate the bundled biphenyl demo (22 atoms, 9 angles)
python scripts/generate_demo_scan.py

# Launch the workbench
python run.py data/demo_scan
```

The `run.py` argument is optional. If supplied, the folder is opened
immediately and the donor/acceptor configuration dialog is shown. Without an
argument, you start with an empty workbench and use **File → Open Scan
Folder…**.

Run the test suite (73 tests):

```bash
QT_QPA_PLATFORM=offscreen pytest tests/ -v
```

## Documentation

- [docs/CHEMISTRY.md](docs/CHEMISTRY.md) — TADF science, why ΔE(S₁−T₁) ≤ 0.2 eV
  matters, what the dominant CI coefficient is, what the IUPAC dihedral measures.
- [docs/USER_GUIDE.md](docs/USER_GUIDE.md) — how to prepare a scan folder, run
  the app, configure the rotatable bond, read each panel.
- [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) — module layout, data flow,
  fragment-detection algorithm, dihedral math.
- [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md) — dev environment, running tests,
  adding new metrics or parsers.
- [docs/BUILD_WINDOWS.md](docs/BUILD_WINDOWS.md) — package a standalone Windows
  `.exe` (via GitHub Actions or locally on Windows).
- [docs/plans/2026-04-24-tadf-angle-scan.md](docs/plans/2026-04-24-tadf-angle-scan.md)
  — the original implementation plan.

## File-naming convention for scan folders

Each Gaussian `.log` or `.out` in a scan folder must contain the dihedral
angle as an integer in its filename. The first integer in the basename wins:

| Filename               | Detected θ |
|------------------------|------------|
| `10.log`               | 10°        |
| `angle_30.log`         | 30°        |
| `scan_d45deg.out`      | 45°        |
| `mol_60_opt.log`       | 60°        |
| `theta-75.log`         | 75°        |

Files without a detectable integer can still be imported one at a time via
**File → Import Log Files…**, with a prompt for the angle.

## License

MIT — see [LICENSE](LICENSE).
