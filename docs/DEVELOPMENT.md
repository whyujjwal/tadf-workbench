# Development Guide

## Set up the dev environment

The project targets Python 3.10+ and uses `uv` (or `pip`) with a venv:

```bash
# from the repo root
uv venv                              # create .venv/
source .venv/bin/activate
uv pip install -e ".[dev]"           # editable install + pytest, pytest-qt
```

If you don't have `uv`, plain `pip install -e ".[dev]"` works the same way.

## Run the app

```bash
python run.py                         # empty workbench
python run.py data/demo_scan          # opens the bundled biphenyl demo
python -m tadf_workbench.main         # equivalent to run.py with no arg
tadf-workbench                        # entry-point script (after install)
```

## Run the test suite

```bash
QT_QPA_PLATFORM=offscreen pytest tests/ -v
```

The `QT_QPA_PLATFORM=offscreen` env var lets Qt tests run without a display
(CI, headless macOS). All 73 tests should pass in well under 5 seconds.

To run one file:

```bash
QT_QPA_PLATFORM=offscreen pytest tests/test_fragments.py -v
```

## Regenerate the demo scan

The bundled `data/demo_scan/*.log` files are biphenyl geometries with
synthetic S₁/T₁/TDM values. Edit `ANGLE_TO_ENERGIES` in
[`scripts/generate_demo_scan.py`](../scripts/generate_demo_scan.py) and
re-run:

```bash
python scripts/generate_demo_scan.py
```

This is the script to copy and adapt if you want to demo a different
donor–acceptor system without running real Gaussian jobs.

## Adding a new metric

The plot widgets, the data table, and the Current State card all read
through [`analysis/tadf_metrics.py`](../src/tadf_workbench/analysis/tadf_metrics.py).
To expose a new per-θ quantity:

1. Add a property on `AnglePoint` ([`core/angle_scan.py`](../src/tadf_workbench/core/angle_scan.py))
   that derives it from `self.molecule.quantum_data`.
2. Add a `*_curve(scan)` helper in `analysis/tadf_metrics.py` that returns
   `(angles, values)`.
3. Add the value to the `summary_row(point)` dict so it shows up in the
   table and CSV export.
4. (Optional) Add a new plot widget in
   [`ui/angle_scan_plots.py`](../src/tadf_workbench/ui/angle_scan_plots.py)
   following the existing `_AnglePlotBase` pattern.
5. (Optional) Add a stat chip and a label in `CurrentStateCard`.

Add tests for steps 1–3 in `test_angle_scan_core.py` and
`test_tadf_metrics.py` before touching UI.

## Adding support for a non-Gaussian log format

`parse_gaussian_file` returns a `Molecule` with a `QuantumData`. Anything
that returns the same shape will work with the rest of the stack. To add a
parser for, say, ORCA:

1. Create `src/tadf_workbench/parsers/orca_parser.py` with a function
   `parse_orca_file(path) → Molecule`.
2. Update `parsers/angle_scan_parser.py` to dispatch by extension, or add a
   new `parse_orca_folder` analogue.
3. Re-export from `parsers/__init__.py`.
4. Add fixture log snippets under `tests/fixtures/` and parser tests.

The existing `parse_angle_folder` only looks for `.log` and `.out`. Extend
the `_SUPPORTED_EXTS` set in
[`parsers/angle_scan_parser.py`](../src/tadf_workbench/parsers/angle_scan_parser.py)
to add new extensions.

## Coding conventions

- **Pure-logic modules stay pure.** `core/fragments.py` and
  `analysis/tadf_metrics.py` import only `numpy` and stdlib. Don't pull Qt
  into them — it makes them harder to test and slows imports.
- **Tests first.** Every chemistry primitive in `core/fragments.py` and
  every metric in `analysis/tadf_metrics.py` has unit tests. Add tests in
  the same commit as the implementation.
- **Comments are sparse.** Code is structured to read top-to-bottom
  without commentary. Add a comment only where the WHY is non-obvious from
  the code itself (a workaround, an algorithmic invariant, a numerical
  caveat).
- **No silent failures.** When a file fails to parse, print a `[scan]`
  warning and move on; never swallow the exception without telling the user.

## Project layout

See [ARCHITECTURE.md](ARCHITECTURE.md) for the module map and data flow.
See [CHEMISTRY.md](CHEMISTRY.md) for the science the workbench is built around.
See [USER_GUIDE.md](USER_GUIDE.md) for the end-user-facing workflow.

## Releasing

There is no release pipeline yet. To make a new version:

1. Bump `version` in [`pyproject.toml`](../pyproject.toml) and
   `__version__` in [`src/tadf_workbench/__init__.py`](../src/tadf_workbench/__init__.py).
2. Update `CHANGELOG.md` (create it if missing).
3. Tag the commit: `git tag -a v0.x.0 -m "..."`.

## Known limitations / future work

- `angle_from_filename` only matches the first integer in the filename
  stem. Wider naming conventions (decimal angles, negative angles, "dihedral=42p5"
  patterns) need richer extraction — see the Open Items section of the
  original plan in [`docs/plans/2026-04-24-tadf-angle-scan.md`](plans/2026-04-24-tadf-angle-scan.md).
- The donor/acceptor dialog uses dropdowns to pick the bond. Click-to-pick
  on the 3D scene would be more ergonomic.
- The 3D viewer snaps to the nearest scanned angle. Continuous "ghost"
  rotation between snapshots could make the UI feel more alive.
- Reading the dihedral on logs whose Standard orientation differs from the
  Z-Matrix orientation might give surprises if the user's optimized
  geometry shifts the atom indices around. We always use the last Standard
  orientation block in the log.

## License

MIT — see [LICENSE](../LICENSE).
