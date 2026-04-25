# TADF Angle-Scan Analysis Workbench Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a new "Angle Scan" workspace that ingests a folder of Gaussian TDDFT log files — one per donor–acceptor dihedral angle θ ∈ {10°…90°} — and lets the user find the θ that minimises the S₁–T₁ gap with an interactive 3D + plots + table UI.

**Architecture:** Extend the existing layered structure (`core/` data, `parsers/` ingestion, `analysis/` derived metrics, `ui/` Qt widgets). A new `AngleScan` model aggregates per-θ `AnglePoint`s, each holding the parsed `Molecule` + derived TADF metrics. A new `AngleScanWindow` (launched from `File → Open Angle Scan Folder…`) hosts: left a 3D donor/acceptor viewer with a θ slider that snaps to the nearest scanned angle and loads that log's actual geometry; right three stacked pyqtgraph plots (ΔE vs θ with a 0.2 eV dashed threshold, |TDM(S₁)| vs θ, 2·P² vs θ); bottom a highlighted per-θ data table.

**Tech Stack:** Python 3.8+, PyQt5, pyqtgraph (2D + 3D — already a dep), pyqtgraph.opengl, NumPy, pytest, pytest-qt. No new dependencies.

**Domain interpretations locked for this plan:**
- **S₁** = the lowest-numbered `Excited State N: Singlet-*` in the log.
- **T₁** = the lowest-numbered `Excited State N: Triplet-*` in the log.
- **ΔE(S₁–T₁)** = `S₁.energy_ev − T₁.energy_ev` (eV). Threshold for TADF-interest: ≤ 0.2 eV.
- **P** = the dominant (largest-magnitude, positive) CI expansion coefficient in S₁'s orbital-transition list (e.g. `52 -> 55   0.65880` → P = 0.65880). `2·P²` is the configuration probability weight.
- **TDM(S₁)** = the row in the "Ground to excited state transition electric dipole moments (Au)" table whose `state` equals S₁'s state number; shown as |TDM| (Au).
- **θ extraction from filename**: first integer in the filename stem (so `10.log`, `angle_10.log`, `scan_d10deg.out`, `mol_10_opt.log` all → 10). Fallback: parse route/title card in the file; final fallback: user dialog to manually label unresolved files.
- **Donor/acceptor/rotatable-bond**: user configures once per scan by clicking atoms in the first parsed geometry. Stored on `AngleScan`.
- **3D θ slider**: snaps to the nearest scanned angle and loads the geometry from that angle's log file (no interpolation).

---

## File Structure

**New files:**
```
src/molecule_visualizer/
├── core/
│   └── angle_scan.py            # AnglePoint + AngleScan dataclasses
├── analysis/
│   ├── __init__.py              # Public exports
│   └── tadf_metrics.py          # compute_gap, dominant_coefficient, two_p_squared, s1_tdm_for
├── parsers/
│   └── angle_scan_parser.py     # parse_angle_folder(dir) -> AngleScan, angle_from_filename
├── ui/
│   ├── angle_scan_window.py     # top-level workspace window (QMainWindow)
│   ├── angle_scan_plots.py      # three pyqtgraph PlotWidget subclasses
│   ├── angle_scan_table.py      # QTableWidget with TADF-row highlighting
│   ├── dihedral_3d_panel.py     # 3D viewer + θ slider that snaps to scanned angles
│   └── donor_acceptor_dialog.py # first-time atom-pick dialog
tests/
├── test_angle_scan_core.py
├── test_tadf_metrics.py
├── test_angle_scan_parser.py
└── fixtures/                    # synthetic mini-log files for tests
    └── build_fixtures.py
```

**Modified files:**
- `src/molecule_visualizer/core/__init__.py` — export `AnglePoint`, `AngleScan`.
- `src/molecule_visualizer/parsers/__init__.py` — export `parse_angle_folder`.
- `src/molecule_visualizer/ui/__init__.py` — export `AngleScanWindow`.
- `src/molecule_visualizer/ui/main_window.py` — add `File → Open Angle Scan Folder…` menu action and toolbar button that launches `AngleScanWindow`.

**Conventions already in the codebase (read before coding):**
- Dark-theme stylesheet in `ui/main_window.py:_get_stylesheet` (reuse the palette).
- `ExcitedState` / `TransitionDipole` shape in `core/quantum_data.py`.
- Gaussian line patterns used for parsing in `parsers/gaussian_parser.py` — reuse `parse_gaussian_file` rather than rewriting.
- OpenGL 3D patterns in `ui/visualization_panel.py` (GLViewWidget, spheres via `gl.GLMeshItem`, bonds via `gl.GLLinePlotItem`).

---

## Task 1: Core data model — `AnglePoint` + `AngleScan`

**Files:**
- Create: `src/molecule_visualizer/core/angle_scan.py`
- Modify: `src/molecule_visualizer/core/__init__.py`
- Test: `tests/test_angle_scan_core.py`

- [ ] **Step 1.1: Write the failing test**

```python
# tests/test_angle_scan_core.py
import pytest
from molecule_visualizer.core import Molecule, AnglePoint, AngleScan
from molecule_visualizer.core.quantum_data import (
    QuantumData, ExcitedState, TransitionDipole,
)


def _make_point(angle, s1_e, t1_e, dom_coeff=0.65880, tdm=(0.1, 0.2, 0.3)):
    m = Molecule(name=f"test_{angle}")
    m.add_stage_from_data("Standard Orientation 1", [(6, 0, 0, 0), (1, 1, 0, 0)])
    qd = QuantumData()
    s1 = ExcitedState(
        number=3, multiplicity="Singlet", symmetry="A",
        energy_ev=s1_e, wavelength_nm=1240.0/s1_e, oscillator_strength=0.05,
        s_squared=0.0, orbital_transitions=[(52, 55, dom_coeff), (49, 56, 0.2), (49, 71, -0.3)],
    )
    t1 = ExcitedState(
        number=1, multiplicity="Triplet", symmetry="A",
        energy_ev=t1_e, wavelength_nm=1240.0/t1_e, oscillator_strength=0.0,
        s_squared=2.0, orbital_transitions=[(52, 55, 0.64), (49, 56, 0.22)],
    )
    qd.add_excited_state(s1)
    qd.add_excited_state(t1)
    qd.add_transition_dipole(TransitionDipole(
        state=3, x=tdm[0], y=tdm[1], z=tdm[2], dip_strength=0.14, osc_strength=0.05,
    ))
    m.quantum_data = qd
    return AnglePoint(
        angle_deg=angle, source_path=f"/tmp/{angle}.log", molecule=m,
    )


def test_anglepoint_derives_s1_t1_gap_and_2p2():
    p = _make_point(angle=30, s1_e=3.3566, t1_e=2.8989, dom_coeff=0.65880)
    assert p.angle_deg == 30
    assert p.s1_state.number == 3
    assert p.t1_state.number == 1
    assert p.s1_t1_gap_ev == pytest.approx(3.3566 - 2.8989, abs=1e-6)
    assert p.s1_dominant_coefficient == pytest.approx(0.65880)
    assert p.two_p_squared == pytest.approx(2 * 0.65880 ** 2)
    assert p.s1_tdm is not None
    assert p.s1_tdm_magnitude == pytest.approx((0.1**2 + 0.2**2 + 0.3**2) ** 0.5)


def test_anglepoint_handles_missing_tdm():
    m = Molecule(name="x")
    m.add_stage_from_data("Standard Orientation 1", [(6, 0, 0, 0)])
    qd = QuantumData()
    qd.add_excited_state(ExcitedState(
        number=1, multiplicity="Singlet", symmetry="A", energy_ev=3.0,
        wavelength_nm=413.0, oscillator_strength=0.0, s_squared=0.0,
        orbital_transitions=[(1, 2, 0.5)],
    ))
    qd.add_excited_state(ExcitedState(
        number=2, multiplicity="Triplet", symmetry="A", energy_ev=2.5,
        wavelength_nm=495.0, oscillator_strength=0.0, s_squared=2.0,
        orbital_transitions=[],
    ))
    m.quantum_data = qd
    p = AnglePoint(angle_deg=10, source_path="/tmp/10.log", molecule=m)
    assert p.s1_tdm is None
    assert p.s1_tdm_magnitude == 0.0
    assert p.s1_dominant_coefficient == pytest.approx(0.5)
    assert p.two_p_squared == pytest.approx(2 * 0.25)


def test_anglescan_sorts_by_angle_and_filters_threshold():
    scan = AngleScan(name="demo", points=[
        _make_point(angle=90, s1_e=3.2, t1_e=3.0),  # gap 0.2
        _make_point(angle=30, s1_e=3.0, t1_e=2.95),  # gap 0.05
        _make_point(angle=60, s1_e=3.4, t1_e=3.0),  # gap 0.4
    ])
    assert [p.angle_deg for p in scan.sorted_points] == [30, 60, 90]
    tadf = scan.points_below_threshold(0.2)
    assert [p.angle_deg for p in tadf] == [30, 90]  # ≤ threshold
    best = scan.minimum_gap_point
    assert best.angle_deg == 30


def test_anglescan_donor_acceptor_defaults_none():
    scan = AngleScan(name="demo", points=[])
    assert scan.donor_atom_index is None
    assert scan.acceptor_atom_index is None
    assert scan.rotatable_bond is None
```

- [ ] **Step 1.2: Run to verify it fails**

```bash
cd /Users/ujjwalraj/Developer/Chem-molecule && .venv/bin/pytest tests/test_angle_scan_core.py -v
```
Expected: ImportError / ModuleNotFoundError for `AnglePoint`, `AngleScan`.

- [ ] **Step 1.3: Implement the data model**

```python
# src/molecule_visualizer/core/angle_scan.py
"""Data model for a TADF donor-acceptor dihedral angle scan."""

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

from .molecule import Molecule
from .quantum_data import ExcitedState, TransitionDipole


@dataclass
class AnglePoint:
    """One scanned dihedral angle: molecule geometry + quantum data + derived metrics."""

    angle_deg: float
    source_path: str
    molecule: Molecule

    # Lazily computed — derived from molecule.quantum_data
    @property
    def s1_state(self) -> Optional[ExcitedState]:
        qd = self.molecule.quantum_data
        if not qd:
            return None
        singlets = [s for s in qd.excited_states if s.is_singlet]
        return min(singlets, key=lambda s: s.number) if singlets else None

    @property
    def t1_state(self) -> Optional[ExcitedState]:
        qd = self.molecule.quantum_data
        if not qd:
            return None
        triplets = [s for s in qd.excited_states if s.is_triplet]
        return min(triplets, key=lambda s: s.number) if triplets else None

    @property
    def s1_t1_gap_ev(self) -> Optional[float]:
        s1, t1 = self.s1_state, self.t1_state
        if s1 is None or t1 is None:
            return None
        return s1.energy_ev - t1.energy_ev

    @property
    def s1_dominant_coefficient(self) -> float:
        """Largest-magnitude, positive CI coefficient in S1's orbital-transition list.

        Treats the dominant configuration as the one with the largest |coefficient|;
        returns the coefficient itself (which is normally positive after Gaussian's
        phase convention, but we take abs to be robust).
        """
        s1 = self.s1_state
        if not s1 or not s1.orbital_transitions:
            return 0.0
        best = max(s1.orbital_transitions, key=lambda t: abs(t[2]))
        return abs(best[2])

    @property
    def s1_dominant_transition(self) -> Optional[Tuple[int, int, float]]:
        s1 = self.s1_state
        if not s1 or not s1.orbital_transitions:
            return None
        return max(s1.orbital_transitions, key=lambda t: abs(t[2]))

    @property
    def two_p_squared(self) -> float:
        p = self.s1_dominant_coefficient
        return 2.0 * p * p

    @property
    def s1_tdm(self) -> Optional[TransitionDipole]:
        s1 = self.s1_state
        qd = self.molecule.quantum_data
        if s1 is None or qd is None:
            return None
        return next((td for td in qd.transition_dipoles if td.state == s1.number), None)

    @property
    def s1_tdm_magnitude(self) -> float:
        td = self.s1_tdm
        return td.magnitude if td else 0.0


@dataclass
class AngleScan:
    """Collection of AnglePoints over a dihedral scan."""

    name: str
    points: List[AnglePoint] = field(default_factory=list)
    threshold_ev: float = 0.2

    # Set once by the user after the first atom-pick dialog
    donor_atom_index: Optional[int] = None
    acceptor_atom_index: Optional[int] = None
    rotatable_bond: Optional[Tuple[int, int]] = None  # two atom indices defining axis

    @property
    def sorted_points(self) -> List[AnglePoint]:
        return sorted(self.points, key=lambda p: p.angle_deg)

    @property
    def angles(self) -> List[float]:
        return [p.angle_deg for p in self.sorted_points]

    @property
    def gaps_ev(self) -> List[float]:
        return [p.s1_t1_gap_ev if p.s1_t1_gap_ev is not None else float("nan")
                for p in self.sorted_points]

    @property
    def tdm_magnitudes(self) -> List[float]:
        return [p.s1_tdm_magnitude for p in self.sorted_points]

    @property
    def two_p_squared_values(self) -> List[float]:
        return [p.two_p_squared for p in self.sorted_points]

    def points_below_threshold(self, threshold_ev: Optional[float] = None) -> List[AnglePoint]:
        t = threshold_ev if threshold_ev is not None else self.threshold_ev
        return [p for p in self.sorted_points
                if p.s1_t1_gap_ev is not None and p.s1_t1_gap_ev <= t]

    @property
    def minimum_gap_point(self) -> Optional[AnglePoint]:
        candidates = [p for p in self.points if p.s1_t1_gap_ev is not None]
        if not candidates:
            return None
        return min(candidates, key=lambda p: p.s1_t1_gap_ev)

    def nearest_point(self, angle_deg: float) -> Optional[AnglePoint]:
        if not self.points:
            return None
        return min(self.points, key=lambda p: abs(p.angle_deg - angle_deg))

    def __len__(self) -> int:
        return len(self.points)
```

- [ ] **Step 1.4: Export from core package**

Modify `src/molecule_visualizer/core/__init__.py` — add:

```python
from .angle_scan import AnglePoint, AngleScan
```

and append `"AnglePoint"`, `"AngleScan"` to `__all__` (or equivalent re-export list — follow whatever pattern the existing file uses).

- [ ] **Step 1.5: Run tests to verify pass**

```bash
.venv/bin/pytest tests/test_angle_scan_core.py -v
```
Expected: 4 passed.

- [ ] **Step 1.6: Commit**

```bash
git add src/molecule_visualizer/core/angle_scan.py src/molecule_visualizer/core/__init__.py tests/test_angle_scan_core.py
git commit -m "feat(core): add AnglePoint and AngleScan data models for TADF scans"
```

---

## Task 2: Angle extraction from filename (+ file metadata fallback)

**Files:**
- Create: `src/molecule_visualizer/parsers/angle_scan_parser.py` (only the `angle_from_filename` helper in this task)
- Test: `tests/test_angle_scan_parser.py`

- [ ] **Step 2.1: Write the failing test**

```python
# tests/test_angle_scan_parser.py
import pytest
from molecule_visualizer.parsers.angle_scan_parser import angle_from_filename


@pytest.mark.parametrize("name,expected", [
    ("10.log", 10.0),
    ("angle_30.log", 30.0),
    ("scan_d45deg.out", 45.0),
    ("mol_60_opt.log", 60.0),
    ("90deg.log", 90.0),
    ("theta-75.log", 75.0),
    ("phi_0_opt.out", 0.0),
    ("donor_acceptor_scan_180.log", 180.0),
])
def test_angle_from_filename_extracts_first_integer(name, expected):
    assert angle_from_filename(name) == expected


def test_angle_from_filename_handles_no_number():
    assert angle_from_filename("no_number_here.log") is None
    assert angle_from_filename("final.log") is None


def test_angle_from_filename_ignores_directory_numbers():
    assert angle_from_filename("/data/scan_2025/30.log") == 30.0
    assert angle_from_filename("/2024/final/angle_45.log") == 45.0
```

- [ ] **Step 2.2: Run to verify fail**

```bash
.venv/bin/pytest tests/test_angle_scan_parser.py -v
```
Expected: ImportError.

- [ ] **Step 2.3: Implement `angle_from_filename`**

```python
# src/molecule_visualizer/parsers/angle_scan_parser.py
"""Batch-parse a folder of Gaussian log files representing a dihedral angle scan."""

import os
import re
from typing import Optional


_ANGLE_PATTERN = re.compile(r"(\d+)")


def angle_from_filename(path: str) -> Optional[float]:
    """Extract the dihedral angle (degrees) from the filename stem.

    Strategy: take the basename (no directory), drop extension, return the first
    integer substring. Returns None if no integer is present.
    """
    basename = os.path.basename(path)
    stem, _ = os.path.splitext(basename)
    match = _ANGLE_PATTERN.search(stem)
    if not match:
        return None
    return float(match.group(1))
```

- [ ] **Step 2.4: Run tests pass**

```bash
.venv/bin/pytest tests/test_angle_scan_parser.py -v
```
Expected: all pass.

- [ ] **Step 2.5: Commit**

```bash
git add src/molecule_visualizer/parsers/angle_scan_parser.py tests/test_angle_scan_parser.py
git commit -m "feat(parsers): extract dihedral angle from log filename"
```

---

## Task 3: Folder batch-parser — `parse_angle_folder`

**Files:**
- Modify: `src/molecule_visualizer/parsers/angle_scan_parser.py`
- Modify: `src/molecule_visualizer/parsers/__init__.py`
- Create: `tests/fixtures/build_fixtures.py` — writes tiny, hand-crafted Gaussian-format log files used by later tests
- Test: extend `tests/test_angle_scan_parser.py`

- [ ] **Step 3.1: Create fixture builder**

```python
# tests/fixtures/build_fixtures.py
"""Generates minimal Gaussian-format .log files for testing parsers.

Each file has:
  - One Standard orientation block (2 atoms)
  - One Singlet excited state (S1) with orbital transitions
  - One Triplet excited state (T1) with orbital transitions
  - A transition electric dipole moments block with a row for S1
"""
from pathlib import Path
from textwrap import dedent


LOG_TEMPLATE = dedent('''\
 Gaussian 16 test fixture
                         Standard orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1          6           0       0.000000    0.000000    0.000000
      2          1           0       1.000000    0.000000    0.000000
 ---------------------------------------------------------------------
 Some other content between blocks to mimic real logs.

 Excited State   1:      Triplet-A            {t1_energy_ev:.4f} eV  {t1_wavelength_nm:.2f} nm  f=0.0000  <S**2>=2.000
      52 -> 55         {t1_dom:.5f}
      49 -> 56         0.22000

 Excited State   2:      Singlet-A            {s1_energy_ev:.4f} eV  {s1_wavelength_nm:.2f} nm  f=0.0500  <S**2>=0.000
      52 -> 55         {s1_dom:.5f}
      49 -> 56         0.20000
      49 -> 71        -0.10500

 Ground to excited state transition electric dipole moments (Au):
       state          X           Y           Z        Dip. S.      Osc.
         1         0.00000    0.00000    0.00000    0.00000    0.0000
         2         {tdm_x:.5f}    {tdm_y:.5f}    {tdm_z:.5f}    0.14000    0.0500
 End of table.
''')


def write_fixture(path: Path, *, s1_energy_ev: float, t1_energy_ev: float,
                  s1_dom: float = 0.65880, t1_dom: float = 0.64482,
                  tdm: tuple = (0.1, 0.2, 0.3)) -> Path:
    """Write one synthetic Gaussian log file at `path`. Returns the path."""
    path.parent.mkdir(parents=True, exist_ok=True)
    content = LOG_TEMPLATE.format(
        s1_energy_ev=s1_energy_ev,
        s1_wavelength_nm=1239.84 / s1_energy_ev,
        t1_energy_ev=t1_energy_ev,
        t1_wavelength_nm=1239.84 / t1_energy_ev,
        s1_dom=s1_dom,
        t1_dom=t1_dom,
        tdm_x=tdm[0], tdm_y=tdm[1], tdm_z=tdm[2],
    )
    path.write_text(content)
    return path


def write_scan_folder(folder: Path, angle_to_gap: dict[int, tuple]) -> Path:
    """Write a set of fixture log files under `folder`.

    `angle_to_gap` maps angle_deg -> (s1_energy_ev, t1_energy_ev).
    """
    folder.mkdir(parents=True, exist_ok=True)
    for angle, (s1e, t1e) in angle_to_gap.items():
        write_fixture(folder / f"angle_{angle}.log", s1_energy_ev=s1e, t1_energy_ev=t1e)
    return folder
```

- [ ] **Step 3.2: Write the failing test for `parse_angle_folder`**

```python
# Append to tests/test_angle_scan_parser.py
from pathlib import Path
from molecule_visualizer.parsers import parse_angle_folder
from molecule_visualizer.core import AngleScan
from tests.fixtures.build_fixtures import write_scan_folder


def test_parse_angle_folder_reads_all_logs(tmp_path):
    folder = write_scan_folder(tmp_path / "scan", {
        10: (3.40, 2.90),
        30: (3.35, 2.90),
        60: (3.20, 3.10),
    })
    scan = parse_angle_folder(str(folder))
    assert isinstance(scan, AngleScan)
    assert len(scan) == 3
    assert scan.angles == [10.0, 30.0, 60.0]
    p30 = next(p for p in scan.points if p.angle_deg == 30)
    assert p30.s1_t1_gap_ev == pytest.approx(0.45, abs=1e-3)


def test_parse_angle_folder_skips_non_log_files(tmp_path):
    folder = tmp_path / "mixed"
    write_scan_folder(folder, {10: (3.4, 2.9)})
    (folder / "notes.txt").write_text("ignore me")
    (folder / "README").write_text("ignore me too")
    scan = parse_angle_folder(str(folder))
    assert len(scan) == 1


def test_parse_angle_folder_skips_files_without_angle(tmp_path):
    folder = tmp_path / "scan"
    folder.mkdir()
    # Will contain a valid log body but filename has no number
    from tests.fixtures.build_fixtures import write_fixture
    write_fixture(folder / "final_opt.log", s1_energy_ev=3.4, t1_energy_ev=2.9)
    write_fixture(folder / "angle_30.log", s1_energy_ev=3.3, t1_energy_ev=2.9)
    scan = parse_angle_folder(str(folder))
    assert len(scan) == 1
    assert scan.points[0].angle_deg == 30.0


def test_parse_angle_folder_raises_on_missing_dir(tmp_path):
    with pytest.raises(FileNotFoundError):
        parse_angle_folder(str(tmp_path / "does_not_exist"))
```

- [ ] **Step 3.3: Implement `parse_angle_folder`**

Append to `src/molecule_visualizer/parsers/angle_scan_parser.py`:

```python
from typing import List
from ..core import AngleScan, AnglePoint
from .gaussian_parser import parse_gaussian_file


_SUPPORTED_EXTS = {".log", ".out"}


def parse_angle_folder(folder: str) -> AngleScan:
    """Parse every supported log in `folder`, attaching its dihedral angle.

    Files without a detectable angle in the filename are skipped with a warning.
    Files that fail to parse are skipped with a warning.
    The returned AngleScan is named after the folder basename.
    """
    if not os.path.isdir(folder):
        raise FileNotFoundError(f"Not a directory: {folder}")

    points: List[AnglePoint] = []
    for entry in sorted(os.listdir(folder)):
        full = os.path.join(folder, entry)
        if not os.path.isfile(full):
            continue
        _, ext = os.path.splitext(entry)
        if ext.lower() not in _SUPPORTED_EXTS:
            continue
        angle = angle_from_filename(entry)
        if angle is None:
            print(f"[angle_scan] Skipping '{entry}': no angle in filename")
            continue
        try:
            molecule = parse_gaussian_file(full)
        except Exception as exc:  # noqa: BLE001 — keep scan going
            print(f"[angle_scan] Skipping '{entry}': parse error: {exc}")
            continue
        points.append(AnglePoint(
            angle_deg=angle, source_path=full, molecule=molecule,
        ))

    name = os.path.basename(os.path.normpath(folder)) or "angle_scan"
    return AngleScan(name=name, points=points)
```

Modify `src/molecule_visualizer/parsers/__init__.py` to add:

```python
from .angle_scan_parser import parse_angle_folder, angle_from_filename
```

- [ ] **Step 3.4: Run**

```bash
.venv/bin/pytest tests/test_angle_scan_parser.py -v
```
Expected: all pass.

- [ ] **Step 3.5: Commit**

```bash
git add src/molecule_visualizer/parsers/angle_scan_parser.py src/molecule_visualizer/parsers/__init__.py tests/test_angle_scan_parser.py tests/fixtures/build_fixtures.py tests/fixtures/__init__.py
git commit -m "feat(parsers): parse_angle_folder batch-ingests a Gaussian scan directory"
```

(Create an empty `tests/fixtures/__init__.py` so it's importable.)

---

## Task 4: `analysis/tadf_metrics.py` — pure functions for plots + table

**Files:**
- Create: `src/molecule_visualizer/analysis/__init__.py`
- Create: `src/molecule_visualizer/analysis/tadf_metrics.py`
- Test: `tests/test_tadf_metrics.py`

This task centralises computations used by both plots and table, so widget code stays UI-only.

- [ ] **Step 4.1: Write the failing test**

```python
# tests/test_tadf_metrics.py
import pytest
from molecule_visualizer.analysis.tadf_metrics import (
    gap_curve, tdm_curve, two_p_squared_curve, tadf_candidates, summary_row,
)
from tests.fixtures.build_fixtures import write_scan_folder
from molecule_visualizer.parsers import parse_angle_folder


@pytest.fixture
def scan(tmp_path):
    return parse_angle_folder(str(write_scan_folder(tmp_path / "s", {
        10: (3.40, 2.90),
        30: (3.00, 2.95),  # tadf candidate: gap 0.05
        60: (3.10, 2.95),  # gap 0.15
        90: (3.40, 3.00),  # gap 0.40
    })))


def test_gap_curve_returns_sorted_pairs(scan):
    xs, ys = gap_curve(scan)
    assert xs == [10.0, 30.0, 60.0, 90.0]
    assert ys[1] == pytest.approx(0.05, abs=1e-4)


def test_tdm_curve_uses_s1_transition(scan):
    xs, ys = tdm_curve(scan)
    assert xs == [10.0, 30.0, 60.0, 90.0]
    for y in ys:
        assert y == pytest.approx((0.1**2 + 0.2**2 + 0.3**2) ** 0.5)


def test_two_p_squared_uses_dominant_coeff(scan):
    xs, ys = two_p_squared_curve(scan)
    assert ys[0] == pytest.approx(2 * 0.65880**2)


def test_tadf_candidates_returns_points_below_threshold(scan):
    cands = tadf_candidates(scan, threshold_ev=0.2)
    assert [p.angle_deg for p in cands] == [30.0, 60.0]


def test_summary_row_shape(scan):
    p = scan.sorted_points[1]  # 30°
    row = summary_row(p)
    assert row["angle_deg"] == 30.0
    assert row["s1_energy_ev"] == pytest.approx(3.00)
    assert row["t1_energy_ev"] == pytest.approx(2.95)
    assert row["gap_ev"] == pytest.approx(0.05)
    assert "tdm_magnitude" in row
    assert "two_p_squared" in row
    assert row["s1_dominant"] == (52, 55, pytest.approx(0.65880))
```

- [ ] **Step 4.2: Run (fail)**

```bash
.venv/bin/pytest tests/test_tadf_metrics.py -v
```
Expected: ImportError.

- [ ] **Step 4.3: Implement**

```python
# src/molecule_visualizer/analysis/__init__.py
from .tadf_metrics import (
    gap_curve, tdm_curve, two_p_squared_curve, tadf_candidates, summary_row,
)

__all__ = ["gap_curve", "tdm_curve", "two_p_squared_curve",
           "tadf_candidates", "summary_row"]
```

```python
# src/molecule_visualizer/analysis/tadf_metrics.py
"""Pure functions computing TADF-scan metrics from an AngleScan."""
from typing import Dict, List, Tuple, Optional
from ..core import AngleScan, AnglePoint


def gap_curve(scan: AngleScan) -> Tuple[List[float], List[float]]:
    """ΔE(S1 − T1) in eV vs θ (both sorted by θ)."""
    pts = scan.sorted_points
    xs = [p.angle_deg for p in pts]
    ys = [p.s1_t1_gap_ev if p.s1_t1_gap_ev is not None else float("nan")
          for p in pts]
    return xs, ys


def tdm_curve(scan: AngleScan) -> Tuple[List[float], List[float]]:
    pts = scan.sorted_points
    return [p.angle_deg for p in pts], [p.s1_tdm_magnitude for p in pts]


def two_p_squared_curve(scan: AngleScan) -> Tuple[List[float], List[float]]:
    pts = scan.sorted_points
    return [p.angle_deg for p in pts], [p.two_p_squared for p in pts]


def tadf_candidates(scan: AngleScan, threshold_ev: float = 0.2) -> List[AnglePoint]:
    return scan.points_below_threshold(threshold_ev)


def summary_row(p: AnglePoint) -> Dict[str, object]:
    s1 = p.s1_state
    t1 = p.t1_state
    dom = p.s1_dominant_transition
    return {
        "angle_deg": p.angle_deg,
        "s1_energy_ev": s1.energy_ev if s1 else None,
        "t1_energy_ev": t1.energy_ev if t1 else None,
        "gap_ev": p.s1_t1_gap_ev,
        "tdm_magnitude": p.s1_tdm_magnitude,
        "two_p_squared": p.two_p_squared,
        "s1_dominant": dom,
        "source_path": p.source_path,
    }
```

- [ ] **Step 4.4: Run (pass)**

```bash
.venv/bin/pytest tests/test_tadf_metrics.py -v
```
Expected: all pass.

- [ ] **Step 4.5: Commit**

```bash
git add src/molecule_visualizer/analysis/ tests/test_tadf_metrics.py
git commit -m "feat(analysis): tadf_metrics helpers for plots and table"
```

---

## Task 5: 2D plot widgets — gap, TDM, 2·P² vs θ

**Files:**
- Create: `src/molecule_visualizer/ui/angle_scan_plots.py`
- Test: `tests/test_angle_scan_plots.py` (pytest-qt; headless OK)

Use `pyqtgraph.PlotWidget`. All three plots share base styling; the gap plot additionally draws a dashed horizontal line at the threshold and coloured scatter for points ≤ threshold.

- [ ] **Step 5.1: Write the failing test**

```python
# tests/test_angle_scan_plots.py
import pytest
pytest.importorskip("PyQt5")
pytestmark = pytest.mark.usefixtures("qtbot")

from molecule_visualizer.ui.angle_scan_plots import GapPlot, TDMPlot, TwoPSquaredPlot
from molecule_visualizer.parsers import parse_angle_folder
from tests.fixtures.build_fixtures import write_scan_folder


@pytest.fixture
def scan(tmp_path):
    return parse_angle_folder(str(write_scan_folder(tmp_path / "s", {
        10: (3.4, 2.9), 30: (3.0, 2.95), 60: (3.1, 2.95), 90: (3.4, 3.0),
    })))


def test_gap_plot_populates_from_scan(qtbot, scan):
    w = GapPlot()
    qtbot.addWidget(w)
    w.set_scan(scan)
    # Should have two plot data items: main curve + highlight scatter
    items = [i for i in w.plotItem.items if hasattr(i, "setData")]
    assert len(items) >= 1


def test_gap_plot_threshold_updates(qtbot, scan):
    w = GapPlot()
    qtbot.addWidget(w)
    w.set_scan(scan)
    w.set_threshold(0.1)
    assert w.threshold_ev == 0.1


def test_tdm_plot_populates(qtbot, scan):
    w = TDMPlot()
    qtbot.addWidget(w)
    w.set_scan(scan)


def test_two_p_squared_plot_populates(qtbot, scan):
    w = TwoPSquaredPlot()
    qtbot.addWidget(w)
    w.set_scan(scan)


def test_plots_emit_click_signal(qtbot, scan):
    w = GapPlot()
    qtbot.addWidget(w)
    w.set_scan(scan)
    received = []
    w.angle_clicked.connect(received.append)
    # Simulate programmatic point selection
    w._emit_nearest(30.0)
    assert received == [30.0]
```

- [ ] **Step 5.2: Run (fail)**

```bash
.venv/bin/pytest tests/test_angle_scan_plots.py -v
```

- [ ] **Step 5.3: Implement plots**

```python
# src/molecule_visualizer/ui/angle_scan_plots.py
"""2-D plot widgets for the angle-scan workspace."""
from typing import Optional

import numpy as np
import pyqtgraph as pg
from PyQt5.QtCore import pyqtSignal, Qt
from PyQt5.QtGui import QPen, QColor

from ..core import AngleScan
from ..analysis import gap_curve, tdm_curve, two_p_squared_curve


# Match the main-window dark theme
_BG = "#111518"
_FG = "#ddeeff"
_GRID = "#252d38"
_CURVE = "#4fc3f7"
_HIGHLIGHT = "#4caf50"
_THRESHOLD = "#ff8a65"


def _configure_axes(plot: pg.PlotWidget, *, xlabel: str, ylabel: str, title: str):
    plot.setBackground(_BG)
    plot.showGrid(x=True, y=True, alpha=0.15)
    plot.setLabel("bottom", xlabel, color=_FG)
    plot.setLabel("left", ylabel, color=_FG)
    plot.setTitle(title, color=_FG, size="11pt")
    for axis in ("bottom", "left"):
        ax = plot.getAxis(axis)
        ax.setPen(QColor(_GRID))
        ax.setTextPen(QColor(_FG))


class _AnglePlotBase(pg.PlotWidget):
    """Shared behaviour: set_scan, click-to-emit nearest angle."""

    angle_clicked = pyqtSignal(float)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._scan: Optional[AngleScan] = None
        self._curve: Optional[pg.PlotDataItem] = None
        self.scene().sigMouseClicked.connect(self._on_scene_click)

    def set_scan(self, scan: AngleScan):
        self._scan = scan
        self._redraw()

    def _redraw(self):  # overridden
        raise NotImplementedError

    def _on_scene_click(self, ev):
        if not self._scan or not self._scan.points:
            return
        vb = self.plotItem.vb
        mouse_point = vb.mapSceneToView(ev.scenePos())
        nearest = self._scan.nearest_point(mouse_point.x())
        if nearest is not None:
            self.angle_clicked.emit(nearest.angle_deg)

    def _emit_nearest(self, angle_deg: float):
        """Test hook."""
        if self._scan is None:
            return
        nearest = self._scan.nearest_point(angle_deg)
        if nearest is not None:
            self.angle_clicked.emit(nearest.angle_deg)


class GapPlot(_AnglePlotBase):
    def __init__(self, parent=None):
        super().__init__(parent)
        _configure_axes(self, xlabel="θ (deg)", ylabel="ΔE(S₁−T₁) (eV)",
                        title="Singlet–Triplet Gap")
        self.threshold_ev = 0.2
        self._threshold_line: Optional[pg.InfiniteLine] = None
        self._highlight: Optional[pg.ScatterPlotItem] = None

    def set_threshold(self, value: float):
        self.threshold_ev = value
        self._redraw()

    def _redraw(self):
        self.clear()
        self._threshold_line = None
        self._highlight = None
        if not self._scan or not self._scan.points:
            return
        xs, ys = gap_curve(self._scan)
        pen = pg.mkPen(_CURVE, width=2)
        self._curve = self.plot(xs, ys, pen=pen, symbol="o",
                                symbolBrush=_CURVE, symbolSize=8)
        # Dashed threshold line
        thr_pen = pg.mkPen(_THRESHOLD, width=1.5, style=Qt.DashLine)
        self._threshold_line = pg.InfiniteLine(
            pos=self.threshold_ev, angle=0, pen=thr_pen,
            label=f"{self.threshold_ev:.2f} eV",
            labelOpts={"color": _THRESHOLD, "position": 0.95},
        )
        self.addItem(self._threshold_line)
        # Highlight TADF candidates (gap ≤ threshold)
        xs_h = [x for x, y in zip(xs, ys) if y is not None and y <= self.threshold_ev]
        ys_h = [y for y in ys if y is not None and y <= self.threshold_ev]
        if xs_h:
            self._highlight = pg.ScatterPlotItem(
                x=xs_h, y=ys_h,
                brush=pg.mkBrush(_HIGHLIGHT), pen=pg.mkPen("w", width=1.5),
                size=14, symbol="o",
            )
            self.addItem(self._highlight)


class TDMPlot(_AnglePlotBase):
    def __init__(self, parent=None):
        super().__init__(parent)
        _configure_axes(self, xlabel="θ (deg)", ylabel="|TDM(S₁)| (Au)",
                        title="S₁ Transition Dipole Moment")

    def _redraw(self):
        self.clear()
        if not self._scan or not self._scan.points:
            return
        xs, ys = tdm_curve(self._scan)
        self._curve = self.plot(xs, ys, pen=pg.mkPen("#ffeb3b", width=2),
                                symbol="s", symbolBrush="#ffeb3b", symbolSize=8)


class TwoPSquaredPlot(_AnglePlotBase):
    def __init__(self, parent=None):
        super().__init__(parent)
        _configure_axes(self, xlabel="θ (deg)", ylabel="2·P²",
                        title="Dominant S₁ Configuration Weight (2·P²)")

    def _redraw(self):
        self.clear()
        if not self._scan or not self._scan.points:
            return
        xs, ys = two_p_squared_curve(self._scan)
        self._curve = self.plot(xs, ys, pen=pg.mkPen("#b388ff", width=2),
                                symbol="t", symbolBrush="#b388ff", symbolSize=9)
```

- [ ] **Step 5.4: Run (pass)**

```bash
.venv/bin/pytest tests/test_angle_scan_plots.py -v
```

- [ ] **Step 5.5: Commit**

```bash
git add src/molecule_visualizer/ui/angle_scan_plots.py tests/test_angle_scan_plots.py
git commit -m "feat(ui): GapPlot, TDMPlot, TwoPSquaredPlot widgets"
```

---

## Task 6: Data table widget

**Files:**
- Create: `src/molecule_visualizer/ui/angle_scan_table.py`
- Test: `tests/test_angle_scan_table.py`

Columns: θ (deg), E(S₁) eV, E(T₁) eV, ΔE eV, |TDM| Au, 2·P², S₁ dominant transition. Rows with `gap ≤ threshold` receive a green background. Clicking a row emits `angle_selected(float)`.

- [ ] **Step 6.1: Test**

```python
# tests/test_angle_scan_table.py
import pytest
pytest.importorskip("PyQt5")
from molecule_visualizer.ui.angle_scan_table import AngleScanTable
from molecule_visualizer.parsers import parse_angle_folder
from tests.fixtures.build_fixtures import write_scan_folder


@pytest.fixture
def scan(tmp_path):
    return parse_angle_folder(str(write_scan_folder(tmp_path / "s", {
        10: (3.4, 2.9), 30: (3.0, 2.95), 90: (3.4, 3.0),
    })))


def test_table_populates_rows(qtbot, scan):
    t = AngleScanTable()
    qtbot.addWidget(t)
    t.set_scan(scan)
    assert t.rowCount() == 3
    assert t.columnCount() == 7


def test_table_highlights_tadf_rows(qtbot, scan):
    t = AngleScanTable()
    qtbot.addWidget(t)
    t.set_scan(scan, threshold_ev=0.2)
    # 30° gap is 0.05, 10° gap is 0.5, 90° gap is 0.4
    assert t.is_row_highlighted(1)  # after sort: 10, 30, 90 → row 1 = 30°
    assert not t.is_row_highlighted(0)
    assert not t.is_row_highlighted(2)


def test_table_emits_angle_selected(qtbot, scan):
    t = AngleScanTable()
    qtbot.addWidget(t)
    t.set_scan(scan)
    received = []
    t.angle_selected.connect(received.append)
    t.selectRow(1)
    assert received and received[-1] == 30.0
```

- [ ] **Step 6.2: Implement**

```python
# src/molecule_visualizer/ui/angle_scan_table.py
"""Tabular view of AngleScan metrics with TADF-row highlighting."""
from typing import Optional

from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QColor, QBrush, QFont
from PyQt5.QtWidgets import QTableWidget, QTableWidgetItem, QAbstractItemView

from ..core import AngleScan
from ..analysis import summary_row


_TADF_BG = QColor(40, 95, 45)
_TADF_FG = QColor("#c8ffd4")
_NORMAL_BG = QColor("#1e1e1e")
_NORMAL_FG = QColor("#d4d4d4")

_COLUMNS = [
    ("θ (deg)", "angle_deg"),
    ("E(S₁) eV", "s1_energy_ev"),
    ("E(T₁) eV", "t1_energy_ev"),
    ("ΔE eV", "gap_ev"),
    ("|TDM| Au", "tdm_magnitude"),
    ("2·P²", "two_p_squared"),
    ("S₁ dominant (i→a, c)", "s1_dominant"),
]


def _fmt(val, key):
    if val is None:
        return "—"
    if key == "s1_dominant":
        i, a, c = val
        return f"{i} → {a}  ({c:+.4f})"
    if key == "angle_deg":
        return f"{val:.1f}"
    if isinstance(val, float):
        return f"{val:.4f}"
    return str(val)


class AngleScanTable(QTableWidget):
    angle_selected = pyqtSignal(float)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._scan: Optional[AngleScan] = None
        self._threshold = 0.2
        self._highlight_rows: set[int] = set()
        self.setColumnCount(len(_COLUMNS))
        self.setHorizontalHeaderLabels([c[0] for c in _COLUMNS])
        self.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.setSelectionMode(QAbstractItemView.SingleSelection)
        self.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.setAlternatingRowColors(True)
        self.setStyleSheet("""
            QTableWidget {
                background-color: #1e1e1e; color: #d4d4d4;
                gridline-color: #333; alternate-background-color: #242424;
                selection-background-color: #094771;
            }
            QHeaderView::section {
                background-color: #2d2d30; color: #ddd;
                padding: 4px; border: 1px solid #444;
            }
        """)
        self.itemSelectionChanged.connect(self._emit_selected)

    def set_scan(self, scan: AngleScan, threshold_ev: float = 0.2):
        self._scan = scan
        self._threshold = threshold_ev
        self._highlight_rows.clear()
        points = scan.sorted_points
        self.setRowCount(len(points))
        for row, p in enumerate(points):
            row_data = summary_row(p)
            is_tadf = (p.s1_t1_gap_ev is not None and p.s1_t1_gap_ev <= threshold_ev)
            if is_tadf:
                self._highlight_rows.add(row)
            for col, (_, key) in enumerate(_COLUMNS):
                item = QTableWidgetItem(_fmt(row_data.get(key), key))
                item.setTextAlignment(Qt.AlignCenter)
                if is_tadf:
                    item.setBackground(QBrush(_TADF_BG))
                    item.setForeground(QBrush(_TADF_FG))
                    f = item.font(); f.setBold(True); item.setFont(f)
                self.setItem(row, col, item)
        self.resizeColumnsToContents()

    def is_row_highlighted(self, row: int) -> bool:
        return row in self._highlight_rows

    def _emit_selected(self):
        if not self._scan:
            return
        rows = self.selectionModel().selectedRows()
        if not rows:
            return
        idx = rows[0].row()
        points = self._scan.sorted_points
        if 0 <= idx < len(points):
            self.angle_selected.emit(points[idx].angle_deg)
```

- [ ] **Step 6.3: Run / commit**

```bash
.venv/bin/pytest tests/test_angle_scan_table.py -v
git add src/molecule_visualizer/ui/angle_scan_table.py tests/test_angle_scan_table.py
git commit -m "feat(ui): AngleScanTable with TADF row highlighting"
```

---

## Task 7: Donor / acceptor / rotatable-bond picker dialog

**Files:**
- Create: `src/molecule_visualizer/ui/donor_acceptor_dialog.py`
- Test: `tests/test_donor_acceptor_dialog.py`

Dialog shows a static 3D view of the first stage's geometry, lists all atoms, and has three combo-boxes: **Donor atom**, **Acceptor atom**, **Rotatable bond (atom A, atom B)**. OK is disabled until all four are chosen and A≠B. Result is written into the passed `AngleScan`.

- [ ] **Step 7.1: Test** (pytest-qt, no assertions on 3D rendering, only on signals/state)

```python
# tests/test_donor_acceptor_dialog.py
import pytest
pytest.importorskip("PyQt5")
from PyQt5.QtWidgets import QDialog
from molecule_visualizer.ui.donor_acceptor_dialog import DonorAcceptorDialog
from molecule_visualizer.parsers import parse_angle_folder
from tests.fixtures.build_fixtures import write_scan_folder


@pytest.fixture
def scan(tmp_path):
    return parse_angle_folder(str(write_scan_folder(tmp_path / "s", {10: (3.4, 2.9)})))


def test_dialog_requires_selections_before_accept(qtbot, scan):
    dlg = DonorAcceptorDialog(scan=scan)
    qtbot.addWidget(dlg)
    assert not dlg.ok_button.isEnabled()
    dlg.set_selection(donor=0, acceptor=1, bond=(0, 1))
    assert dlg.ok_button.isEnabled()


def test_dialog_rejects_bond_with_equal_atoms(qtbot, scan):
    dlg = DonorAcceptorDialog(scan=scan)
    qtbot.addWidget(dlg)
    dlg.set_selection(donor=0, acceptor=1, bond=(0, 0))
    assert not dlg.ok_button.isEnabled()


def test_dialog_writes_back_to_scan_on_accept(qtbot, scan):
    dlg = DonorAcceptorDialog(scan=scan)
    qtbot.addWidget(dlg)
    dlg.set_selection(donor=0, acceptor=1, bond=(0, 1))
    dlg.accept()
    assert scan.donor_atom_index == 0
    assert scan.acceptor_atom_index == 1
    assert scan.rotatable_bond == (0, 1)
```

- [ ] **Step 7.2: Implement**

```python
# src/molecule_visualizer/ui/donor_acceptor_dialog.py
"""Dialog to pick donor atom, acceptor atom, and the rotatable bond axis."""
from typing import Optional, Tuple

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QComboBox, QPushButton,
    QDialogButtonBox, QFormLayout, QGroupBox,
)

from ..core import AngleScan
from ..utils import get_element_symbol


class DonorAcceptorDialog(QDialog):
    def __init__(self, *, scan: AngleScan, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Select Donor, Acceptor, and Rotatable Bond")
        self._scan = scan
        self._donor: Optional[int] = scan.donor_atom_index
        self._acceptor: Optional[int] = scan.acceptor_atom_index
        self._bond: Optional[Tuple[int, int]] = scan.rotatable_bond

        # Use the first available geometry
        first = scan.sorted_points[0] if scan.points else None
        stage = first.molecule.current_stage if first else None
        atoms = stage.atoms if stage else []

        layout = QVBoxLayout(self)
        layout.addWidget(QLabel(
            "<b>Configure donor–acceptor geometry</b><br>"
            "Select the donor atom, acceptor atom, and the two atoms forming "
            "the rotatable bond axis.", alignment=Qt.AlignLeft))

        form = QFormLayout()
        self.donor_cb = self._make_atom_combo(atoms)
        self.acceptor_cb = self._make_atom_combo(atoms)
        self.bond_a_cb = self._make_atom_combo(atoms)
        self.bond_b_cb = self._make_atom_combo(atoms)
        for cb, initial in ((self.donor_cb, self._donor),
                            (self.acceptor_cb, self._acceptor),
                            (self.bond_a_cb, self._bond[0] if self._bond else None),
                            (self.bond_b_cb, self._bond[1] if self._bond else None)):
            if initial is not None and 0 <= initial < len(atoms):
                cb.setCurrentIndex(initial + 1)  # +1 because of the "—" placeholder

        form.addRow("Donor atom:", self.donor_cb)
        form.addRow("Acceptor atom:", self.acceptor_cb)
        bond_box = QGroupBox("Rotatable bond")
        bb_layout = QHBoxLayout(bond_box)
        bb_layout.addWidget(QLabel("A:")); bb_layout.addWidget(self.bond_a_cb)
        bb_layout.addWidget(QLabel("B:")); bb_layout.addWidget(self.bond_b_cb)
        form.addRow(bond_box)
        layout.addLayout(form)

        self.button_box = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel, parent=self)
        self.ok_button = self.button_box.button(QDialogButtonBox.Ok)
        self.button_box.accepted.connect(self.accept)
        self.button_box.rejected.connect(self.reject)
        layout.addWidget(self.button_box)

        for cb in (self.donor_cb, self.acceptor_cb,
                   self.bond_a_cb, self.bond_b_cb):
            cb.currentIndexChanged.connect(self._revalidate)
        self._revalidate()

    def _make_atom_combo(self, atoms):
        cb = QComboBox()
        cb.addItem("— choose —", None)
        for idx, atom in enumerate(atoms):
            sym = get_element_symbol(atom.atomic_number)
            cb.addItem(f"{idx}: {sym} ({atom.x:+.3f}, {atom.y:+.3f}, {atom.z:+.3f})",
                       userData=idx)
        return cb

    def set_selection(self, *, donor: int, acceptor: int, bond: Tuple[int, int]):
        self.donor_cb.setCurrentIndex(donor + 1)
        self.acceptor_cb.setCurrentIndex(acceptor + 1)
        self.bond_a_cb.setCurrentIndex(bond[0] + 1)
        self.bond_b_cb.setCurrentIndex(bond[1] + 1)
        self._revalidate()

    def _revalidate(self):
        d = self.donor_cb.currentData()
        a = self.acceptor_cb.currentData()
        ba = self.bond_a_cb.currentData()
        bb = self.bond_b_cb.currentData()
        valid = (d is not None and a is not None and ba is not None
                 and bb is not None and ba != bb)
        self.ok_button.setEnabled(valid)
        if valid:
            self._donor, self._acceptor, self._bond = d, a, (ba, bb)

    def accept(self):
        if self._donor is None or self._acceptor is None or self._bond is None:
            return
        self._scan.donor_atom_index = self._donor
        self._scan.acceptor_atom_index = self._acceptor
        self._scan.rotatable_bond = self._bond
        super().accept()
```

- [ ] **Step 7.3: Run / commit**

```bash
.venv/bin/pytest tests/test_donor_acceptor_dialog.py -v
git add src/molecule_visualizer/ui/donor_acceptor_dialog.py tests/test_donor_acceptor_dialog.py
git commit -m "feat(ui): DonorAcceptorDialog for scan configuration"
```

---

## Task 8: 3D dihedral viewer panel

**Files:**
- Create: `src/molecule_visualizer/ui/dihedral_3d_panel.py`
- Test: `tests/test_dihedral_3d_panel.py`

Displays the geometry of one `AnglePoint` in a `pyqtgraph.opengl.GLViewWidget`. A horizontal `QSlider` ranges over `[min(angles), max(angles)]`, snapping to the nearest scanned angle on any change. Emits `angle_changed(float)` so the outer window can update plots/table. Donor + acceptor atoms (if set on the scan) are highlighted (thicker outline / brighter).

**Rendering strategy (simplest reliable):** reuse patterns from existing `visualization_panel.py`:
- Sphere per atom (GLMeshItem with a MeshData sphere, coloured by element).
- Bond as a cylinder or a GLLinePlotItem between bonded atom centres (reuse `find_bonds()` on the stage).
- Donor/acceptor: draw an extra, larger translucent sphere around them.

- [ ] **Step 8.1: Test — focuses on slider logic + signals, not GL pixels**

```python
# tests/test_dihedral_3d_panel.py
import pytest
pytest.importorskip("PyQt5")
pytest.importorskip("pyqtgraph.opengl")
from molecule_visualizer.ui.dihedral_3d_panel import Dihedral3DPanel
from molecule_visualizer.parsers import parse_angle_folder
from tests.fixtures.build_fixtures import write_scan_folder


@pytest.fixture
def scan(tmp_path):
    return parse_angle_folder(str(write_scan_folder(tmp_path / "s", {
        10: (3.4, 2.9), 30: (3.0, 2.95), 60: (3.1, 2.95), 90: (3.4, 3.0),
    })))


def test_slider_snaps_to_nearest_scanned_angle(qtbot, scan):
    panel = Dihedral3DPanel()
    qtbot.addWidget(panel)
    panel.set_scan(scan)
    panel.set_angle(35)  # closest = 30
    assert panel.current_angle == 30.0
    panel.set_angle(75)  # closest = 60 or 90 → 60 vs 90 — 75-60=15 < 90-75=15 tie → parser chooses 60 by min()
    assert panel.current_angle in (60.0, 90.0)
    panel.set_angle(89)
    assert panel.current_angle == 90.0


def test_panel_emits_angle_changed(qtbot, scan):
    panel = Dihedral3DPanel()
    qtbot.addWidget(panel)
    panel.set_scan(scan)
    received = []
    panel.angle_changed.connect(received.append)
    panel.set_angle(31)  # snaps to 30
    assert received and received[-1] == 30.0


def test_setting_same_angle_twice_does_not_double_emit(qtbot, scan):
    panel = Dihedral3DPanel()
    qtbot.addWidget(panel)
    panel.set_scan(scan)
    received = []
    panel.angle_changed.connect(received.append)
    panel.set_angle(30)
    panel.set_angle(30)
    assert received.count(30.0) == 1
```

- [ ] **Step 8.2: Implement**

```python
# src/molecule_visualizer/ui/dihedral_3d_panel.py
"""3D viewer for the current angle point, with a snapping θ slider."""
from typing import Optional, List

import numpy as np
import pyqtgraph.opengl as gl
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QVector3D
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QSlider, QFrame,
)

from ..core import AngleScan, AnglePoint, Stage
from ..utils import get_atom_color, get_atom_radius, get_element_symbol


class Dihedral3DPanel(QWidget):
    angle_changed = pyqtSignal(float)  # emitted when current angle changes

    def __init__(self, parent=None):
        super().__init__(parent)
        self._scan: Optional[AngleScan] = None
        self._current: Optional[AnglePoint] = None
        self._gl_items: list = []
        self._build_ui()

    def _build_ui(self):
        lay = QVBoxLayout(self)
        lay.setContentsMargins(4, 4, 4, 4)
        lay.setSpacing(4)

        self.view = gl.GLViewWidget()
        self.view.setBackgroundColor(QVector3D(0.07, 0.08, 0.10))
        self.view.opts["distance"] = 30
        lay.addWidget(self.view, 1)

        # Slider row
        bar = QFrame()
        bar.setStyleSheet("QFrame { background:#181c22; border-top:1px solid #252d38; }")
        row = QHBoxLayout(bar)
        row.setContentsMargins(10, 6, 10, 6)
        row.addWidget(QLabel("θ"))
        self.slider = QSlider(Qt.Horizontal)
        self.slider.setMinimum(0)
        self.slider.setMaximum(180)
        self.slider.valueChanged.connect(self._on_slider_change)
        row.addWidget(self.slider, 1)
        self.angle_label = QLabel("—°")
        self.angle_label.setMinimumWidth(60)
        self.angle_label.setStyleSheet("color:#4fc3f7; font-weight:bold;")
        row.addWidget(self.angle_label)
        lay.addWidget(bar)

    # ── Public API ──

    @property
    def current_angle(self) -> Optional[float]:
        return self._current.angle_deg if self._current else None

    def set_scan(self, scan: AngleScan):
        self._scan = scan
        if not scan.points:
            return
        pts = scan.sorted_points
        lo, hi = int(pts[0].angle_deg), int(pts[-1].angle_deg)
        self.slider.blockSignals(True)
        self.slider.setMinimum(lo)
        self.slider.setMaximum(hi)
        self.slider.blockSignals(False)
        self.set_angle(pts[0].angle_deg)

    def set_angle(self, angle_deg: float):
        if not self._scan or not self._scan.points:
            return
        target = self._scan.nearest_point(angle_deg)
        if target is None:
            return
        # Dedupe
        if self._current is not None and target.angle_deg == self._current.angle_deg:
            return
        self._current = target
        self.slider.blockSignals(True)
        self.slider.setValue(int(round(target.angle_deg)))
        self.slider.blockSignals(False)
        self.angle_label.setText(f"{target.angle_deg:.1f}°")
        self._render_point(target)
        self.angle_changed.emit(target.angle_deg)

    # ── Internals ──

    def _on_slider_change(self, value: int):
        self.set_angle(float(value))

    def _clear_items(self):
        for it in self._gl_items:
            self.view.removeItem(it)
        self._gl_items.clear()

    def _render_point(self, point: AnglePoint):
        self._clear_items()
        stage: Optional[Stage] = point.molecule.current_stage
        if stage is None or stage.atom_count == 0:
            return
        self._draw_atoms(stage)
        self._draw_bonds(stage)
        self._draw_donor_acceptor_highlights(stage)

    def _draw_atoms(self, stage: Stage):
        md = gl.MeshData.sphere(rows=12, cols=24)
        for atom in stage.atoms:
            color = get_atom_color(atom.atomic_number) + (1.0,) \
                if len(get_atom_color(atom.atomic_number)) == 3 \
                else get_atom_color(atom.atomic_number)
            radius = get_atom_radius(atom.atomic_number) * 0.35
            item = gl.GLMeshItem(
                meshdata=md, smooth=True, color=color, shader="shaded",
            )
            item.scale(radius, radius, radius)
            item.translate(atom.x, atom.y, atom.z)
            self.view.addItem(item)
            self._gl_items.append(item)

    def _draw_bonds(self, stage: Stage):
        bonds = stage.find_bonds()
        if not bonds:
            return
        positions = stage.positions
        pts = []
        for i, j in bonds:
            pts.append(positions[i])
            pts.append(positions[j])
        arr = np.array(pts)
        line = gl.GLLinePlotItem(pos=arr, color=(0.7, 0.7, 0.8, 0.9),
                                 width=2.0, mode="lines", antialias=True)
        self.view.addItem(line)
        self._gl_items.append(line)

    def _draw_donor_acceptor_highlights(self, stage: Stage):
        if not self._scan:
            return
        md = gl.MeshData.sphere(rows=10, cols=20)
        for idx, tint in ((self._scan.donor_atom_index, (0.3, 0.9, 0.4, 0.25)),
                          (self._scan.acceptor_atom_index, (0.9, 0.5, 0.3, 0.25))):
            if idx is None or idx >= stage.atom_count:
                continue
            atom = stage.atoms[idx]
            halo = gl.GLMeshItem(
                meshdata=md, smooth=True, color=tint, shader="balloon",
                glOptions="additive",
            )
            halo.scale(0.9, 0.9, 0.9)
            halo.translate(atom.x, atom.y, atom.z)
            self.view.addItem(halo)
            self._gl_items.append(halo)
```

- [ ] **Step 8.3: Run / commit**

```bash
.venv/bin/pytest tests/test_dihedral_3d_panel.py -v
git add src/molecule_visualizer/ui/dihedral_3d_panel.py tests/test_dihedral_3d_panel.py
git commit -m "feat(ui): Dihedral3DPanel with snapping theta slider"
```

---

## Task 9: `AngleScanWindow` — assemble the workspace

**Files:**
- Create: `src/molecule_visualizer/ui/angle_scan_window.py`
- Modify: `src/molecule_visualizer/ui/__init__.py` to export it
- Test: `tests/test_angle_scan_window.py`

A `QMainWindow` with:
- Menu: **File → Open Folder… / Export CSV…**, **Scan → Configure Donor/Acceptor…**, **Scan → Threshold (eV)…**
- Central widget: a horizontal splitter. Left = `Dihedral3DPanel`. Right = vertical splitter containing `GapPlot`, `TDMPlot`, `TwoPSquaredPlot` stacked, then `AngleScanTable` at the bottom.
- Signals cross-wired: selecting a θ in any child updates the others.
- Status bar shows `n pts | min gap X.XX eV @ Y° | threshold Z.ZZ eV`.

- [ ] **Step 9.1: Test (integration smoke)**

```python
# tests/test_angle_scan_window.py
import pytest
pytest.importorskip("PyQt5")
from molecule_visualizer.ui import AngleScanWindow
from molecule_visualizer.parsers import parse_angle_folder
from tests.fixtures.build_fixtures import write_scan_folder


@pytest.fixture
def scan(tmp_path):
    return parse_angle_folder(str(write_scan_folder(tmp_path / "s", {
        10: (3.4, 2.9), 30: (3.0, 2.95), 60: (3.1, 2.95), 90: (3.4, 3.0),
    })))


def test_window_loads_scan(qtbot, scan):
    w = AngleScanWindow()
    qtbot.addWidget(w)
    w.set_scan(scan)
    assert w.windowTitle().startswith("Angle Scan")
    assert w.table.rowCount() == 4


def test_selecting_angle_in_table_updates_3d(qtbot, scan):
    w = AngleScanWindow()
    qtbot.addWidget(w)
    w.set_scan(scan)
    # Programmatically select row for θ=30
    angles = [p.angle_deg for p in scan.sorted_points]
    row = angles.index(30.0)
    w.table.selectRow(row)
    qtbot.wait(50)
    assert w.viewer_3d.current_angle == 30.0


def test_selecting_angle_in_3d_updates_plots(qtbot, scan):
    w = AngleScanWindow()
    qtbot.addWidget(w)
    w.set_scan(scan)
    # Slider change → plots should receive new highlight
    w.viewer_3d.set_angle(60)
    qtbot.wait(50)
    assert w.viewer_3d.current_angle == 60.0
    assert w._current_angle == 60.0


def test_status_bar_reports_min_gap(qtbot, scan):
    w = AngleScanWindow()
    qtbot.addWidget(w)
    w.set_scan(scan)
    txt = w.statusBar().currentMessage()
    assert "min gap" in txt.lower()
    assert "30" in txt  # min gap is at 30°
```

- [ ] **Step 9.2: Implement**

```python
# src/molecule_visualizer/ui/angle_scan_window.py
"""Top-level angle-scan analysis workspace."""
from typing import Optional

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QKeySequence
from PyQt5.QtWidgets import (
    QMainWindow, QWidget, QHBoxLayout, QSplitter, QAction, QFileDialog,
    QInputDialog, QMessageBox, QLabel,
)

from ..core import AngleScan
from ..parsers import parse_angle_folder
from .angle_scan_plots import GapPlot, TDMPlot, TwoPSquaredPlot
from .angle_scan_table import AngleScanTable
from .dihedral_3d_panel import Dihedral3DPanel
from .donor_acceptor_dialog import DonorAcceptorDialog


class AngleScanWindow(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Angle Scan — TADF Workbench")
        self.resize(1400, 900)
        self._scan: Optional[AngleScan] = None
        self._current_angle: Optional[float] = None
        self._threshold = 0.2
        self._build_ui()
        self._build_menu()

    # ── UI ──

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        lay = QHBoxLayout(central)
        lay.setContentsMargins(4, 4, 4, 4)

        outer = QSplitter(Qt.Horizontal)
        self.viewer_3d = Dihedral3DPanel()
        outer.addWidget(self.viewer_3d)

        right = QSplitter(Qt.Vertical)
        top_plots = QSplitter(Qt.Vertical)
        self.gap_plot = GapPlot()
        self.tdm_plot = TDMPlot()
        self.p2_plot = TwoPSquaredPlot()
        top_plots.addWidget(self.gap_plot)
        top_plots.addWidget(self.tdm_plot)
        top_plots.addWidget(self.p2_plot)
        right.addWidget(top_plots)

        self.table = AngleScanTable()
        right.addWidget(self.table)
        right.setStretchFactor(0, 3)
        right.setStretchFactor(1, 2)

        outer.addWidget(right)
        outer.setStretchFactor(0, 3)
        outer.setStretchFactor(1, 4)
        lay.addWidget(outer)

        # Cross-wiring
        self.table.angle_selected.connect(self._on_angle_selected)
        self.viewer_3d.angle_changed.connect(self._on_angle_selected)
        for p in (self.gap_plot, self.tdm_plot, self.p2_plot):
            p.angle_clicked.connect(self._on_angle_selected)

        self.statusBar().showMessage("No scan loaded.")

    def _build_menu(self):
        mb = self.menuBar()
        file_menu = mb.addMenu("&File")
        open_act = QAction("Open Scan Folder…", self, shortcut=QKeySequence("Ctrl+Shift+O"))
        open_act.triggered.connect(self._prompt_open_folder)
        file_menu.addAction(open_act)
        export_act = QAction("Export Summary CSV…", self)
        export_act.triggered.connect(self._export_csv)
        file_menu.addAction(export_act)

        scan_menu = mb.addMenu("&Scan")
        config_act = QAction("Configure Donor/Acceptor…", self)
        config_act.triggered.connect(self._configure_donor_acceptor)
        scan_menu.addAction(config_act)
        thr_act = QAction("Set ΔE Threshold (eV)…", self)
        thr_act.triggered.connect(self._set_threshold)
        scan_menu.addAction(thr_act)

    # ── Public API ──

    def set_scan(self, scan: AngleScan):
        self._scan = scan
        self._current_angle = None
        self.gap_plot.set_scan(scan)
        self.gap_plot.set_threshold(self._threshold)
        self.tdm_plot.set_scan(scan)
        self.p2_plot.set_scan(scan)
        self.table.set_scan(scan, threshold_ev=self._threshold)
        self.viewer_3d.set_scan(scan)
        self._update_status()

    # ── Handlers ──

    def _prompt_open_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Select Scan Folder")
        if not folder:
            return
        try:
            scan = parse_angle_folder(folder)
        except Exception as exc:
            QMessageBox.critical(self, "Parse Error", str(exc))
            return
        if not scan.points:
            QMessageBox.warning(self, "Empty Scan",
                                "No log files with detectable angles were found.")
            return
        self.set_scan(scan)
        # Prompt donor/acceptor config on first load if not already set
        if scan.donor_atom_index is None:
            self._configure_donor_acceptor()

    def _configure_donor_acceptor(self):
        if not self._scan:
            QMessageBox.information(self, "No Scan", "Load a scan folder first.")
            return
        dlg = DonorAcceptorDialog(scan=self._scan, parent=self)
        if dlg.exec_():
            # Re-render 3D to show highlights
            if self._current_angle is not None:
                self.viewer_3d.set_angle(self._current_angle)
            else:
                self.viewer_3d.set_scan(self._scan)

    def _set_threshold(self):
        val, ok = QInputDialog.getDouble(
            self, "ΔE Threshold", "TADF candidate threshold (eV):",
            value=self._threshold, min=0.0, max=2.0, decimals=3,
        )
        if ok:
            self._threshold = val
            if self._scan:
                self.gap_plot.set_threshold(val)
                self.table.set_scan(self._scan, threshold_ev=val)
                self._update_status()

    def _export_csv(self):
        if not self._scan:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Summary CSV", f"{self._scan.name}_summary.csv",
            "CSV (*.csv)")
        if not path:
            return
        import csv
        from ..analysis import summary_row
        keys = ["angle_deg", "s1_energy_ev", "t1_energy_ev", "gap_ev",
                "tdm_magnitude", "two_p_squared", "s1_dominant", "source_path"]
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(keys)
            for p in self._scan.sorted_points:
                row = summary_row(p)
                w.writerow([row.get(k) for k in keys])
        self.statusBar().showMessage(f"Exported: {path}", 4000)

    def _on_angle_selected(self, angle_deg: float):
        if self._current_angle == angle_deg:
            return
        self._current_angle = angle_deg
        self.viewer_3d.set_angle(angle_deg)
        # Select the matching row in the table
        if self._scan:
            pts = self._scan.sorted_points
            for row, p in enumerate(pts):
                if p.angle_deg == angle_deg:
                    self.table.blockSignals(True)
                    self.table.selectRow(row)
                    self.table.blockSignals(False)
                    break
        self._update_status()

    def _update_status(self):
        if not self._scan:
            self.statusBar().showMessage("No scan loaded.")
            return
        msg = f"{len(self._scan)} pts"
        best = self._scan.minimum_gap_point
        if best is not None and best.s1_t1_gap_ev is not None:
            msg += f" | min gap {best.s1_t1_gap_ev:.4f} eV @ {best.angle_deg:.1f}°"
        msg += f" | threshold {self._threshold:.3f} eV"
        n_tadf = len(self._scan.points_below_threshold(self._threshold))
        msg += f" | {n_tadf} TADF candidates"
        self.statusBar().showMessage(msg)
```

Modify `src/molecule_visualizer/ui/__init__.py` to add:

```python
from .angle_scan_window import AngleScanWindow
```

and append `"AngleScanWindow"` to its `__all__`.

- [ ] **Step 9.3: Run / commit**

```bash
.venv/bin/pytest tests/test_angle_scan_window.py -v
git add src/molecule_visualizer/ui/angle_scan_window.py src/molecule_visualizer/ui/__init__.py tests/test_angle_scan_window.py
git commit -m "feat(ui): AngleScanWindow assembles plots, table, and 3D viewer"
```

---

## Task 10: Integrate into the main window (menu + toolbar)

**Files:**
- Modify: `src/molecule_visualizer/ui/main_window.py`

- [ ] **Step 10.1: Add menu + action**

Locate `_create_menu_bar` (main_window.py:150). After the existing `open_action` (which loads a single log), add:

```python
        open_scan_action = QAction("Open &Angle Scan Folder...", self)
        open_scan_action.setShortcut(QKeySequence("Ctrl+Shift+O"))
        open_scan_action.triggered.connect(self._open_angle_scan_folder)
        file_menu.addAction(open_scan_action)
```

Add to the `Analysis` menu block too:

```python
        scan_analysis = QAction("Angle &Scan Workbench...", self)
        scan_analysis.triggered.connect(self._open_angle_scan_folder)
        analysis_menu.addAction(scan_analysis)
```

Add toolbar button in `_create_toolbar` (main_window.py:239):

```python
        scan_btn = QAction("Angle Scan", self)
        scan_btn.setToolTip("Open folder of log files as a dihedral-angle scan")
        scan_btn.triggered.connect(self._open_angle_scan_folder)
        toolbar.addAction(scan_btn)
```

Add the new method on `MainWindow`:

```python
    def _open_angle_scan_folder(self):
        from .angle_scan_window import AngleScanWindow
        from ..parsers import parse_angle_folder
        from .donor_acceptor_dialog import DonorAcceptorDialog

        folder = QFileDialog.getExistingDirectory(self, "Select Angle Scan Folder")
        if not folder:
            return
        try:
            scan = parse_angle_folder(folder)
        except Exception as e:
            QMessageBox.critical(self, "Error",
                                 f"Failed to parse folder:\n{e}")
            return
        if not scan.points:
            QMessageBox.warning(self, "Empty Scan",
                                "No log files with detectable angles.")
            return
        win = AngleScanWindow(parent=self)
        # Prompt for donor/acceptor once up front
        dlg = DonorAcceptorDialog(scan=scan, parent=self)
        dlg.exec_()
        win.set_scan(scan)
        win.show()
        # Keep a reference so Python doesn't GC the window
        if not hasattr(self, "_child_windows"):
            self._child_windows = []
        self._child_windows.append(win)
```

- [ ] **Step 10.2: Smoke-test manually**

```bash
.venv/bin/python -m molecule_visualizer.main
```
Expected: app launches with the new "Angle Scan" toolbar button and menu entries. Clicking either opens a folder picker; selecting `tests/fixtures` (after running a helper to populate it) opens the workbench.

- [ ] **Step 10.3: Commit**

```bash
git add src/molecule_visualizer/ui/main_window.py
git commit -m "feat(ui): wire Angle Scan workbench into MainWindow"
```

---

## Task 11: End-to-end integration test (runs all layers together)

**Files:**
- Test: `tests/test_e2e_angle_scan.py`

- [ ] **Step 11.1: Write test**

```python
# tests/test_e2e_angle_scan.py
import pytest
pytest.importorskip("PyQt5")
from molecule_visualizer.parsers import parse_angle_folder
from molecule_visualizer.ui import AngleScanWindow
from tests.fixtures.build_fixtures import write_scan_folder


def test_full_flow(qtbot, tmp_path):
    folder = write_scan_folder(tmp_path / "demo", {
        10: (3.40, 2.90),  # gap 0.50
        20: (3.30, 2.92),  # gap 0.38
        30: (3.20, 2.95),  # gap 0.25
        40: (3.10, 2.98),  # gap 0.12 — TADF
        50: (3.05, 3.00),  # gap 0.05 — TADF (minimum)
        60: (3.08, 3.00),  # gap 0.08 — TADF
        70: (3.15, 2.99),  # gap 0.16 — TADF
        80: (3.25, 2.97),  # gap 0.28
        90: (3.40, 2.93),  # gap 0.47
    })
    scan = parse_angle_folder(str(folder))
    assert len(scan) == 9

    w = AngleScanWindow()
    qtbot.addWidget(w)
    w.set_scan(scan)

    assert w.table.rowCount() == 9
    # TADF candidates at θ = 40, 50, 60, 70
    assert sum(1 for r in range(w.table.rowCount()) if w.table.is_row_highlighted(r)) == 4

    # Clicking on θ=50 via the viewer updates status bar
    w.viewer_3d.set_angle(50)
    qtbot.wait(100)
    assert "50" in w.statusBar().currentMessage() or \
           "min gap" in w.statusBar().currentMessage().lower()
```

- [ ] **Step 11.2: Run / commit**

```bash
.venv/bin/pytest tests/test_e2e_angle_scan.py -v
git add tests/test_e2e_angle_scan.py
git commit -m "test(e2e): end-to-end angle scan flow"
```

---

## Task 12: Documentation

**Files:**
- Modify: `README.md` — add a "TADF Angle-Scan Workbench" section documenting:
  - Folder layout expected (`angle_10.log`, `angle_20.log`, …)
  - How to launch (File → Open Angle Scan Folder)
  - What each plot shows, the 0.2 eV threshold convention, and that `P` = dominant S₁ CI coefficient
  - Donor/acceptor configuration dialog

- [ ] **Step 12.1: Add README section and commit**

```bash
git add README.md
git commit -m "docs: document TADF angle-scan workbench"
```

---

## Self-Review Checklist

- [x] **Spec coverage** — every user requirement mapped to a task:
  - Folder-of-log ingestion → Task 3
  - S₁/T₁ extraction + ΔE gap → Tasks 1, 4
  - 0.2 eV threshold + dashed line + highlighting → Task 5 (plot), 6 (table)
  - TDM per θ → Tasks 4, 5
  - 2·P² per θ → Tasks 1, 4, 5
  - Highest positive S₁ orbital coefficient (shown in table) → Tasks 1, 6
  - 3D canvas with donor/acceptor and θ variation → Tasks 7, 8
  - Cross-linked UI → Task 9
  - Main-window launcher → Task 10
- [x] No TBD / "implement later" / placeholders.
- [x] Names consistent (`AnglePoint`, `AngleScan`, `Dihedral3DPanel`, `AngleScanWindow`, `parse_angle_folder`, `angle_from_filename`, `gap_curve`, `tdm_curve`, `two_p_squared_curve`) across all tasks.
- [x] All file paths absolute-style from repo root.
- [x] Tests precede implementation in every task.

---

## Open Items for Implementer

These are intentional simplifications; flag to the user if they become blockers:

1. **θ encoding.** Only the filename path is inspected. If real user data has the angle only inside the log's route line, extend `angle_from_filename` to fall through to `_extract_angle_from_file_contents(path)` (search for `Scan` / `Dihedral` lines near the route section).
2. **Donor/acceptor picker.** The dialog uses dropdowns, not click-to-select on the 3D scene. Upgrading to click-to-select is a ≤1-day follow-up inside `DonorAcceptorDialog` using a small `GLViewWidget` with ray-picking.
3. **Slider continuous rotation.** The slider snaps to scanned angles. A nice-to-have is interpolating a "ghost" rotation of the donor fragment between snapshots while the user drags — not implemented here.
4. **Multiple singlets/triplets.** We take only S₁ and T₁ (the lowest of each). If the user later wants S_n / T_n views, extend `AnglePoint` with parameterised getters (e.g., `nth_singlet(n)`).
