# Architecture

## Module layout

```
src/tadf_workbench/
├── core/
│   ├── stage.py           # Atom, Stage (one geometry snapshot)
│   ├── molecule.py        # Molecule (multiple stages + quantum_data)
│   ├── quantum_data.py    # ExcitedState, TransitionDipole, MullikenCharges, …
│   ├── angle_scan.py      # AnglePoint, AngleScan (one row of the dashboard)
│   └── fragments.py       # bond graph, BFS split, IUPAC dihedral, plane fit
├── parsers/
│   ├── gaussian_parser.py # parses ONE Gaussian .log → Molecule + QuantumData
│   └── angle_scan_parser.py
│                          # parse_angle_folder() + angle_from_filename()
├── analysis/
│   └── tadf_metrics.py    # gap_curve, tdm_curve, two_p_squared_curve, summary_row
├── ui/
│   ├── angle_scan_window.py     # top-level QMainWindow assembling everything
│   ├── angle_scan_plots.py      # GapPlot, TDMPlot, TwoPSquaredPlot
│   ├── angle_scan_table.py      # AngleScanTable
│   ├── dihedral_3d_panel.py     # 3D viewer + slider + dihedral readout
│   └── donor_acceptor_dialog.py # rotatable-bond + fragment configurator
├── utils/
│   └── atom_properties.py # element symbols, colours, radii, bond thresholds
└── main.py                # entry point
```

## Data flow (single scan load)

```
folder/  ──▶  parse_angle_folder()                               (parsers)
                │
                │   for each *.log/*.out in folder:
                │     angle = angle_from_filename(path)
                │     molecule = parse_gaussian_file(path)        (parsers)
                │     point = AnglePoint(angle, path, molecule)   (core)
                │
                ▼
            AngleScan(name, points)                              (core)
                │
                │   user picks rotatable bond in dialog:
                │     scan.configure_fragments(stage)              (core)
                │       └─ split_at_bond → donor / acceptor sets
                │       └─ pick_reference_atom → N₁, N₂
                │
                ▼
        AngleScanWindow.set_scan(scan)                           (ui)
            ├─▶ GapPlot.set_scan(scan)
            ├─▶ TDMPlot.set_scan(scan)            ┐
            ├─▶ TwoPSquaredPlot.set_scan(scan)    ├─ all use analysis/*_curve()
            ├─▶ AngleScanTable.set_scan(scan)     ┘   from a sorted view
            └─▶ Dihedral3DPanel.set_scan(scan)
                    └─ renders donor/acceptor planes, dihedral arc
```

Selection events flow back the other way: a click on a plot, a row pick on
the table, or a slider drag all reach `AngleScanWindow._on_angle_selected`,
which in turn updates every other panel.

## Key abstractions

### `AnglePoint` ([`core/angle_scan.py`](../src/tadf_workbench/core/angle_scan.py))

Frozen-by-convention dataclass holding one scan point. All derived metrics
are properties that read from `self.molecule.quantum_data`:

| Property                     | What it computes                                     |
|------------------------------|------------------------------------------------------|
| `s1_state`                  | Lowest-numbered `Singlet-*` excited state            |
| `t1_state`                  | Lowest-numbered `Triplet-*` excited state            |
| `s1_t1_gap_ev`              | `E(S₁) − E(T₁)` in eV                                |
| `s1_dominant_coefficient`   | Largest \|c\| in S₁'s orbital-transition list (= P)  |
| `s1_dominant_transition`    | The full `(i, a, c)` tuple                           |
| `two_p_squared`             | `2 · P²`                                             |
| `s1_tdm`                    | TransitionDipole row matching `s1_state.number`      |
| `s1_tdm_magnitude`          | \|TDM\| in atomic units                              |

Returning `None` for any of the above is normal when the underlying log
lacks the data; the UI handles it gracefully.

### `AngleScan` ([`core/angle_scan.py`](../src/tadf_workbench/core/angle_scan.py))

Aggregates the points and the donor/acceptor configuration:

```python
@dataclass
class AngleScan:
    name: str
    points: List[AnglePoint]
    threshold_ev: float = 0.2
    donor_atom_index: Optional[int] = None
    acceptor_atom_index: Optional[int] = None
    rotatable_bond: Optional[Tuple[int, int]] = None
    donor_fragment: Optional[Set[int]] = None
    acceptor_fragment: Optional[Set[int]] = None
    donor_reference_atom: Optional[int] = None
    acceptor_reference_atom: Optional[int] = None
```

Helpers: `sorted_points`, `gaps_ev`, `tdm_magnitudes`,
`two_p_squared_values`, `points_below_threshold(t)`, `minimum_gap_point`,
`nearest_point(θ)`, `find_point(θ)`, `add_or_replace(point)`,
`configure_fragments(stage)`.

### Fragment detection ([`core/fragments.py`](../src/tadf_workbench/core/fragments.py))

```python
build_adjacency(stage)           # → dict[int, set[int]]
split_at_bond(stage, a, b)       # BFS from each side; raises BondNotRotatableError on rings
pick_reference_atom(stage, frag, anchor)
                                 # heaviest in-fragment neighbour, prefers high-degree
dihedral_angle(p1, p2, p3, p4)   # IUPAC, signed, in (-180°, 180°]
fit_plane(positions)             # → (centroid, unit_normal) via SVD
plane_quad(centroid, normal,     # 4-corner square in that plane
           half_extent, in_plane_hint)
```

These five primitives are the entire chemistry layer. They are pure (no Qt,
no I/O) and individually unit-tested in
[`tests/test_fragments.py`](../tests/test_fragments.py).

### Parsers

`parse_gaussian_file(path) → Molecule` is the existing single-log parser,
unchanged from the upstream project. It extracts every Standard /
Input / Z-Matrix orientation block, every Excited State block, every
transition-dipole block, dipole/quadrupole moments, and Mulliken charges.

`parse_angle_folder(folder) → AngleScan` walks a directory, calls
`angle_from_filename` per file, parses each via the single-log parser, and
collects into an `AngleScan`. Files without a detectable angle or with
parse errors are skipped (printed to stdout, scan continues).

`angle_from_filename(path) → Optional[float]` — current implementation
returns the **first integer in the filename stem**. See `parsers/angle_scan_parser.py`.

### Analysis

`analysis/tadf_metrics.py` is a thin layer that converts an `AngleScan`
into plot-ready tuples and table-ready dicts. Keeping it separate from the
UI panels means the UI has no math in it — the panels just call
`gap_curve(scan)` and plot the result.

## UI architecture

### `AngleScanWindow` ([`ui/angle_scan_window.py`](../src/tadf_workbench/ui/angle_scan_window.py))

Top-level QMainWindow. Lays out:

```
QVBoxLayout
├── stat-chip bar (QHBoxLayout of StatChip)
├── QSplitter(horizontal)
│   ├── 3D viewer (Dihedral3DPanel, hero — wide)
│   └── sidebar (QVBoxLayout)
│       ├── CurrentStateCard
│       └── QSplitter(vertical)
│           ├── GapPlot
│           ├── TDMPlot
│           └── TwoPSquaredPlot
└── AngleScanTable (compact, full width)
```

Cross-wires every `angle_clicked` / `angle_selected` / `angle_changed` signal
to a single `_on_angle_selected(angle)` slot that in turn updates all
other panels. The `_refresh_chips` method is the single source of truth
for the stat-chip bar and the Current State card.

### `Dihedral3DPanel` ([`ui/dihedral_3d_panel.py`](../src/tadf_workbench/ui/dihedral_3d_panel.py))

Wraps `pyqtgraph.opengl.GLViewWidget`. The render pipeline per scan point:

1. `_fit_camera_to_stage` (first frame only) — sets the camera distance and
   target from the molecular bounding box.
2. `_draw_atoms` — one `GLMeshItem` sphere per atom, colour blended toward
   the green or orange fragment tint depending on which side the atom is on.
3. `_draw_bonds` — every non-rotatable bond as a single `GLLinePlotItem`
   line set; the rotatable bond is excluded so it can be drawn separately.
4. `_draw_rotatable_axis` — thick yellow segment for the bond, plus a
   translucent extended axis line.
5. `_draw_fragment_plane(donor) / (acceptor)` — fits a plane via SVD and
   draws a translucent two-triangle mesh quad through it.
6. `_draw_dihedral` — computes the IUPAC dihedral, updates the **MEASURED**
   label, and draws an arc + guide lines in the plane perpendicular to the
   bond axis.

The slider snaps to the nearest scanned angle via `AngleScan.nearest_point`
and emits `angle_changed`.

### `DonorAcceptorDialog` ([`ui/donor_acceptor_dialog.py`](../src/tadf_workbench/ui/donor_acceptor_dialog.py))

Two-step UX: pick the rotatable bond → confirm/swap fragment assignment.
Updates `scan.rotatable_bond`, `donor_fragment`, `acceptor_fragment`,
`donor_atom_index`, `acceptor_atom_index`, `donor_reference_atom`,
`acceptor_reference_atom` on accept. Has a backwards-compatibility hook
`set_selection(donor, acceptor, bond)` used by older tests with tiny
stub molecules where bond perception fails.

### `AngleScanTable`, `GapPlot`, `TDMPlot`, `TwoPSquaredPlot`

Thin presentational widgets. Each takes an `AngleScan`, calls into
`analysis.tadf_metrics`, and exposes one signal:
`angle_selected(float)` (table) or `angle_clicked(float)` (plots).

## Tests

73 tests across 11 files, all run headless under `QT_QPA_PLATFORM=offscreen`:

| File                                | Coverage                                    |
|-------------------------------------|---------------------------------------------|
| `test_angle_scan_core.py`           | AnglePoint properties, AngleScan helpers    |
| `test_fragments.py`                 | Bond graph, BFS split, dihedral, plane fit  |
| `test_angle_scan_parser.py`         | Filename → angle, folder ingestion          |
| `test_tadf_metrics.py`              | gap_curve, tdm_curve, two_p_squared_curve   |
| `test_angle_scan_plots.py`          | Plot widget population + click signal       |
| `test_angle_scan_table.py`          | Table population + TADF-row highlighting    |
| `test_donor_acceptor_dialog.py`     | Selection validation + backwards compat     |
| `test_dihedral_3d_panel.py`         | Slider snap, signal dedup                   |
| `test_angle_scan_window.py`         | Cross-wiring, status bar                    |
| `test_e2e_angle_scan.py`            | Full end-to-end flow with 9-angle scan      |
| `tests/fixtures/build_fixtures.py`  | Synthetic Gaussian-format log generator     |

## Threading

Single-threaded. Every operation runs on the Qt main thread. Parsing the
9-point demo scan is a few tens of milliseconds; full scans of dozens of
points with rich logs are still well under a second on a laptop. If that
ever stops being true, route `parse_angle_folder` through `QThreadPool` —
it's a pure function that can run off-thread without changes.

## Conventions

- Public API is exported from each subpackage's `__init__.py`. Internal
  helpers (`_underscored`) stay module-local.
- pure-logic modules (`core/fragments.py`, `analysis/tadf_metrics.py`)
  import only `numpy` and stdlib — no Qt. Easy to test and reuse.
- UI modules import core and analysis but never the parsers directly,
  except `AngleScanWindow` which orchestrates loading.
- Type hints on public functions; missing on small private helpers is fine.
