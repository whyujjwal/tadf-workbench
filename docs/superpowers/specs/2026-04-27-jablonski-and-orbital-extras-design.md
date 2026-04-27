# Design — Jablonski Diagram + Orbital-Transition Extras

Date: 2026-04-27
Branch: main

## Goal

Extend the existing TADF angle-scan workbench with three additions on top of
what already ships in the dashboard:

1. **Jablonski diagram per θ** — energy-level diagram (S₀, all singlets, all
   triplets) with **S₁**, **f-max singlet** ("Sₙ"), and **T₁** highlighted.
   Available as a live dock that updates with the selected angle, and as an
   all-angles grid window that shows every θ side-by-side.
2. **Top-2 orbital transitions for S₁** — currently the workbench shows only
   the dominant `i→a` jump for S₁; we add the **second-highest** as well, in
   both the per-θ table and the Jablonski labels.
3. **Two modes for the 2·P² plot** — keep the existing `P = dominant CI
   coefficient` view (configuration weight) and add a `P = S₁ energy in eV`
   view, switchable from a small toggle, with an ⓘ tooltip explaining both.

## What's already done (context)

The repo already parses Gaussian logs, picks S₁/T₁ as the first Singlet-A /
first Triplet-A, plots θ vs ΔE / TDM / 2·P², shows a 3D viewer with the
acceptor fixed and donor rotated, and ships a per-θ table with CSV export.

This spec only describes the additions. The existing files we touch are
limited to wiring (`AngleScanWindow`), a tweak to `TwoPSquaredPlot`, and one
new column on `AngleScanTable`.

## Approach

Approach 2 from brainstorming: **new `analysis/jablonski.py` + new UI files
under `ui/`**, wired into `AngleScanWindow` via a dock and a menu. Leaves the
existing dashboard layout intact and keeps new code in isolated, testable
units (per `docs/ARCHITECTURE.md`'s separation of `parsers/` / `analysis/` /
`ui/`).

## Component breakdown

### Data layer — `core/angle_scan.py`

Add three lazy properties to `AnglePoint`:

```python
@property
def s1_top2_transitions(self) -> List[Tuple[int, int, float]]:
    """Top 2 orbital transitions for S₁ ranked by |coefficient|, descending.
    Returns up to 2 entries; empty list if no S₁ or no transitions."""

@property
def fmax_singlet_state(self) -> Optional[ExcitedState]:
    """Singlet (any number) with the highest oscillator strength.
    May coincide with S₁."""

@property
def fmax_singlet_top2_transitions(self) -> List[Tuple[int, int, float]]:
    """Top 2 transitions for the f-max singlet, same shape as s1_top2_transitions."""
```

These are pure derivations from `molecule.quantum_data.excited_states` and
have no side effects.

### Analysis — `analysis/jablonski.py` (new file)

Pure functions that compute layout data; no Qt imports.

```python
@dataclass
class Level:
    energy_ev: float
    multiplicity: str        # "Singlet" or "Triplet"
    state_number: int
    label: str               # "" for unhighlighted, "S₁: 96→101 65.4% / 92→103 6.0%" otherwise
    is_highlighted: bool
    column: str              # "singlet" or "triplet"

@dataclass
class Arrow:
    from_energy_ev: float
    to_energy_ev: float
    kind: str                # "absorption" | "isc"
    label: str               # "f=0.358" or "ΔE_ST=0.04 eV"
    column: str              # which column the arrow lives in

@dataclass
class JablonskiLayout:
    levels: List[Level]
    arrows: List[Arrow]
    title: str               # e.g. "θ = 30°"

def jablonski_layout(point: AnglePoint, mode: Literal["key", "all"]) -> JablonskiLayout: ...

def percent_2p2(coeff: float) -> float:
    """2·P² × 100 for percent display."""
    return 2.0 * coeff * coeff * 100.0
```

In `mode="key"` only S₀, S₁, the f-max singlet, and T₁ are emitted as
`Level`s. In `mode="all"` every singlet and triplet is emitted, with the
three above marked `is_highlighted=True` and labelled.

Two arrows always: `S₀→S₁` (label `f=…` from S₁'s oscillator strength) and
`S₀→Sₙ` (`f=…` from the f-max singlet). Plus an `S₁→T₁` ISC arrow labelled
with `ΔE_ST=…`.

### UI — `ui/jablonski_panel.py` (new)

`JablonskiPanel(QWidget)` — a single-angle Jablonski:

- `pyqtgraph.PlotWidget` background, two-column layout (singlets left,
  triplets right) on a shared eV y-axis.
- Top of panel: a small `Key only / All states` `QButtonGroup` of radio
  buttons; default = key-only. Re-renders on toggle.
- Slot `set_point(point: Optional[AnglePoint])`. With `None`, clears.
- Highlighted levels drawn as thick coloured bars with the multi-line label
  rendered next to them (singlet column = blue, triplet column = orange — same
  palette as the existing window's `SINGLET`/`TRIPLET` colours).
- Faint horizontal lines for non-highlighted states in `all` mode.
- Arrows rendered as `pg.ArrowItem` with the label as a `pg.TextItem` next to
  the arrow midpoint.

### UI — `ui/jablonski_grid_window.py` (new)

`JablonskiGridWindow(QMainWindow)` — opens non-modally from the View menu:

- Scrollable `QGridLayout` of small `JablonskiPanel`s, one per θ in the
  current scan.
- Each tile has a header label `θ=10° · ΔE=0.34 eV · f_max=0.36` above the
  panel.
- Tiles are click-to-select: panel emits `angle_picked(float)`; the grid
  window relays it; the main window connects it to `_on_angle_selected`.
- Defaults to `key`-mode for tiles (the "all" mode is too dense at small
  size).
- Tile size auto-fits to a 3×N grid by default (3 columns); user can resize
  the window.

### `ui/angle_scan_plots.py` — `TwoPSquaredPlot` revision

Add a `Mode` enum / string with two values: `"ci"` (default — current
behaviour) and `"energy"` (P = S₁ energy in eV). Tiny `QHBoxLayout` above the
plot containing two `QRadioButton`s and a small `ⓘ` `QToolButton` whose
tooltip reads:

> **2·P² modes**
>
> *CI coefficient (default):* P is the largest |coefficient| of S₁'s
> orbital transition list. 2·P² is the standard "weight of the dominant
> configuration" — dimensionless, ~0–1, used widely in TADF papers.
>
> *S₁ energy (eV):* P is S₁'s vertical excitation energy in eV. 2·P² has
> units of eV² and is non-standard; provided on user request.

Recompute the curve on toggle.

Note: this requires wrapping the current `pg.PlotWidget` subclass into a
`QWidget` containing the toggle + plot. Keep the same `set_scan` /
`angle_clicked` API so `AngleScanWindow` doesn't notice.

### `ui/angle_scan_table.py` — new column

Insert a new column `"S₁ 2nd jump"` between the existing `"S₁ dominant"`
column and the next column. Format `"93→105 7.8%"` (using `percent_2p2`).
Empty string when there's no second transition.

### `ui/angle_scan_window.py` — wiring

Three small additions:

1. In `_build_ui`, construct `self.jablonski_panel = JablonskiPanel()` and
   add it as a `QDockWidget` titled "Jablonski" on the right edge of the
   main window. Default visible.
2. In `_on_angle_selected`, after the existing `set_angle` call, also call
   `self.jablonski_panel.set_point(point_for_angle)`.
3. In `_build_menu`, add a `View` menu with two actions:
   - `Jablonski Panel` (checkable, toggles dock visibility).
   - `Jablonski Grid…` — opens `JablonskiGridWindow(self._scan, parent=self)`
     non-modally; connects its `angle_picked` signal to
     `self._on_angle_selected`.

## File touch list

New:
- `src/tadf_workbench/analysis/jablonski.py`
- `src/tadf_workbench/ui/jablonski_panel.py`
- `src/tadf_workbench/ui/jablonski_grid_window.py`

Edit:
- `src/tadf_workbench/core/angle_scan.py` (3 properties on `AnglePoint`)
- `src/tadf_workbench/analysis/__init__.py` (re-export jablonski helpers)
- `src/tadf_workbench/ui/__init__.py` (re-export new UI classes)
- `src/tadf_workbench/ui/angle_scan_plots.py` (`TwoPSquaredPlot` toggle)
- `src/tadf_workbench/ui/angle_scan_table.py` (new column)
- `src/tadf_workbench/ui/angle_scan_window.py` (dock + View menu)

## Out of scope

- New tests (per user direction — existing 73 tests must keep passing,
  but we don't write new ones).
- Changes to the parser / quantum_data — all data needed is already
  extracted.
- Changes to the 3D viewer or fragments logic.
- Persistence of the Jablonski mode toggle / dock visibility across
  sessions.
