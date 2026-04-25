# User Guide

## Step 1 — prepare a scan folder

The workbench ingests a folder containing one Gaussian `.log` (or `.out`)
file per scanned dihedral angle. The **first integer in each filename** is
the angle in degrees:

```
my_donor_acceptor_scan/
├── angle_10.log
├── angle_20.log
├── angle_30.log
├── angle_40.log
├── angle_50.log
├── angle_60.log
├── angle_70.log
├── angle_80.log
└── angle_90.log
```

All these patterns work too: `10.log`, `scan_d45deg.out`, `mol_60_opt.log`,
`theta-75.log`. Files without a detectable integer are skipped (use **File
→ Import Log Files…** to add them with a manual angle prompt).

Each log must contain at minimum:
1. A **Standard orientation** block (the molecular geometry).
2. At least one **`Excited State N: Singlet-*`** block.
3. At least one **`Excited State N: Triplet-*`** block.
4. A **Ground to excited state transition electric dipole moments (Au)**
   block (so we can read the S₁ TDM).

A `td=(50-50,nstates=10) b3lyp/...` route at the right dihedral generally
produces all of the above.

## Step 2 — launch

```bash
# from the project root
python run.py /path/to/my_donor_acceptor_scan
```

If no folder is supplied, the workbench opens empty — use **File → Open
Scan Folder…** (Ctrl+Shift+O).

## Step 3 — configure the rotatable bond

The first time a scan is loaded, the **Configure Rotatable Bond** dialog
appears. Two steps:

1. **Pick the rotatable bond** — choose the two atoms (B₁, B₂) that form
   the single bond between donor and acceptor. The dialog uses the first
   scan point's geometry.
2. **Confirm donor / acceptor assignment.** As soon as both atoms are
   selected, the dialog reports

   > ✓ Bond is rotatable — molecule splits into 11 + 11 atoms.

   and shows two side-cards with the atom count, formula (`C6H5`), anchor
   atom and reference atom for each side. If the workbench guessed the
   wrong side as donor, click **⇄ Swap donor / acceptor** to flip them.

If you pick a bond that's part of a ring, the dialog says

> ✗ bond X-Y is part of a ring; cutting it does not split the molecule.

Pick a different (single, non-ring) bond.

You can revisit this dialog later via **Scan → Configure Donor/Acceptor…**.

## Step 4 — read the dashboard

```
┌─ Stat chips ───────────────────────────────────────────────────┐
│ CURRENT θ │ ΔE AT θ │ TADF │ MIN ΔE @ θ │ |TDM| │ 2·P²        │
├──────────────────────────────────┬─────────────────────────────┤
│                                  │ Current State Card          │
│                                  ├─────────────────────────────┤
│   3D molecule viewer             │ ΔE(S₁−T₁) vs θ  (hero plot) │
│   + θ slider                     ├─────────────────────────────┤
│   + measured dihedral readout    │ |TDM(S₁)| vs θ              │
│                                  ├─────────────────────────────┤
│                                  │ 2·P² vs θ                   │
├──────────────────────────────────┴─────────────────────────────┤
│ Per-θ data table (TADF rows highlighted green)                 │
└────────────────────────────────────────────────────────────────┘
```

### 3D viewer (left)

- **Donor atoms** are blended toward green; **acceptor atoms** toward
  orange. Element colour is preserved underneath, so a carbon is still
  recognisable.
- The **rotatable bond is yellow** with a translucent extended axis so you
  can see the rotation axis even when the rings are nearly coplanar.
- Two **semi-transparent quads** (green for donor, orange for acceptor) are
  the best-fit planes through each fragment.
- The **yellow arc** between the donor and acceptor reference atoms shows
  the dihedral; its sweep visually equals the angle.
- The bottom bar shows two angles:
  - **SCAN θ** (blue) — the scanned value from the filename.
  - **MEASURED N₁–B₁–B₂–N₂** (yellow) — the IUPAC dihedral computed live
    from atom positions.

The slider at the bottom snaps to the nearest scanned angle. Mouse: left-
drag rotates the camera, right-drag pans, wheel zooms.

### Stat chips (top)

Six compact chips:

| Chip            | Shows                                                       |
|-----------------|-------------------------------------------------------------|
| `CURRENT θ`     | The selected scan angle                                     |
| `ΔE AT θ`       | The S₁−T₁ gap at the selected angle (green if TADF)         |
| `TADF`          | YES / no — whether the current ΔE is ≤ threshold            |
| `MIN ΔE @ θ`    | The smallest ΔE in the whole scan and the angle at which it occurs |
| `\|TDM\|`        | S₁ transition dipole magnitude at the current angle         |
| `2·P²`          | Dominant CI weight for S₁ at the current angle              |

### Current State card

Side-by-side **S₁** and **T₁** boxes (energy in eV, wavelength in nm,
oscillator strength for S₁, ⟨S²⟩ for T₁), a coloured **ΔE banner** (green ≤
threshold, amber ≤ 2·threshold, red otherwise), the **dominant MO
transition** for S₁ (`MO i → MO a (coefficient)`), the |TDM(S₁)|, and 2·P².

### Plots (right column)

- **ΔE(S₁−T₁) vs θ** — the hero plot. Solid line through all scan points;
  dashed orange line at the threshold; large green markers on points at or
  below the threshold. Click any point to select that angle.
- **|TDM(S₁)| vs θ** — secondary plot in yellow. Larger TDM → brighter
  emission once you escape T₁.
- **2·P² vs θ** — secondary plot in violet. Indicates how
  single-configuration the S₁ state is at each angle. A drop here usually
  flags a state with strong configurational mixing.

### Data table (bottom)

One row per scanned angle with: θ, E(S₁), E(T₁), ΔE, |TDM|, 2·P², and the
dominant `MO i → MO a (c)`. **TADF rows are highlighted green and bold.**
Clicking a row selects that angle across the whole UI.

## Step 5 — find the best angle

The workflow:

1. Look at the gap plot. The bowl-shaped minimum tells you the rough
   region.
2. Read **MIN ΔE @ θ** in the top stat bar — it points you straight at the
   best scanned angle.
3. Click that θ in the table or drag the slider there. Confirm the **TADF
   = YES** badge is green.
4. Check the |TDM| plot at the same angle: a TADF candidate with vanishing
   TDM emits weakly. The trade-off is real and the workbench surfaces both
   axes.
5. Use **Scan → Set ΔE Threshold…** to tighten or relax the cutoff.
6. **File → Export Summary CSV…** dumps the full table for downstream
   reporting.

## Importing log files one at a time

**File → Import Log Files…** (Ctrl+I) opens a multi-select picker. For
each file:
- The angle is auto-extracted from the filename if possible, otherwise the
  app prompts you with a "what θ is this?" dialog.
- If the angle is already in the current scan, the app asks once whether
  to **replace this and all later conflicts** or **skip all later
  conflicts** for the rest of the batch.
- After import, the status bar reports `Import: N added · M replaced · K
  skipped`.

If no scan is loaded yet, the imported files become a fresh scan and the
configuration dialog opens automatically.

## Tips

- The workbench never modifies your log files. Everything is in-memory.
- Camera not framed well? Resize the window slightly — the camera auto-fits
  on the first render after a fresh scan load.
- For a clean slate: **File → Open Scan Folder…** loads a new folder
  without merging.
- The bundled biphenyl demo (`data/demo_scan/`) is regenerated by
  `python scripts/generate_demo_scan.py`. Edit `ANGLE_TO_ENERGIES` in that
  script to play with different gap profiles.
