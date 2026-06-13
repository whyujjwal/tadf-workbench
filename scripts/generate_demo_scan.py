"""Generate a chemistry-realistic angle-scan demo using biphenyl geometry.

Biphenyl (two phenyl rings joined by a single C–C bond) is the textbook
example of a rotatable inter-ring dihedral. The two rings ARE the donor and
acceptor fragments; the C₁–C₇ bond IS the rotatable single bond.

For each θ in the scan, we:
  1. Take the planar biphenyl reference geometry.
  2. Rotate the donor ring (and its H's) around the C₁–C₇ axis by (θ − 50°).
  3. Synthesise S₁/T₁ energies in a TADF-like profile (gap minimum near 50°).
  4. Emit a Gaussian-format .log with Standard orientation + excited states.

Output: data/demo_scan/angle_{N}.log for N ∈ {10, 20, …, 90}.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tadf_workbench.core.stage import Atom  # noqa: E402

OUT_DIR = ROOT / "data" / "demo_scan"

# Reference rotation = 0° → biphenyl is planar at this rotation; scan θ then
# IS (up to sign) the measured N₁-B₁-B₂-N₂ dihedral angle.
REFERENCE_ANGLE = 0.0

# ΔE(S₁−T₁) profile: minimum near 50°, bowls up on either side
ANGLE_TO_ENERGIES: dict[int, tuple[float, float]] = {
    10: (3.40, 2.88),
    20: (3.32, 2.91),
    30: (3.22, 2.95),
    40: (3.10, 2.98),
    50: (3.03, 3.00),
    60: (3.08, 2.99),
    70: (3.17, 2.97),
    80: (3.28, 2.94),
    90: (3.40, 2.90),
}


# ── Biphenyl construction ────────────────────────────────────────────────────

def build_biphenyl() -> tuple[list[Atom], tuple[int, int], list[int]]:
    """Build a planar biphenyl geometry.

    Returns (atoms, rotatable_bond, donor_atom_indices) where
    rotatable_bond = (C_acc, C_don) and donor_atom_indices are the indices of
    every atom that should rotate together with the donor ring.

    Atom layout (22 atoms total):
      0..5  : acceptor ring carbons (C-A1 .. C-A6); index 0 is the bond C
      6..11 : donor ring carbons (C-D1 .. C-D6);   index 6 is the bond C
      12..16: H's on acceptor ring (one per non-bonded C, indices 1..5)
      17..21: H's on donor ring    (one per non-bonded C, indices 7..11)
    """
    R = 1.40           # aromatic ring radius (C-C ~1.40 Å)
    BOND = 1.48        # inter-ring single C-C bond
    H_LEN = 1.08
    atoms: list[Atom] = []

    def hexagon(center: np.ndarray, start_deg: float) -> list[np.ndarray]:
        return [
            center + R * np.array([
                math.cos(math.radians(start_deg + 60 * k)),
                math.sin(math.radians(start_deg + 60 * k)),
                0.0,
            ])
            for k in range(6)
        ]

    # Acceptor ring: hexagon centered at origin; bond C at +x (vertex 0)
    acc_center = np.zeros(3)
    acc_carbons = hexagon(acc_center, start_deg=0.0)

    # Donor ring: shifted along +x by (R + BOND + R), bond C is closest (vertex 0
    # of the donor ring rotated 180° so its first vertex points back toward acceptor)
    don_center = np.array([R + BOND + R, 0.0, 0.0])
    don_carbons = hexagon(don_center, start_deg=180.0)

    for c in acc_carbons:
        atoms.append(Atom(atomic_number=6, x=c[0], y=c[1], z=c[2]))
    for c in don_carbons:
        atoms.append(Atom(atomic_number=6, x=c[0], y=c[1], z=c[2]))

    # H's — one per non-bonded carbon (skip vertex 0 in each ring), pointing
    # radially outward.
    h_indices: list[int] = []
    for i, c in enumerate(acc_carbons):
        if i == 0:
            continue  # vertex 0 is the bond C — no H here
        d = c - acc_center
        d = d / np.linalg.norm(d)
        h = c + H_LEN * d
        atoms.append(Atom(atomic_number=1, x=h[0], y=h[1], z=h[2]))
    for i, c in enumerate(don_carbons):
        if i == 0:
            continue
        d = c - don_center
        d = d / np.linalg.norm(d)
        h = c + H_LEN * d
        atoms.append(Atom(atomic_number=1, x=h[0], y=h[1], z=h[2]))

    # Re-index atoms
    for idx, a in enumerate(atoms):
        a.index = idx

    rotatable_bond = (0, 6)  # acceptor-bond-C ↔ donor-bond-C
    # Donor side: donor ring carbons (6..11) + their H's (17..21)
    donor_indices = list(range(6, 12)) + list(range(17, 22))
    return atoms, rotatable_bond, donor_indices


# ── Geometry helpers ─────────────────────────────────────────────────────────

def rodrigues_rotate(points: np.ndarray, axis_p1: np.ndarray,
                     axis_p2: np.ndarray, angle_rad: float) -> np.ndarray:
    axis = axis_p2 - axis_p1
    axis = axis / np.linalg.norm(axis)
    translated = points - axis_p1
    c, s = math.cos(angle_rad), math.sin(angle_rad)
    rotated = (translated * c
               + np.cross(axis, translated) * s
               + np.outer(translated @ axis, axis) * (1 - c))
    return rotated + axis_p1


# ── Gaussian log emission ────────────────────────────────────────────────────

def format_orientation_block(atoms: list[Atom]) -> str:
    lines = [
        "                         Standard orientation:                         ",
        " ---------------------------------------------------------------------",
        " Center     Atomic      Atomic             Coordinates (Angstroms)",
        " Number     Number       Type             X           Y           Z",
        " ---------------------------------------------------------------------",
    ]
    for atom in atoms:
        lines.append(
            f"{atom.index + 1:>7}{atom.atomic_number:>11}{0:>12}"
            f"{atom.x:>16.6f}{atom.y:>12.6f}{atom.z:>12.6f}"
        )
    lines.append(" ---------------------------------------------------------------------")
    return "\n".join(lines) + "\n"


def format_quantum_block(s1_energy_ev: float, t1_energy_ev: float,
                         tdm: tuple[float, float, float]) -> str:
    s1_nm = 1239.84 / s1_energy_ev
    t1_nm = 1239.84 / t1_energy_ev
    return (
        "\n"
        f" Excited State   1:      Triplet-A            {t1_energy_ev:.4f} eV  "
        f"{t1_nm:.2f} nm  f=0.0000  <S**2>=2.000\n"
        f"      52 -> 55         0.64482\n"
        f"      49 -> 56         0.22000\n"
        "\n"
        f" Excited State   2:      Singlet-A            {s1_energy_ev:.4f} eV  "
        f"{s1_nm:.2f} nm  f=0.0500  <S**2>=0.000\n"
        f"      52 -> 55         0.65880\n"
        f"      49 -> 56         0.20000\n"
        f"      49 -> 71        -0.10500\n"
        "\n"
        " Ground to excited state transition electric dipole moments (Au):\n"
        "       state          X           Y           Z        Dip. S.      Osc.\n"
        "         1         0.00000    0.00000    0.00000    0.00000    0.0000\n"
        f"         2         {tdm[0]:.5f}    {tdm[1]:.5f}    {tdm[2]:.5f}    "
        "0.14000    0.0500\n"
        " End of table.\n"
    )


def tdm_for_angle(angle_deg: float) -> tuple[float, float, float]:
    mag = 0.25 + 0.30 * abs(math.cos(math.radians(angle_deg - 50)))
    return (mag * 0.6, mag * 0.5, mag * 0.4)


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    base_atoms, bond, donor_idx = build_biphenyl()
    print(f"biphenyl built: {len(base_atoms)} atoms, "
          f"rotatable bond ({bond[0]}, {bond[1]}), "
          f"donor fragment {len(donor_idx)} atoms ({donor_idx[0]}..{donor_idx[-1]})")

    base_pos = np.array([(a.x, a.y, a.z) for a in base_atoms])
    axis_p1 = base_pos[bond[0]]
    axis_p2 = base_pos[bond[1]]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for path in OUT_DIR.glob("*.log"):
        path.unlink()

    for angle_deg, (s1_e, t1_e) in ANGLE_TO_ENERGIES.items():
        delta = math.radians(angle_deg - REFERENCE_ANGLE)
        positions = base_pos.copy()
        positions[donor_idx] = rodrigues_rotate(
            base_pos[donor_idx], axis_p1, axis_p2, delta,
        )
        rotated_atoms = [
            Atom(atomic_number=base_atoms[i].atomic_number,
                 x=float(positions[i, 0]),
                 y=float(positions[i, 1]),
                 z=float(positions[i, 2]),
                 index=i)
            for i in range(len(base_atoms))
        ]
        content = (
            f" Gaussian 16 synthetic biphenyl scan point θ={angle_deg}°\n"
            " (donor = ring2 / atoms 6-11+17-21; acceptor = ring1 / atoms 0-5+12-16)\n"
            " ---------------------------------------------------------------------\n"
            + format_orientation_block(rotated_atoms)
            + format_quantum_block(s1_e, t1_e, tdm_for_angle(angle_deg))
        )
        out = OUT_DIR / f"angle_{angle_deg}.log"
        # Explicit UTF-8: the content contains θ/°, which crash on platforms
        # whose default text encoding isn't UTF-8 (e.g. Windows cp1252). The
        # parser reads these files as UTF-8, so write them the same way.
        out.write_text(content, encoding="utf-8")
        gap = s1_e - t1_e
        tag = " TADF" if gap <= 0.2 else ""
        print(f"  wrote {out.name}  ΔE={gap:+.3f} eV{tag}")
    print(f"done. scan folder: {OUT_DIR}")


if __name__ == "__main__":
    main()
