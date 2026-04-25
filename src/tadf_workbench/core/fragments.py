"""Fragment detection + dihedral geometry for donor–acceptor scans.

A donor–acceptor molecule has one rotatable single bond. Cutting that bond
splits the molecular graph into two connected components — the donor side
and the acceptor side. The dihedral angle θ is measured between a chosen
reference atom on each side, around the rotatable bond axis.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from .stage import Stage


# ── Bond graph ────────────────────────────────────────────────────────────────

def build_adjacency(stage: Stage, tolerance: float = 0.3) -> Dict[int, Set[int]]:
    """Adjacency dict: atom_index -> set of bonded atom_indices."""
    adj: Dict[int, Set[int]] = {i: set() for i in range(stage.atom_count)}
    for i, j in stage.find_bonds(tolerance=tolerance):
        adj[i].add(j)
        adj[j].add(i)
    return adj


# ── Fragment splitting ────────────────────────────────────────────────────────

class BondNotRotatableError(ValueError):
    """Raised when the requested bond is in a ring (cutting it doesn't split
    the molecule into two disconnected components)."""


def split_at_bond(stage: Stage, atom_a: int, atom_b: int,
                  tolerance: float = 0.3) -> Tuple[Set[int], Set[int]]:
    """Split the molecular graph by removing the bond (atom_a, atom_b).

    Returns (side_a_atoms, side_b_atoms) — atoms reachable from `atom_a`
    (excluding the path through `atom_b`) and from `atom_b` respectively.
    Each side INCLUDES its own anchor atom.

    Raises BondNotRotatableError if the two ends remain connected after the
    cut (i.e. the bond is part of a ring).
    """
    if atom_a == atom_b:
        raise ValueError("rotatable bond requires two distinct atoms")
    adj = build_adjacency(stage, tolerance=tolerance)
    if atom_b not in adj.get(atom_a, set()):
        # Not actually bonded — still allow but warn via a clear failure:
        raise BondNotRotatableError(
            f"atoms {atom_a} and {atom_b} are not bonded in the current geometry"
        )

    def reach_from(start: int, blocked_edge: Tuple[int, int]) -> Set[int]:
        visited = {start}
        stack = [start]
        while stack:
            u = stack.pop()
            for v in adj[u]:
                # Skip the cut edge in either direction
                if {u, v} == set(blocked_edge):
                    continue
                if v not in visited:
                    visited.add(v)
                    stack.append(v)
        return visited

    cut = (atom_a, atom_b)
    side_a = reach_from(atom_a, cut)
    side_b = reach_from(atom_b, cut)

    if side_a & side_b:
        raise BondNotRotatableError(
            f"bond {atom_a}-{atom_b} is part of a ring; "
            "cutting it does not split the molecule"
        )
    return side_a, side_b


# ── Reference-atom picking for dihedral ───────────────────────────────────────

def pick_reference_atom(stage: Stage, fragment: Set[int], anchor: int,
                        tolerance: float = 0.3) -> Optional[int]:
    """Pick a heavy-atom (Z > 1) neighbor of `anchor` that lies inside
    `fragment`, preferring sp²/aromatic-looking environments (degree ≥ 2).

    Returns None only if no neighbor exists inside the fragment.
    """
    adj = build_adjacency(stage, tolerance=tolerance)
    neighbors = [n for n in adj[anchor] if n in fragment and n != anchor]
    if not neighbors:
        return None
    heavy = [n for n in neighbors if stage.atoms[n].atomic_number > 1]
    pool = heavy or neighbors
    # Prefer the neighbor with the highest in-fragment degree (more "planar-y")
    return max(pool, key=lambda n: len(adj[n] & fragment))


# ── Dihedral math ─────────────────────────────────────────────────────────────

def dihedral_angle(p1: np.ndarray, p2: np.ndarray,
                   p3: np.ndarray, p4: np.ndarray) -> float:
    """IUPAC dihedral angle (degrees) of points p1–p2–p3–p4.

    The bond p2→p3 is the central axis; the angle is measured from the
    plane (p1,p2,p3) to the plane (p2,p3,p4), signed by the right-hand rule
    around p2→p3. Returns a value in (-180, 180].
    """
    b1 = np.asarray(p2) - np.asarray(p1)
    b2 = np.asarray(p3) - np.asarray(p2)
    b3 = np.asarray(p4) - np.asarray(p3)
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    m1 = np.cross(n1, b2 / np.linalg.norm(b2))
    x = np.dot(n1, n2)
    y = np.dot(m1, n2)
    return float(np.degrees(np.arctan2(y, x)))


# ── Plane fitting (for the visual donor/acceptor planes) ──────────────────────

def fit_plane(positions: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Best-fit plane through the points (rows are 3D coords).

    Returns (centroid, unit_normal). For < 3 points returns the centroid and
    a default normal of (0, 0, 1).
    """
    pts = np.asarray(positions, dtype=float)
    if len(pts) < 3:
        return (pts.mean(axis=0) if len(pts) else np.zeros(3),
                np.array([0.0, 0.0, 1.0]))
    centroid = pts.mean(axis=0)
    centred = pts - centroid
    # SVD: smallest singular vector ≈ plane normal
    _, _, vh = np.linalg.svd(centred, full_matrices=False)
    normal = vh[-1]
    n = normal / np.linalg.norm(normal)
    return centroid, n


def plane_quad(centroid: np.ndarray, normal: np.ndarray,
               half_extent: float = 2.0,
               in_plane_hint: Optional[np.ndarray] = None) -> np.ndarray:
    """Return 4 corners of a square quad lying in the plane (centroid, normal).

    If `in_plane_hint` is provided, the quad's first edge is aligned with the
    component of that vector in-plane (so the quad rotates with the fragment).
    """
    n = normal / np.linalg.norm(normal)
    if in_plane_hint is not None:
        u = np.asarray(in_plane_hint, dtype=float)
        u = u - np.dot(u, n) * n
        if np.linalg.norm(u) < 1e-6:
            u = None
        else:
            u = u / np.linalg.norm(u)
    if in_plane_hint is None or u is None:
        # Pick any vector not parallel to n
        ref = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        u = ref - np.dot(ref, n) * n
        u = u / np.linalg.norm(u)
    v = np.cross(n, u)
    return np.array([
        centroid + half_extent * (u + v),
        centroid + half_extent * (-u + v),
        centroid + half_extent * (-u - v),
        centroid + half_extent * (u - v),
    ])
