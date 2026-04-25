"""Tests for fragment detection and dihedral geometry."""
import numpy as np
import pytest

from tadf_workbench.core import Stage
from tadf_workbench.core.fragments import (
    build_adjacency, split_at_bond, BondNotRotatableError,
    pick_reference_atom, dihedral_angle, fit_plane, plane_quad,
)


def _linear_chain_stage() -> Stage:
    """C-C-C-C-C linear chain — all single bonds, no rings."""
    atoms = [(6, float(i) * 1.5, 0.0, 0.0) for i in range(5)]
    return Stage("linear", atoms)


def _two_rings_via_bridge_stage() -> Stage:
    """Two equilateral C-triangles bridged by a single C–C bond.

    Spacing chosen so no stray cross-bonds form (all non-bonded distances > 2.1 Å).

    Atoms 0,1,2 form ring A (donor-like).
    Atoms 3,4,5 form ring B (acceptor-like).
    Atom 2 is bonded to atom 3 — the rotatable bridge.
    """
    atoms = [
        (6, 0.0,  0.0, 0.0),
        (6, 1.5,  0.0, 0.0),
        (6, 0.75, 1.3, 0.0),     # closes ring A
        (6, 0.75, 2.8, 0.0),     # bonded to atom 2 (1.5 Å away — the bridge)
        (6, 2.25, 2.8, 0.0),
        (6, 1.5,  4.1, 0.0),     # closes ring B
    ]
    return Stage("bridge", atoms)


def test_build_adjacency_finds_chain_bonds():
    stage = _linear_chain_stage()
    adj = build_adjacency(stage)
    assert adj[0] == {1}
    assert adj[2] == {1, 3}
    assert adj[4] == {3}


def test_split_at_bond_breaks_chain_into_two_halves():
    stage = _linear_chain_stage()
    side_a, side_b = split_at_bond(stage, 1, 2)
    assert side_a == {0, 1}
    assert side_b == {2, 3, 4}


def test_split_at_bond_separates_two_rings():
    stage = _two_rings_via_bridge_stage()
    side_a, side_b = split_at_bond(stage, 2, 3)
    assert side_a == {0, 1, 2}
    assert side_b == {3, 4, 5}
    assert side_a.isdisjoint(side_b)


def test_split_at_bond_raises_on_ring_bond():
    stage = _two_rings_via_bridge_stage()
    # The 0-1 bond is part of triangle 0-1-2; cutting it doesn't separate
    with pytest.raises(BondNotRotatableError):
        split_at_bond(stage, 0, 1)


def test_split_at_bond_raises_on_non_bond():
    stage = _linear_chain_stage()
    with pytest.raises(BondNotRotatableError):
        split_at_bond(stage, 0, 4)  # not bonded


def test_split_at_bond_rejects_same_atom():
    stage = _linear_chain_stage()
    with pytest.raises(ValueError):
        split_at_bond(stage, 1, 1)


def test_pick_reference_atom_returns_in_fragment_neighbor():
    stage = _two_rings_via_bridge_stage()
    side_a, _ = split_at_bond(stage, 2, 3)
    ref = pick_reference_atom(stage, side_a, anchor=2)
    assert ref in side_a
    assert ref != 2


def test_pick_reference_atom_handles_isolated_anchor():
    """Anchor with no neighbours in fragment returns None."""
    stage = _linear_chain_stage()
    # Build a 1-atom fragment manually
    ref = pick_reference_atom(stage, fragment={0}, anchor=0)
    assert ref is None


def test_dihedral_planar_zero():
    """All four points coplanar in the XY plane → dihedral 0°."""
    p1 = np.array([0, 1, 0])
    p2 = np.array([0, 0, 0])
    p3 = np.array([1, 0, 0])
    p4 = np.array([1, 1, 0])
    assert dihedral_angle(p1, p2, p3, p4) == pytest.approx(0.0, abs=1e-6)


def test_dihedral_orthogonal_90():
    """p4 lifted into z makes the dihedral 90°."""
    p1 = np.array([0, 1, 0])
    p2 = np.array([0, 0, 0])
    p3 = np.array([1, 0, 0])
    p4 = np.array([1, 0, 1])
    assert abs(dihedral_angle(p1, p2, p3, p4)) == pytest.approx(90.0, abs=1e-6)


def test_dihedral_anti_180():
    """p4 on the opposite side of the bond axis from p1 → ±180°."""
    p1 = np.array([0, 1, 0])
    p2 = np.array([0, 0, 0])
    p3 = np.array([1, 0, 0])
    p4 = np.array([1, -1, 0])
    assert abs(dihedral_angle(p1, p2, p3, p4)) == pytest.approx(180.0, abs=1e-6)


def test_fit_plane_recovers_xy_plane_normal():
    pts = np.array([
        [0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [0.5, 0.5, 0],
    ])
    centroid, normal = fit_plane(pts)
    assert np.allclose(centroid, [0.5, 0.5, 0.0])
    assert abs(abs(normal[2]) - 1.0) < 1e-6  # normal points along ±z


def test_fit_plane_handles_few_points():
    centroid, normal = fit_plane(np.array([[0, 0, 0]]))
    assert np.allclose(centroid, [0, 0, 0])
    assert np.allclose(normal, [0, 0, 1])


def test_plane_quad_lies_in_plane():
    centroid = np.array([0.0, 0.0, 0.0])
    normal = np.array([0.0, 0.0, 1.0])
    quad = plane_quad(centroid, normal, half_extent=1.0)
    # All corners on z=0 plane
    assert np.allclose(quad[:, 2], 0)
    # All at distance √2 from centroid
    for corner in quad:
        assert np.linalg.norm(corner - centroid) == pytest.approx(np.sqrt(2), abs=1e-6)
