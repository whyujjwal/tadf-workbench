# tests/test_angle_scan_core.py
import pytest
from tadf_workbench.core import Molecule, AnglePoint, AngleScan
from tadf_workbench.core.quantum_data import (
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
        _make_point(angle=90, s1_e=3.15, t1_e=2.95),  # gap ≈ 0.2 (just under, due to FP)
        _make_point(angle=30, s1_e=3.0, t1_e=2.95),  # gap 0.05
        _make_point(angle=60, s1_e=3.4, t1_e=3.0),  # gap 0.4
    ])
    assert [p.angle_deg for p in scan.sorted_points] == [30, 60, 90]
    tadf = scan.points_below_threshold(0.2)
    assert [p.angle_deg for p in tadf] == [30, 90]  # ≤ threshold
    best = scan.minimum_gap_point
    assert best.angle_deg == 30


def test_anglescan_add_or_replace_appends_new_angle():
    scan = AngleScan(name="demo", points=[_make_point(angle=10, s1_e=3.4, t1_e=2.9)])
    new_pt = _make_point(angle=20, s1_e=3.3, t1_e=2.95)
    replaced = scan.add_or_replace(new_pt)
    assert replaced is False
    assert len(scan) == 2
    assert scan.angles == [10.0, 20.0]


def test_anglescan_add_or_replace_overwrites_existing_angle():
    original = _make_point(angle=30, s1_e=3.4, t1_e=2.9)
    scan = AngleScan(name="demo", points=[original])
    updated = _make_point(angle=30, s1_e=3.0, t1_e=2.95)  # different energies
    replaced = scan.add_or_replace(updated)
    assert replaced is True
    assert len(scan) == 1
    assert scan.find_point(30).s1_state.energy_ev == pytest.approx(3.0)


def test_anglescan_find_point_returns_none_when_missing():
    scan = AngleScan(name="demo", points=[_make_point(angle=10, s1_e=3.4, t1_e=2.9)])
    assert scan.find_point(99) is None


def test_anglescan_donor_acceptor_defaults_none():
    scan = AngleScan(name="demo", points=[])
    assert scan.donor_atom_index is None
    assert scan.acceptor_atom_index is None
    assert scan.rotatable_bond is None
