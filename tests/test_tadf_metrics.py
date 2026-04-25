# tests/test_tadf_metrics.py
import pytest
from tadf_workbench.analysis.tadf_metrics import (
    gap_curve, tdm_curve, two_p_squared_curve, tadf_candidates, summary_row,
)
from tests.fixtures.build_fixtures import write_scan_folder
from tadf_workbench.parsers import parse_angle_folder


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
