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


def two_p_squared_energy_curve(scan: AngleScan) -> Tuple[List[float], List[float]]:
    """Alternative 2·P² curve where P = S₁ vertical excitation energy in eV."""
    pts = scan.sorted_points
    return [p.angle_deg for p in pts], [p.two_p_squared_energy_ev2 for p in pts]


def tadf_candidates(scan: AngleScan, threshold_ev: float = 0.2) -> List[AnglePoint]:
    return scan.points_below_threshold(threshold_ev)


def summary_row(p: AnglePoint) -> Dict[str, object]:
    s1 = p.s1_state
    t1 = p.t1_state
    dom = p.s1_dominant_transition
    top2 = p.s1_top2_transitions
    second = top2[1] if len(top2) >= 2 else None
    return {
        "angle_deg": p.angle_deg,
        "s1_energy_ev": s1.energy_ev if s1 else None,
        "t1_energy_ev": t1.energy_ev if t1 else None,
        "gap_ev": p.s1_t1_gap_ev,
        "tdm_magnitude": p.s1_tdm_magnitude,
        "two_p_squared": p.two_p_squared,
        "s1_dominant": dom,
        "s1_second": second,
        "source_path": p.source_path,
    }
