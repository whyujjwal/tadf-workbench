"""Pure-function layout for a per-θ Jablonski energy-level diagram.

Computes positions, labels, and arrows; contains no Qt code so it can be
unit-tested and reused by both the live single-angle panel and the
all-angles grid window.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Literal, Optional, Tuple

from ..core import AnglePoint, ExcitedState

Mode = Literal["key", "all"]


def percent_2p2(coeff: float) -> float:
    """2·P² × 100 (percent display of the dominant-configuration weight)."""
    return 2.0 * coeff * coeff * 100.0


def format_jump(t: Tuple[int, int, float]) -> str:
    """'96→101 65.4%' for a single (from_orb, to_orb, coeff) tuple."""
    i, a, c = t
    return f"{i}→{a} {percent_2p2(c):.1f}%"


def format_top2(top2: List[Tuple[int, int, float]]) -> str:
    """'96→101 65.4% / 92→103 5.8%' (or just the first if there's only one)."""
    if not top2:
        return ""
    parts = [format_jump(t) for t in top2[:2]]
    return " / ".join(parts)


@dataclass(frozen=True)
class Level:
    """One horizontal energy level in the diagram."""
    energy_ev: float
    multiplicity: str          # "Singlet" or "Triplet"
    state_number: int          # 0 = ground, otherwise the state #
    label: str                 # full multi-line label, "" when not highlighted
    is_highlighted: bool
    column: str                # "singlet" or "triplet"
    role: str = ""             # "S0" | "S1" | "Sn" | "T1" | "" (for plain levels)


@dataclass(frozen=True)
class Arrow:
    """One arrow on the diagram (absorption, ISC, etc.)."""
    from_energy_ev: float
    to_energy_ev: float
    kind: str                  # "absorption" | "isc"
    label: str
    column: str                # "singlet" | "triplet" | "between"


@dataclass(frozen=True)
class JablonskiLayout:
    levels: List[Level]
    arrows: List[Arrow]
    title: str
    angle_deg: float
    mode: Mode
    has_data: bool = True

    @property
    def max_energy_ev(self) -> float:
        if not self.levels:
            return 1.0
        return max(lv.energy_ev for lv in self.levels)


_S0 = Level(
    energy_ev=0.0, multiplicity="Singlet", state_number=0,
    label="S₀", is_highlighted=True, column="singlet", role="S0",
)


def _label_for_role(role: str, state: ExcitedState,
                    top2: List[Tuple[int, int, float]]) -> str:
    """Compact label that goes next to a highlighted level.

    Mirrors the user's whiteboard sketch: heading + dominant orbital jump
    with 2·P²%, optionally a second jump on its own line. Energy lives on
    the y-axis already, and `f` is on the absorption arrow, so we don't
    repeat them.
    """
    head = {
        "S1": "S₁",
        "Sn": "Sₙ",
        "T1": "T₁",
    }.get(role, f"{state.multiplicity[0]}{state.number}")

    if not top2:
        return head

    first = format_jump(top2[0])
    if len(top2) == 1:
        return f"{head}\n{first}"
    second = format_jump(top2[1])
    return f"{head}\n{first}\n{second}"


def empty_layout(angle_deg: float, mode: Mode = "key",
                 reason: str = "no data") -> JablonskiLayout:
    return JablonskiLayout(
        levels=[], arrows=[], title=f"θ = {angle_deg:.0f}° — {reason}",
        angle_deg=angle_deg, mode=mode, has_data=False,
    )


def jablonski_layout(point: Optional[AnglePoint], mode: Mode = "key") -> JablonskiLayout:
    """Compute the diagram for one AnglePoint.

    `mode="key"` emits S₀, S₁, Sₙ (f-max singlet), T₁ only.
    `mode="all"` adds every other Singlet/Triplet as faint, unlabelled bars.
    """
    if point is None:
        return empty_layout(0.0, mode, "no point selected")

    s1 = point.s1_state
    t1 = point.t1_state
    sn = point.fmax_singlet_state  # f-max singlet — may equal s1

    if s1 is None and t1 is None:
        return empty_layout(point.angle_deg, mode, "no excited states")

    levels: List[Level] = [_S0]
    arrows: List[Arrow] = []

    # S₁
    if s1 is not None:
        s1_top2 = point.s1_top2_transitions
        levels.append(Level(
            energy_ev=s1.energy_ev, multiplicity="Singlet",
            state_number=s1.number, column="singlet",
            label=_label_for_role("S1", s1, s1_top2),
            is_highlighted=True, role="S1",
        ))
        arrows.append(Arrow(
            from_energy_ev=0.0, to_energy_ev=s1.energy_ev,
            kind="absorption", column="singlet",
            label=f"f = {s1.oscillator_strength:.4f}",
        ))

    # Sₙ — only emit as a separate level if it differs from S₁
    if sn is not None and (s1 is None or sn.number != s1.number):
        sn_top2 = point.fmax_singlet_top2_transitions
        levels.append(Level(
            energy_ev=sn.energy_ev, multiplicity="Singlet",
            state_number=sn.number, column="singlet",
            label=_label_for_role("Sn", sn, sn_top2),
            is_highlighted=True, role="Sn",
        ))
        arrows.append(Arrow(
            from_energy_ev=0.0, to_energy_ev=sn.energy_ev,
            kind="absorption", column="singlet",
            label=f"f = {sn.oscillator_strength:.4f}  (max)",
        ))

    # T₁
    if t1 is not None:
        t1_top2 = sorted(
            t1.orbital_transitions or [], key=lambda x: abs(x[2]), reverse=True
        )[:2]
        levels.append(Level(
            energy_ev=t1.energy_ev, multiplicity="Triplet",
            state_number=t1.number, column="triplet",
            label=_label_for_role("T1", t1, t1_top2),
            is_highlighted=True, role="T1",
        ))

    # ISC arrow S₁ ↔ T₁ labelled with ΔE_ST
    if s1 is not None and t1 is not None and point.s1_t1_gap_ev is not None:
        arrows.append(Arrow(
            from_energy_ev=s1.energy_ev, to_energy_ev=t1.energy_ev,
            kind="isc", column="between",
            label=f"ΔE_ST = {point.s1_t1_gap_ev:+.3f} eV",
        ))

    # In "all" mode add every other singlet/triplet as a faint, unlabelled level.
    if mode == "all":
        qd = point.molecule.quantum_data
        if qd is not None:
            highlighted = {(lv.multiplicity, lv.state_number) for lv in levels}
            for st in qd.excited_states:
                key = (st.multiplicity, st.number)
                if key in highlighted:
                    continue
                levels.append(Level(
                    energy_ev=st.energy_ev,
                    multiplicity=st.multiplicity,
                    state_number=st.number,
                    column="singlet" if st.is_singlet else "triplet",
                    label="", is_highlighted=False,
                ))

    title = f"θ = {point.angle_deg:.0f}°"
    if point.s1_t1_gap_ev is not None:
        title += f"   ΔE = {point.s1_t1_gap_ev:+.3f} eV"
    if sn is not None:
        title += f"   f_max = {sn.oscillator_strength:.3f}"

    return JablonskiLayout(
        levels=levels, arrows=arrows, title=title,
        angle_deg=point.angle_deg, mode=mode,
    )
