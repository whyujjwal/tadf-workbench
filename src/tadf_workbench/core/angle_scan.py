"""Data model for a TADF donor-acceptor dihedral angle scan."""

from dataclasses import dataclass, field
from typing import List, Optional, Set, Tuple

from .molecule import Molecule
from .quantum_data import ExcitedState, TransitionDipole


@dataclass
class AnglePoint:
    """One scanned dihedral angle: molecule geometry + quantum data + derived metrics."""

    angle_deg: float
    source_path: str
    molecule: Molecule

    # Lazily computed — derived from molecule.quantum_data
    @property
    def s1_state(self) -> Optional[ExcitedState]:
        qd = self.molecule.quantum_data
        if not qd:
            return None
        singlets = [s for s in qd.excited_states if s.is_singlet]
        return min(singlets, key=lambda s: s.number) if singlets else None

    @property
    def t1_state(self) -> Optional[ExcitedState]:
        qd = self.molecule.quantum_data
        if not qd:
            return None
        triplets = [s for s in qd.excited_states if s.is_triplet]
        return min(triplets, key=lambda s: s.number) if triplets else None

    @property
    def s1_t1_gap_ev(self) -> Optional[float]:
        s1, t1 = self.s1_state, self.t1_state
        if s1 is None or t1 is None:
            return None
        return s1.energy_ev - t1.energy_ev

    @property
    def s1_dominant_coefficient(self) -> float:
        """Largest-magnitude, positive CI coefficient in S1's orbital-transition list.

        Treats the dominant configuration as the one with the largest |coefficient|;
        returns the coefficient itself (which is normally positive after Gaussian's
        phase convention, but we take abs to be robust).
        """
        s1 = self.s1_state
        if not s1 or not s1.orbital_transitions:
            return 0.0
        best = max(s1.orbital_transitions, key=lambda t: abs(t[2]))
        return abs(best[2])

    @property
    def s1_dominant_transition(self) -> Optional[Tuple[int, int, float]]:
        s1 = self.s1_state
        if not s1 or not s1.orbital_transitions:
            return None
        return max(s1.orbital_transitions, key=lambda t: abs(t[2]))

    @property
    def s1_top2_transitions(self) -> List[Tuple[int, int, float]]:
        """Top 2 orbital transitions for S₁, ranked by |coefficient| desc."""
        s1 = self.s1_state
        if not s1 or not s1.orbital_transitions:
            return []
        ordered = sorted(s1.orbital_transitions, key=lambda t: abs(t[2]), reverse=True)
        return ordered[:2]

    @property
    def fmax_singlet_state(self) -> Optional[ExcitedState]:
        """Singlet (any state number) with the highest oscillator strength.
        Often differs from S₁ — this is the bright state in the Jablonski diagram.
        """
        qd = self.molecule.quantum_data
        if not qd:
            return None
        singlets = [s for s in qd.excited_states if s.is_singlet]
        if not singlets:
            return None
        return max(singlets, key=lambda s: s.oscillator_strength)

    @property
    def fmax_singlet_top2_transitions(self) -> List[Tuple[int, int, float]]:
        s = self.fmax_singlet_state
        if not s or not s.orbital_transitions:
            return []
        ordered = sorted(s.orbital_transitions, key=lambda t: abs(t[2]), reverse=True)
        return ordered[:2]

    @property
    def two_p_squared(self) -> float:
        p = self.s1_dominant_coefficient
        return 2.0 * p * p

    @property
    def two_p_squared_energy_ev2(self) -> float:
        """User-requested non-standard view: P = S₁ energy in eV → 2·P² in eV²."""
        s1 = self.s1_state
        if s1 is None:
            return 0.0
        return 2.0 * s1.energy_ev * s1.energy_ev

    @property
    def s1_tdm(self) -> Optional[TransitionDipole]:
        s1 = self.s1_state
        qd = self.molecule.quantum_data
        if s1 is None or qd is None:
            return None
        return next((td for td in qd.transition_dipoles if td.state == s1.number), None)

    @property
    def s1_tdm_magnitude(self) -> float:
        td = self.s1_tdm
        return td.magnitude if td else 0.0


@dataclass
class AngleScan:
    """Collection of AnglePoints over a dihedral scan."""

    name: str
    points: List[AnglePoint] = field(default_factory=list)
    threshold_ev: float = 0.2

    # Set once by the user after the first atom-pick dialog
    donor_atom_index: Optional[int] = None
    acceptor_atom_index: Optional[int] = None
    rotatable_bond: Optional[Tuple[int, int]] = None  # two atom indices defining axis

    # Fragment membership (set after the rotatable bond is configured).
    # By convention donor_atom_index ∈ donor_fragment, acceptor_atom_index ∈ acceptor_fragment.
    donor_fragment: Optional[Set[int]] = None
    acceptor_fragment: Optional[Set[int]] = None

    # Reference atoms used to compute the IUPAC dihedral N₁-B₁-B₂-N₂.
    # N₁ is a neighbour of donor-side bond atom inside the donor fragment.
    # N₂ is a neighbour of acceptor-side bond atom inside the acceptor fragment.
    donor_reference_atom: Optional[int] = None
    acceptor_reference_atom: Optional[int] = None

    def configure_fragments(self, stage) -> None:
        """Populate fragment + reference-atom fields from the rotatable bond.

        Requires `rotatable_bond`, `donor_atom_index`, `acceptor_atom_index`
        to be set. Uses the supplied Stage's bond connectivity. Raises
        BondNotRotatableError if the bond is in a ring.
        """
        from .fragments import split_at_bond, pick_reference_atom
        if self.rotatable_bond is None:
            raise ValueError("rotatable_bond must be set before configuring fragments")
        b1, b2 = self.rotatable_bond
        side_1, side_2 = split_at_bond(stage, b1, b2)

        # Decide which side is donor based on donor_atom_index. If donor_atom_index
        # isn't set, default to side_1 = donor.
        if self.donor_atom_index is not None and self.donor_atom_index in side_2:
            donor_side, acc_side = side_2, side_1
            donor_anchor, acc_anchor = b2, b1
        else:
            donor_side, acc_side = side_1, side_2
            donor_anchor, acc_anchor = b1, b2

        self.donor_fragment = donor_side
        self.acceptor_fragment = acc_side
        if self.donor_atom_index is None or self.donor_atom_index not in donor_side:
            self.donor_atom_index = donor_anchor
        if self.acceptor_atom_index is None or self.acceptor_atom_index not in acc_side:
            self.acceptor_atom_index = acc_anchor

        self.donor_reference_atom = pick_reference_atom(
            stage, donor_side, donor_anchor)
        self.acceptor_reference_atom = pick_reference_atom(
            stage, acc_side, acc_anchor)

    @property
    def sorted_points(self) -> List[AnglePoint]:
        return sorted(self.points, key=lambda p: p.angle_deg)

    @property
    def angles(self) -> List[float]:
        return [p.angle_deg for p in self.sorted_points]

    @property
    def gaps_ev(self) -> List[float]:
        return [p.s1_t1_gap_ev if p.s1_t1_gap_ev is not None else float("nan")
                for p in self.sorted_points]

    @property
    def tdm_magnitudes(self) -> List[float]:
        return [p.s1_tdm_magnitude for p in self.sorted_points]

    @property
    def two_p_squared_values(self) -> List[float]:
        return [p.two_p_squared for p in self.sorted_points]

    @property
    def two_p_squared_energy_values(self) -> List[float]:
        """θ-sorted list of 2·E(S₁)² in eV² — the alternative-mode 2·P² curve."""
        return [p.two_p_squared_energy_ev2 for p in self.sorted_points]

    def points_below_threshold(self, threshold_ev: Optional[float] = None) -> List[AnglePoint]:
        t = threshold_ev if threshold_ev is not None else self.threshold_ev
        return [p for p in self.sorted_points
                if p.s1_t1_gap_ev is not None and p.s1_t1_gap_ev <= t]

    @property
    def minimum_gap_point(self) -> Optional[AnglePoint]:
        candidates = [p for p in self.points if p.s1_t1_gap_ev is not None]
        if not candidates:
            return None
        return min(candidates, key=lambda p: p.s1_t1_gap_ev)

    def nearest_point(self, angle_deg: float) -> Optional[AnglePoint]:
        if not self.points:
            return None
        return min(self.points, key=lambda p: abs(p.angle_deg - angle_deg))

    def find_point(self, angle_deg: float) -> Optional[AnglePoint]:
        """Exact match (within FP tolerance) for a given angle."""
        for p in self.points:
            if abs(p.angle_deg - angle_deg) < 1e-6:
                return p
        return None

    def add_or_replace(self, point: AnglePoint) -> bool:
        """Add a point. If one already exists at the same angle, replace it.

        Returns True if a previous point was replaced, False if appended fresh.
        """
        existing = self.find_point(point.angle_deg)
        if existing is not None:
            self.points[self.points.index(existing)] = point
            return True
        self.points.append(point)
        return False

    def __len__(self) -> int:
        return len(self.points)
