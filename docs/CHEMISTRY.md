# Chemistry Background

This document explains the science the workbench is built around, so the
metrics on the dashboard mean something concrete to anyone reading the code
or extending it.

## Why TADF, and why a dihedral scan?

**Thermally Activated Delayed Fluorescence** (TADF) lets organic emitters
recover the 75 % triplet population that ordinary fluorophores throw away.
The mechanism is **reverse intersystem crossing** (RISC) from T₁ back up to
S₁, followed by S₁ → S₀ radiative decay. The RISC rate scales as

$$
k_\mathrm{RISC} \propto \exp\!\left(-\frac{\Delta E_\mathrm{ST}}{k_\mathrm{B}T}\right)
$$

so the design goal is to make **ΔE_ST = E(S₁) − E(T₁)** as small as possible
— ideally below the room-temperature thermal energy scale, ~0.025 eV, but in
practice ≤ 0.2 eV is the threshold the field uses to flag a candidate
worth taking seriously.

The classic recipe is a **donor–acceptor (D–A) molecule**: an
electron-donating fragment (carbazole, diphenylamine, phenoxazine, …)
covalently linked to an electron-accepting fragment (triazine,
benzophenone, dicyanobenzene, …) by a single σ bond. The HOMO localises on
the donor, the LUMO on the acceptor. When the two fragments are **twisted
out of coplanarity**, the spatial overlap of HOMO and LUMO drops, the
exchange integral *K* drops with it, and ΔE_ST = 2*K* drops.

The angle that controls all of this is the **dihedral angle θ between the
donor plane and the acceptor plane** about the linker bond. A computational
**rigid scan** evaluates the molecule at θ ∈ {10°, 20°, …, 90°} (or a denser
grid) and records the excited-state energies at each point. This workbench
ingests the resulting log files, locates S₁ and T₁ in each, plots
ΔE_ST(θ), and tells you which angle gives the smallest gap.

## What S₁ and T₁ are in a Gaussian log

Gaussian's TDDFT output emits one block per excited state. The first
**Singlet-***-labelled state in the list is S₁; the first **Triplet-***-
labelled state is T₁. Each block looks like:

```
 Excited State   2:      Singlet-A      3.0300 eV  409.19 nm  f=0.0500  <S**2>=0.000
      52 -> 55         0.65880
      49 -> 56         0.20000
      49 -> 71        -0.10500
```

The energy after `Excited State N:` is in eV, the wavelength in nm, the
oscillator strength in `f=`, and the spin-squared in `<S**2>=`. The
indented lines that follow are the **CI expansion**: each row is one
configuration `i → a` with its expansion coefficient *c*. For closed-shell
references Gaussian normalises so that Σ 2 *c*² = 1 over all configurations
in the state, so each `2·c²` is the **probability weight** of that
configuration.

We pick the configuration with the **largest |c|** as the dominant
single-particle excitation for that state. We call its coefficient **P**;
the workbench reports **2·P²** for S₁ at every θ.

## What the workbench reports per scan point

| Symbol                      | Meaning                                                                                    | Source line in the log                            |
|-----------------------------|--------------------------------------------------------------------------------------------|---------------------------------------------------|
| **θ**                       | Donor–acceptor dihedral angle for this scan point                                           | First integer in filename                         |
| **E(S₁)**                  | Energy of the lowest singlet excited state (eV)                                            | First `Singlet-*` Excited State block             |
| **E(T₁)**                  | Energy of the lowest triplet excited state (eV)                                            | First `Triplet-*` Excited State block             |
| **ΔE(S₁ − T₁)**             | Singlet–triplet gap in eV; the TADF figure of merit. Smaller = better RISC.                 | E(S₁) − E(T₁)                                     |
| **\|TDM(S₁)\|**             | Magnitude of the ground→S₁ transition electric dipole moment (atomic units)                 | Row in *Ground to excited state transition electric dipole moments* |
| **2·P²**                    | Probability weight of the dominant CI configuration for S₁                                  | Largest \|c\| in S₁'s orbital-transition list     |
| **S₁ dominant transition**  | The MO indices `i → a` with the largest CI coefficient (e.g. HOMO→LUMO would be the top row)| Same as above                                     |

## How the IUPAC dihedral is measured

Given the rotatable bond **B₁–B₂** and one chosen reference atom on each
side (**N₁** in the donor fragment, **N₂** in the acceptor fragment), the
**IUPAC dihedral** is the angle between the half-planes (N₁, B₁, B₂) and
(B₁, B₂, N₂), signed by the right-hand rule about B₁→B₂. We compute it as

```
b1 = p_B1 − p_N1
b2 = p_B2 − p_B1
b3 = p_N2 − p_B2
n1 = b1 × b2
n2 = b2 × b3
m1 = n1 × b̂2
θ = atan2( m1·n2 , n1·n2 )      # in (-180°, 180°]
```

implemented in [`core/fragments.py`](../src/tadf_workbench/core/fragments.py)
as `dihedral_angle(p1, p2, p3, p4)`. The result is what the **MEASURED
N₁–B₁–B₂–N₂** label in the 3D viewer shows, independently of the scan-θ
read off the filename.

## How the donor and acceptor are identified

The donor and acceptor fragments are the **two connected components of the
molecular graph after the rotatable bond is removed**. We build the bond
graph from the same distance-based bond perception used elsewhere in the
parser, then BFS from each end of the bond. If the two reachable sets
overlap the bond is part of a ring and is not a valid rotatable bond — the
dialog rejects it with `BondNotRotatableError: bond X-Y is part of a ring`.

This is how the configuration dialog can reduce a 22-atom biphenyl to "11 +
11 atoms, formula C6H5 + C6H5" with one click on the rotatable bond.

## Choice of reference atoms for the dihedral

We pick the **highest-degree heavy-atom neighbour of B₁ that lies inside the
donor fragment** (and analogously for B₂ in the acceptor fragment). This
prefers sp²/aromatic-looking environments and avoids picking a hydrogen,
which would give a noisy dihedral. See `pick_reference_atom` in
[`core/fragments.py`](../src/tadf_workbench/core/fragments.py).

For biphenyl with the bond at C₀–C₆, the references end up as one of the
*ortho* aromatic carbons in each ring — exactly what you'd pick by hand.

## Best-fit fragment planes

For visual context we draw a translucent quad through each fragment's
**best-fit plane**. The plane is found by SVD on the centred atom positions
of the fragment: the smallest singular vector is the plane normal. The
quad's in-plane orientation uses the bond-midpoint-to-fragment-centroid
vector as a hint, so the quad rotates with the fragment as you scrub the
slider. See `fit_plane` and `plane_quad`.

## What "TADF candidate" means in this app

A scan point is flagged as a TADF candidate when **ΔE(S₁ − T₁) ≤ threshold**
(default **0.2 eV**, configurable via **Scan → Set ΔE Threshold…**). Such
points get green highlighting in the table, a green dot in the gap plot, a
green tint on the Current State card's gap banner, and a `✓ YES` in the
TADF stat chip. The minimum-gap point across the whole scan is reported
separately so you can see the best the scan offers regardless of the
threshold.

## What this workbench does NOT do

- It does **not run** the quantum-chemistry calculations. You must produce
  the Gaussian log files yourself (typically a `td=(50-50,nstates=10)`
  TDDFT job at each fixed dihedral angle).
- It does **not** model spin–orbit coupling or compute k_RISC directly.
  ΔE_ST is the proxy.
- It does **not** check that your geometries are local minima or transition
  states. A rigid scan is a coarse exploratory tool — confirm with full
  optimisations before publication.

## Further reading

- Uoyama, H. et al. *Highly efficient organic light-emitting diodes from
  delayed fluorescence.* Nature 492, 234–238 (2012).
- Yang, Z. et al. *Recent advances in organic thermally activated delayed
  fluorescence materials.* Chem. Soc. Rev. 46, 915–1016 (2017).
- Chen, X.-K., Kim, D., Brédas, J.-L. *Thermally Activated Delayed
  Fluorescence (TADF) Path toward Efficient Electroluminescence in Purely
  Organic Materials: Molecular Level Insight.* Acc. Chem. Res. 51, 2215–2224 (2018).
