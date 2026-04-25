"""Generates minimal Gaussian-format .log files for testing parsers.

Each file has:
  - One Standard orientation block (2 atoms)
  - One Singlet excited state (S1) with orbital transitions
  - One Triplet excited state (T1) with orbital transitions
  - A transition electric dipole moments block with a row for S1
"""
from pathlib import Path
from textwrap import dedent


LOG_TEMPLATE = dedent('''\
 Gaussian 16 test fixture
                         Standard orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1          6           0       0.000000    0.000000    0.000000
      2          1           0       1.000000    0.000000    0.000000
 ---------------------------------------------------------------------
 Some other content between blocks to mimic real logs.

 Excited State   1:      Triplet-A            {t1_energy_ev:.4f} eV  {t1_wavelength_nm:.2f} nm  f=0.0000  <S**2>=2.000
      52 -> 55         {t1_dom:.5f}
      49 -> 56         0.22000

 Excited State   2:      Singlet-A            {s1_energy_ev:.4f} eV  {s1_wavelength_nm:.2f} nm  f=0.0500  <S**2>=0.000
      52 -> 55         {s1_dom:.5f}
      49 -> 56         0.20000
      49 -> 71        -0.10500

 Ground to excited state transition electric dipole moments (Au):
       state          X           Y           Z        Dip. S.      Osc.
         1         0.00000    0.00000    0.00000    0.00000    0.0000
         2         {tdm_x:.5f}    {tdm_y:.5f}    {tdm_z:.5f}    0.14000    0.0500
 End of table.
''')


def write_fixture(path: Path, *, s1_energy_ev: float, t1_energy_ev: float,
                  s1_dom: float = 0.65880, t1_dom: float = 0.64482,
                  tdm: tuple = (0.1, 0.2, 0.3)) -> Path:
    """Write one synthetic Gaussian log file at `path`. Returns the path."""
    path.parent.mkdir(parents=True, exist_ok=True)
    content = LOG_TEMPLATE.format(
        s1_energy_ev=s1_energy_ev,
        s1_wavelength_nm=1239.84 / s1_energy_ev,
        t1_energy_ev=t1_energy_ev,
        t1_wavelength_nm=1239.84 / t1_energy_ev,
        s1_dom=s1_dom,
        t1_dom=t1_dom,
        tdm_x=tdm[0], tdm_y=tdm[1], tdm_z=tdm[2],
    )
    path.write_text(content)
    return path


def write_scan_folder(folder: Path, angle_to_gap: dict) -> Path:
    """Write a set of fixture log files under `folder`.

    `angle_to_gap` maps angle_deg -> (s1_energy_ev, t1_energy_ev).
    """
    folder.mkdir(parents=True, exist_ok=True)
    for angle, (s1e, t1e) in angle_to_gap.items():
        write_fixture(folder / f"angle_{angle}.log", s1_energy_ev=s1e, t1_energy_ev=t1e)
    return folder
