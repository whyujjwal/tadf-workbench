"""
Parser for Gaussian computational chemistry log files.
"""

from typing import List, Optional
import re
from ..core import Stage, Molecule


class GaussianParser:
    """Parser for Gaussian log files (.log, .out)."""
    
    def __init__(self, filename: str):
        """
        Initialize parser with a filename.
        
        Args:
            filename: Path to Gaussian log file
        """
        self.filename = filename
        self._lines: List[str] = []
        self._molecule: Optional[Molecule] = None
    
    def parse(self) -> Molecule:
        """
        Parse the Gaussian log file and extract all stages.
        
        Returns:
            Molecule object containing all parsed stages and quantum data
        """
        self._read_file()
        self._molecule = Molecule(name=self._extract_molecule_name())
        self._extract_stages()
        self._extract_quantum_data()
        return self._molecule
    
    def _read_file(self):
        """Read the log file into memory."""
        try:
            with open(self.filename, 'r', encoding='utf-8', errors='ignore') as f:
                self._lines = f.readlines()
        except FileNotFoundError:
            raise FileNotFoundError(f"Gaussian log file not found: {self.filename}")
        except Exception as e:
            raise IOError(f"Error reading file {self.filename}: {e}")
    
    def _extract_molecule_name(self) -> str:
        """Extract molecule name from the log file."""
        # Look for title card
        for i, line in enumerate(self._lines):
            if '-----' in line and i + 1 < len(self._lines):
                next_line = self._lines[i + 1].strip()
                # Check if it looks like a title (not a separator)
                if next_line and '-----' not in next_line and len(next_line) < 100:
                    # Filter out warnings and system messages
                    if any(x in next_line for x in ["Warning", "Error", "Gaussian", "Link 1", "copyright"]):
                         continue
                    return next_line
        
        # Fallback to filename
        import os
        return os.path.splitext(os.path.basename(self.filename))[0]
    
    def _extract_stages(self):
        """Extract all coordinate stages from the log file."""
        i = 0
        input_counter = 0
        standard_counter = 0
        
        while i < len(self._lines):
            line = self._lines[i]
            
            # Check for Input orientation (flexible matching)
            if 'Input orientation' in line:
                atoms = self._parse_orientation_block(i)
                if atoms:
                    input_counter += 1
                    stage_name = f"Input Orientation {input_counter}"
                    self._molecule.add_stage_from_data(stage_name, atoms)
                    print(f"Found {stage_name} with {len(atoms)} atoms")
            
            # Check for Standard orientation (flexible matching)
            elif 'Standard orientation' in line:
                atoms = self._parse_orientation_block(i)
                if atoms:
                    standard_counter += 1
                    stage_name = f"Standard Orientation {standard_counter}"
                    self._molecule.add_stage_from_data(stage_name, atoms)
                    print(f"Found {stage_name} with {len(atoms)} atoms")
            
            # Check for Z-Matrix orientation
            elif 'Z-Matrix orientation' in line:
                atoms = self._parse_orientation_block(i)
                if atoms:
                    stage_name = f"Z-Matrix Orientation"
                    self._molecule.add_stage_from_data(stage_name, atoms)
                    print(f"Found {stage_name} with {len(atoms)} atoms")
            
            i += 1
    
    def _parse_orientation_block(self, start_index: int) -> List[tuple]:
        """
        Parse a coordinate orientation block.
        
        Args:
            start_index: Line index where orientation block starts
            
        Returns:
            List of (atomic_number, x, y, z) tuples
        """
        atoms = []
        
        # Find the first separator line (dashes)
        i = start_index + 1
        while i < len(self._lines) and i < start_index + 10:
            if '---' in self._lines[i]:
                # Found first separator, skip to next line (starts header or atoms)
                i += 1
                # Skip header lines (Center, Number, etc.)
                while i < len(self._lines) and 'Center' in self._lines[i] or 'Number' in self._lines[i]:
                    i += 1
                # Skip second separator if present
                if i < len(self._lines) and '---' in self._lines[i]:
                    i += 1
                break
            i += 1
        
        # Now parse atom coordinates
        while i < len(self._lines):
            line = self._lines[i].strip()
            
            # End of coordinate block (separator line or empty line)
            if '---' in line or not line:
                break
            
            # Skip lines that don't start with a number (header remnants)
            parts = line.split()
            if len(parts) < 6:
                i += 1
                continue
            
            # Parse coordinate line - must have numeric first column
            try:
                center_num = int(parts[0])
                atomic_number = int(parts[1])
                atomic_type = int(parts[2])
                x = float(parts[3])
                y = float(parts[4])
                z = float(parts[5])
                
                atoms.append((atomic_number, x, y, z))
            except (ValueError, IndexError):
                # Skip non-coordinate lines
                pass
            
            i += 1
        
        return atoms
    
    def _extract_quantum_data(self):
        """Extract all quantum chemistry data from the log file."""
        from ..core import QuantumData, DipoleMoment, QuadrupoleMoment, ExcitedState, MullikenCharges, TransitionDipole

        quantum_data = QuantumData()

        # Extract dipole moment
        dipole = self._extract_dipole_moment()
        if dipole:
            quantum_data.dipole_moment = dipole
            print(f"Found dipole moment: {dipole.total:.4f} D")

        # Extract quadrupole moment
        quadrupole = self._extract_quadrupole_moment()
        if quadrupole:
            quantum_data.quadrupole_moment = quadrupole

        # Extract excited states
        excited_states = self._extract_excited_states()
        for state in excited_states:
            quantum_data.add_excited_state(state)
        if excited_states:
            print(f"Found {len(excited_states)} excited states")

        # Extract transition dipole moments
        transition_dipoles = self._extract_transition_dipoles()
        for td in transition_dipoles:
            quantum_data.add_transition_dipole(td)
        if transition_dipoles:
            print(f"Found {len(transition_dipoles)} transition dipole moments")

        # Extract Mulliken charges
        charges = self._extract_mulliken_charges()
        if charges:
            quantum_data.mulliken_charges = charges
            print(f"Found Mulliken charges for {len(charges.charges)} atoms")

        self._molecule.quantum_data = quantum_data
    
    def _extract_dipole_moment(self):
        """Extract ground state dipole moment."""
        from ..core import DipoleMoment
        
        for i, line in enumerate(self._lines):
            if 'Dipole moment' in line and 'Debye' in line:
                # Next line should have X= Y= Z= Tot=
                if i + 1 < len(self._lines):
                    data_line = self._lines[i + 1]
                    try:
                        # Parse: X=  value   Y=  value   Z=  value  Tot= value
                        parts = data_line.split()
                        x_idx = parts.index('X=') + 1 if 'X=' in parts else -1
                        y_idx = parts.index('Y=') + 1 if 'Y=' in parts else -1
                        z_idx = parts.index('Z=') + 1 if 'Z=' in parts else -1
                        tot_idx = parts.index('Tot=') + 1 if 'Tot=' in parts else -1
                        
                        if all(idx > 0 for idx in [x_idx, y_idx, z_idx, tot_idx]):
                            x = float(parts[x_idx])
                            y = float(parts[y_idx])
                            z = float(parts[z_idx])
                            tot = float(parts[tot_idx])
                            return DipoleMoment(x, y, z, tot)
                    except (ValueError, IndexError):
                        pass
        return None

    def _extract_transition_dipoles(self):
        """Extract ground to excited state transition electric dipole moments,
        velocity dipole moments, magnetic dipole moments, and rotatory strengths."""
        from ..core import TransitionDipole

        transition_dipoles = []
        vel_data = {}  # state -> (x, y, z, dip_s, osc)
        mag_data = {}  # state -> (x, y, z)
        rot_vel_data = {}  # state -> R(velocity)
        rot_len_data = {}  # state -> R(length)

        for i, line in enumerate(self._lines):
            # Electric dipole moments
            if 'Ground to excited state transition electric dipole moments (Au):' in line:
                j = i + 2
                while j < len(self._lines):
                    data_line = self._lines[j].strip()
                    if not data_line or 'Ground to excited state transition' in data_line:
                        break
                    try:
                        parts = data_line.split()
                        if len(parts) >= 6:
                            state = int(parts[0])
                            td = TransitionDipole(
                                state=state,
                                x=float(parts[1]),
                                y=float(parts[2]),
                                z=float(parts[3]),
                                dip_strength=float(parts[4]),
                                osc_strength=float(parts[5])
                            )
                            transition_dipoles.append(td)
                    except (ValueError, IndexError):
                        pass
                    j += 1

            # Velocity dipole moments
            elif 'Ground to excited state transition velocity dipole moments (Au):' in line:
                j = i + 2
                while j < len(self._lines):
                    data_line = self._lines[j].strip()
                    if not data_line or 'Ground to excited state transition' in data_line:
                        break
                    try:
                        parts = data_line.split()
                        if len(parts) >= 6:
                            vel_data[int(parts[0])] = (
                                float(parts[1]), float(parts[2]), float(parts[3]),
                                float(parts[4]), float(parts[5])
                            )
                    except (ValueError, IndexError):
                        pass
                    j += 1

            # Magnetic dipole moments
            elif 'Ground to excited state transition magnetic dipole moments (Au):' in line:
                j = i + 2
                while j < len(self._lines):
                    data_line = self._lines[j].strip()
                    if not data_line or 'Ground to excited state transition' in data_line:
                        break
                    try:
                        parts = data_line.split()
                        if len(parts) >= 4:
                            mag_data[int(parts[0])] = (
                                float(parts[1]), float(parts[2]), float(parts[3])
                            )
                    except (ValueError, IndexError):
                        pass
                    j += 1

            # Rotatory strengths (velocity gauge)
            elif 'Rotatory Strengths (R) in cgs' in line and 'R(velocity)' in self._lines[i + 1] if i + 1 < len(self._lines) else False:
                j = i + 2
                while j < len(self._lines):
                    data_line = self._lines[j].strip()
                    if not data_line or '1/2[' in data_line or 'Rotatory Strengths' in data_line:
                        break
                    try:
                        parts = data_line.split()
                        if len(parts) >= 5:
                            rot_vel_data[int(parts[0])] = float(parts[4])
                    except (ValueError, IndexError):
                        pass
                    j += 1

            # Rotatory strengths (length gauge)
            elif 'Rotatory Strengths (R) in cgs' in line and 'R(length)' in self._lines[i + 1] if i + 1 < len(self._lines) else False:
                j = i + 2
                while j < len(self._lines):
                    data_line = self._lines[j].strip()
                    if not data_line or 'Excited State' in data_line:
                        break
                    try:
                        parts = data_line.split()
                        if len(parts) >= 5:
                            rot_len_data[int(parts[0])] = float(parts[4])
                    except (ValueError, IndexError):
                        pass
                    j += 1

        # Merge velocity, magnetic, rotatory data into TransitionDipole objects
        for td in transition_dipoles:
            if td.state in vel_data:
                vx, vy, vz, vds, vosc = vel_data[td.state]
                td.vel_x, td.vel_y, td.vel_z = vx, vy, vz
                td.vel_dip_strength, td.vel_osc_strength = vds, vosc
            if td.state in mag_data:
                td.mag_x, td.mag_y, td.mag_z = mag_data[td.state]
            if td.state in rot_vel_data:
                td.rotatory_velocity = rot_vel_data[td.state]
            if td.state in rot_len_data:
                td.rotatory_length = rot_len_data[td.state]

        return transition_dipoles

    def _extract_quadrupole_moment(self):
        """Extract ground state quadrupole moment."""
        from ..core import QuadrupoleMoment
        
        for i, line in enumerate(self._lines):
            if 'Traceless Quadrupole moment' in line and 'Debye-Ang' in line:
                try:
                    # Next two lines have XX= YY= ZZ= and XY= XZ= YZ=
                    if i + 2 < len(self._lines):
                        line1 = self._lines[i + 1].split()
                        line2 = self._lines[i + 2].split()
                        
                        xx = float(line1[line1.index('XX=') + 1])
                        yy = float(line1[line1.index('YY=') + 1])
                        zz = float(line1[line1.index('ZZ=') + 1])
                        xy = float(line2[line2.index('XY=') + 1])
                        xz = float(line2[line2.index('XZ=') + 1])
                        yz = float(line2[line2.index('YZ=') + 1])
                        
                        return QuadrupoleMoment(xx, yy, zz, xy, xz, yz)
                except (ValueError, IndexError):
                    pass
        return None
    
    def _extract_excited_states(self):
        """Extract all excited states."""
        from ..core import ExcitedState
        
        states = []
        i = 0
        
        while i < len(self._lines):
            line = self._lines[i]
            
            # Look for "Excited State" lines
            if 'Excited State' in line:
                try:
                    # Parse line: Excited State   5:      Singlet-A      3.7489 eV  330.72 nm  f=0.1609  <S**2>=0.000
                    parts = line.split()
                    
                    # Find key indices
                    state_idx = parts.index('State') + 1
                    number = int(parts[state_idx].rstrip(':'))
                    
                    # Multiplicity and symmetry
                    mult_sym = parts[state_idx + 1].split('-')
                    multiplicity = mult_sym[0]
                    symmetry = mult_sym[1] if len(mult_sym) > 1 else 'A'
                    
                    # Energy
                    ev_idx = [j for j, p in enumerate(parts) if p == 'eV'][0]
                    energy_ev = float(parts[ev_idx - 1])
                    
                    # Wavelength
                    nm_idx = [j for j, p in enumerate(parts) if p == 'nm'][0]
                    wavelength_nm = float(parts[nm_idx - 1])
                    
                    # Oscillator strength
                    f_str = [p for p in parts if p.startswith('f=')][0]
                    osc_strength = float(f_str.split('=')[1])
                    
                    # S squared
                    s2_str = [p for p in parts if '<S**2>' in p or 'S**2' in p][0]
                    s_squared = float(s2_str.split('=')[1])
                    
                    # Parse orbital transitions (next few lines)
                    orbital_transitions = []
                    j = i + 1
                    while j < len(self._lines) and j < i + 10:
                        trans_line = self._lines[j].strip()
                        if not trans_line or 'Excited State' in trans_line:
                            break
                        # Parse: 50 -> 54        -0.10260
                        if '->' in trans_line:
                            trans_parts = trans_line.split()
                            try:
                                from_orb = int(trans_parts[0])
                                to_orb = int(trans_parts[2])
                                coeff = float(trans_parts[3])
                                orbital_transitions.append((from_orb, to_orb, coeff))
                            except (ValueError, IndexError):
                                pass
                        j += 1
                    
                    state = ExcitedState(
                        number=number,
                        multiplicity=multiplicity,
                        symmetry=symmetry,
                        energy_ev=energy_ev,
                        wavelength_nm=wavelength_nm,
                        oscillator_strength=osc_strength,
                        s_squared=s_squared,
                        orbital_transitions=orbital_transitions
                    )
                    states.append(state)
                    
                except (ValueError, IndexError) as e:
                    print(f"Warning: Could not parse excited state: {e}")
            
            i += 1
        
        return states
    
    def _extract_mulliken_charges(self):
        """Extract Mulliken atomic charges."""
        from ..core import MullikenCharges
        
        for i, line in enumerate(self._lines):
            if 'Mulliken charges:' in line:
                charges = []
                indices = []
                
                j = i + 2  # Skip header line
                while j < len(self._lines):
                    charge_line = self._lines[j].strip()
                    if not charge_line or 'Sum of Mulliken' in charge_line:
                        break
                    
                    try:
                        parts = charge_line.split()
                        if len(parts) >= 3:
                            atom_idx = int(parts[0])
                            charge = float(parts[2])
                            indices.append(atom_idx)
                            charges.append(charge)
                    except (ValueError, IndexError):
                        pass
                    
                    j += 1
                
                if charges:
                    return MullikenCharges(indices, charges)
        
        return None
    
    @staticmethod
    def is_gaussian_file(filename: str) -> bool:
        """
        Check if a file appears to be a Gaussian log file.
        
        Args:
            filename: Path to file
            
        Returns:
            True if file appears to be a Gaussian log
        """
        try:
            with open(filename, 'r', encoding='utf-8', errors='ignore') as f:
                # Read first few lines
                header = ''.join(f.readlines()[:20])
                
                # Check for Gaussian markers
                gaussian_markers = [
                    'Gaussian',
                    'G16',
                    'G09',
                    'G03',
                    'orientation:',
                    'Entering Link 1',
                ]
                
                return any(marker in header for marker in gaussian_markers)
        except Exception:
            return False


def parse_gaussian_file(filename: str) -> Molecule:
    """
    Convenience function to parse a Gaussian log file.
    
    Args:
        filename: Path to Gaussian log file
        
    Returns:
        Molecule object with parsed stages
    """
    parser = GaussianParser(filename)
    return parser.parse()
