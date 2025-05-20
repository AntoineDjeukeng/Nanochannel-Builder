import numpy as np
import re
from collections import defaultdict

from md_simulation.io.gro_file import GroFile

from md_simulation.io.lmp_line import lmpline
import copy
class LmpFile:
    """
    Represents a LAMMPS data file.

    Attributes:
        atoms (list of dict): List of atoms.
        atom_types (list): List of atom types.
        box (np.ndarray): Box dimensions, shape (3, 2) for x, y, z.
        bonds (list of tuple): List of bonds (id, type, atom1, atom2).
        angles (list of tuple): List of angles (id, type, atom1, atom2, atom3).
        top_system_data (dict): Info like molecule types and counts.
    """

    def __init__(self):
        self.title =""
        self.atoms = []
        self.natoms = 0
        self.atom_types = {}
        self.box = np.zeros((3, 2))
        self.gro_box = []
        self.bonds = []
        self.angles = []
        self.top_system_data = {}

    @classmethod
    def from_gro(cls, gro_obj):
        """
        Build a LmpFile from a GroFile instance.

        Parameters:
            gro_obj (GroFile): Parsed GroFile instance.

        Returns:
            LmpFile instance
        """
        obj = cls()
        obj.atoms = gro_obj.atoms
        obj.atom_types = gro_obj.atom_types
        obj.natoms = gro_obj.natoms
        obj.box = gro_obj.lmp_box
        obj.bonds = gro_obj.bonds
        obj.angles = gro_obj.angles
        obj.top_system_data = gro_obj.top_system_data
        return obj

    def write(self, filename):
        """
        Write the LAMMPS data file in full atom style,
        grouping atoms by molecule and remapping indices for bonds/angles.
        """
        with open(filename, 'w') as f:
            f.write("LAMMPS data file via LmpFile\n\n")
            f.write(f"{self.natoms} atoms\n")
            f.write(f"{len(self.atom_types)} atom types\n")
            f.write(f"{len(self.bonds)} bonds\n")
            f.write(f"1 bond types\n")
            f.write(f"{len(self.angles)} angles\n")
            f.write(f"1 angle types\n\n")

            f.write(f"{10*self.box[0][0]:.4f} {10*self.box[0][1]:.4f} xlo xhi\n")
            f.write(f"{10*self.box[1][0]:.4f} {10*self.box[1][1]:.4f} ylo yhi\n")
            f.write(f"{10*self.box[2][0]:.4f} {10*self.box[2][1]:.4f} zlo zhi\n\n")

            # # Prepare index remapping
            next_atom_id = 1
            next_mol_id = 1

            f.write("Atoms\n\n")
            mol_data = {}
            old_data = {}
            for atom in self.atoms:
                mol_name = atom['mol_name']
                mol_index = atom['mol_index']
                if mol_name not in mol_data:
                    mol_data[mol_name] = {}
                if mol_index not in mol_data[mol_name]:
                    mol_data[mol_name][mol_index] = []
                
                mol_data[mol_name][mol_index].append(atom)
            

            for mol_name, mol_indices in mol_data.items():
                for mol_index, mol_atoms in mol_indices.items():
                    for atom in mol_atoms:
                        atom_type = self.atom_types.get(atom['atom_name'], 1)
                        
                        charge = 0.0
                        if 'atom_charge' in atom:
                            charge = atom['atom_charge']

                        data = {
                            "atom_index": next_atom_id,
                            "mol_index": next_mol_id,
                            "atom_type": atom_type,
                            "atom_charge": charge,
                            "x": 10 * atom["x"],
                            "y": 10 * atom["y"],
                            "z": 10 * atom["z"]
                        }
                        line = lmpline(data, 'write')
                        f.write(line)
                        old_data[atom['atom_index']] = next_atom_id
                        next_atom_id += 1
                    next_mol_id += 1
            
            # Write bonds
            if self.bonds:
                f.write("\nBonds\n\n")
                n=0
                for bond in self.bonds:
                    n+=1
                    f.write(f"{n} 1  {old_data[bond['i']]} {old_data[bond['j']]}\n")
            if self.angles:
                f.write("\nAngles\n\n")
                n=0
                for angle in self.angles:
                    n+=1
                    f.write(f"{n} 1 {old_data[angle['i']]} {old_data[angle['j']]} {old_data[angle['k']]}\n")


    @classmethod
    def read(cls, filename):
        """
        Read a LAMMPS data file written by this class and return a LmpFile instance.

        This method assumes:
        - The header section lists counts of atoms, atom types, bonds, bond types, angles, angle types.
        - Box bounds are in the format "xlo xhi", "ylo yhi", "zlo zhi".
        - Sections "Atoms", "Bonds", and "Angles" appear in order.
        - Atoms lines: id mol_id atom_type charge x y z
        - Bonds lines: id bond_type atom1 atom2
        - Angles lines: id angle_type atom1 atom2 atom3

        Returns:
            LmpFile instance with populated attributes.
        """
        obj = cls()
        
        with open(filename, 'r') as f:
            lines = f.readlines()
        obj.title = lines[0].strip()  # First line is title
        # Regex helpers for numeric lines
        natoms = 0
        nbond = 0
        nangle = 0
        atom_types_count = 0
        bond_types_count = 0
        angle_types_count = 0

        # Parse header
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if line.endswith("atoms"):
                natoms = int(line.split()[0])
            elif line.endswith("atom types"):
                atom_types_count = int(line.split()[0])
            elif line.endswith("bonds"):
                nbond = int(line.split()[0])
            elif line.endswith("bond types"):
                bond_types_count = int(line.split()[0])
            elif line.endswith("angles"):
                nangle = int(line.split()[0])
            elif line.endswith("angle types"):
                angle_types_count = int(line.split()[0])
            elif re.match(r"^[\d\.\-eE\+]+\s+[\d\.\-eE\+]+\s+xlo xhi", line):
                vals = line.split()
                obj.box[0, 0] = float(vals[0]) / 10.0  # Convert back to nm if needed
                obj.box[0, 1] = float(vals[1]) / 10.0
            elif re.match(r"^[\d\.\-eE\+]+\s+[\d\.\-eE\+]+\s+ylo yhi", line):
                vals = line.split()
                obj.box[1, 0] = float(vals[0]) / 10.0
                obj.box[1, 1] = float(vals[1]) / 10.0
            elif re.match(r"^[\d\.\-eE\+]+\s+[\d\.\-eE\+]+\s+zlo zhi", line):
                vals = line.split()
                obj.box[2, 0] = float(vals[0]) / 10.0
                obj.box[2, 1] = float(vals[1]) / 10.0

            elif line == "Atoms":
                i += 2  # Skip the blank line after "Atoms"
                # Read atoms
                atoms = []
                for _ in range(natoms):
                    # parts = lines[i].strip().split()
                    # # Expected format: id mol_id atom_type charge x y z
                    # atom_id = int(parts[0])
                    # mol_id = int(parts[1])
                    # atom_type = int(parts[2])
                    # charge = float(parts[3])  # can ignore if not used
                    # x = float(parts[4]) / 10.0  # Convert back from Angstrom to nm
                    # y = float(parts[5]) / 10.0
                    # z = float(parts[6]) / 10.0

                    # atom = {
                    #     'atom_index': atom_id,  # original atom id
                    #     'mol_index': mol_id,
                    #     'atom_type': atom_type,
                    #     'charge': charge,
                    #     'x': x,
                    #     'y': y,
                    #     'z': z,
                    #     # We'll fill mol_name and atom_name later or guess dummy for now
                    #     'mol_name': f'mol{mol_id}',
                    #     'atom_name': f'at{atom_type}',
                    # }
                    atom = lmpline(lines[i], 'read')
                    atoms.append(atom)
                    i += 1

                obj.atoms = atoms
                obj.natoms = natoms

            elif line == "Bonds":
                i += 2  # Skip the blank line after "Bonds"
                bonds = []
                for _ in range(nbond):
                    parts = lines[i].strip().split()
                    # Expected format: id bond_type atom1 atom2
                    bond_id = int(parts[0])
                    bond_type = int(parts[1])
                    a1 = int(parts[2])
                    a2 = int(parts[3])
                    bonds.append({'id': bond_id, 'type': bond_type, 'i': a1, 'j': a2})
                    i += 1
                obj.bonds = bonds

            elif line == "Angles":
                i += 2  # Skip the blank line after "Angles"
                angles = []
                for _ in range(nangle):
                    parts = lines[i].strip().split()
                    # Expected format: id angle_type atom1 atom2 atom3
                    angle_id = int(parts[0])
                    angle_type = int(parts[1])
                    a1 = int(parts[2])
                    a2 = int(parts[3])
                    a3 = int(parts[4])
                    angles.append({'id': angle_id, 'type': angle_type, 'i': a1, 'j': a2, 'k': a3})
                    i += 1
                obj.angles = angles

            i += 1
        a=obj.box[0][1] - obj.box[0][0]
        b=obj.box[1][1] - obj.box[1][0]
        c=obj.box[2][1] - obj.box[2][0]
        # obj.gro_box = [str(a:.5f), str(b:.5f), str(c:.5f)]
        obj.gro_box=[a,b,c]
        # Build atom_types dict from atoms (assuming atom_type maps to atom_name)
        atom_types = {}
        for atom in obj.atoms:
            # Use atom_name from atom_type or fallback dummy
            atype = atom.get('atom_type', 1)
            if atype not in atom_types:
                atom_types[atype] = atom.get('atom_name', f'at{atype}')
        obj.atom_types = atom_types

        # Optionally fill top_system_data based on molecule IDs present
        mol_counts = {}
        for atom in obj.atoms:
            mol_name = atom.get('mol_name', f"mol{atom['mol_index']}")
            mol_counts[mol_name] = mol_counts.get(mol_name, 0) + 1
        obj.top_system_data = mol_counts

        return obj

    def to_gro(self, mol_atom_map: dict) -> GroFile:
        """
        Convert a LAMMPS data structure into a GroFile using a mapping from
        molecule names to atom names and atom types.

        Args:
            mol_atom_map (dict): Mapping like
                {
                    "SOL": {"OW": 1, "HW1": 2, "HW2": 3},
                    "Na": {"NA": 4},
                    "Cl": {"Cl": 5}
                }

        Returns:
            GroFile: Constructed GroFile object
        """
        atoms_out = []
        atom_types_out = set()
        top_system_data = {}
        mol_atom_counts = {}

        # Reverse lookup: atom_type -> list of (mol_name, atom_name)
        atom_type_to_names = {}
        for mol_name, atom_dict in mol_atom_map.items():
            for atom_name, atom_type in atom_dict.items():
                atom_type_to_names.setdefault(atom_type, []).append((mol_name, atom_name))

        # Sort to ensure deterministic behavior
        for lst in atom_type_to_names.values():
            lst.sort()

        # Track mol_id → mol_name to avoid reselecting every time
        mol_id_to_type = {}

        for atom in self.atoms:
            atom_type = atom['atom_type']
            mol_id = atom['mol_index']

            name_matches = atom_type_to_names.get(atom_type)
            if not name_matches:
                raise ValueError(f"Atom type {atom_type} not found in mapping.")

            # Get or determine molecule type
            if mol_id not in mol_id_to_type:
                mol_id_to_type[mol_id] = name_matches[0][0]  # pick first match
            mol_type = mol_id_to_type[mol_id]

            # Count how many of each atom_name used in this molecule
            count = mol_atom_counts.setdefault(mol_id, defaultdict(int))

            # Get valid atom names for this molecule type
            valid_names = [name for mname, name in name_matches if mname == mol_type]
            used_names = set(count.keys())

            # Try unused names first, fallback to round-robin
            available = [name for name in valid_names if name not in used_names]
            if available:
                atom_name = available[0]
            else:
                idx = count[valid_names[0]] % len(valid_names)
                atom_name = valid_names[idx]

            count[atom_name] += 1

            atoms_out.append({
                'mol_name': mol_type,
                'mol_index': mol_id,
                'atom_name': atom_name,
                'atom_index': atom['atom_index'],
                'atom_type': atom_type,
                'charge': atom['charge'],
                'x': atom['x'],
                'y': atom['y'],
                'z': atom['z'],
            })
            atom_types_out.add(atom_type)

            # Track atom count per mol_type
            top_system_data.setdefault(mol_type, 0)
            top_system_data[mol_type] += 1

        temp_lmp = copy.deepcopy(self)
        temp_lmp.atoms = atoms_out
        temp_lmp.atom_types = sorted(atom_types_out)
        temp_lmp.top_system_data = top_system_data

        # Use the GroFile classmethod to build the gro object
        return GroFile.from_lmp(temp_lmp)