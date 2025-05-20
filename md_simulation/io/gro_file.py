from md_simulation.io.gro_line import groline
from md_simulation.topology.build_topology import build_bonds_angles  # You must create this
import numpy as np

class GroFile:
    """
    Represents a GROMACS .gro coordinate file.

    Attributes:
        title (str): Title/comment line from the .gro file.
        natoms (int): Number of atoms.
        atoms (list of dict): Atom records.
        box (str): Box vector line.
        bonds (list of tuples): Pairs of bonded atom indices.
        angles (list of tuples): Triplets of angle atom indices.
    """

    def __init__(self):
        self.title = ""
        self.natoms = 0
        self.atoms = []
        self.atom_types = {}
        self.box = []
        self.temp_top = None
        self.lmp_box = np.zeros((3, 2))
        self.bonds = []
        self.angles = []
        self.top_system_data = {}

    @classmethod
    def from_lmp(cls, lmp_obj):
        """
        Build a GroFile from a LmpFile instance.

        Parameters:
            lmp_obj (LmpFile): Parsed LmpFile instance.

        Returns:
            GroFile: A new instance containing GRO-compatible data.
        """
        obj = cls()

        # Copy relevant fields from LmpFile object
        obj.title = getattr(lmp_obj, 'title', "Generated from LAMMPS data")
        obj.atoms = lmp_obj.atoms.copy() if hasattr(lmp_obj, 'atoms') else []
        obj.atom_types = lmp_obj.atom_types.copy() if hasattr(lmp_obj, 'atom_types') else []
        obj.natoms = len(obj.atoms)
        obj.box = lmp_obj.gro_box if hasattr(lmp_obj, 'gro_box') else None
        obj.bonds = lmp_obj.bonds.copy() if hasattr(lmp_obj, 'bonds') and lmp_obj.bonds else []
        obj.angles = lmp_obj.angles.copy() if hasattr(lmp_obj, 'angles') and lmp_obj.angles else []
        obj.top_system_data = lmp_obj.top_system_data.copy() if hasattr(lmp_obj, 'top_system_data') else {}

        return obj


    @classmethod
    def read(cls, filename, with_topology=True):
        """
        Read a .gro file and return a GroFile instance.

        Parameters:
            filename (str): Path to the .gro file.
            with_topology (bool): Whether to build bonds/angles.

        Returns:
            GroFile: instance populated with data from the file.
        """
        obj = cls()
        with open(filename, 'r') as f:
            obj.title = f.readline().rstrip('\n')
            obj.natoms = int(f.readline().strip())
            for _ in range(obj.natoms):
                line = f.readline()
                atom = groline(line, 'read')
                if atom is None:
                    raise ValueError(f"Failed to parse atom line: {line}")
                obj.atoms.append(atom)
            box_line = f.readline().split()
            obj.box = [float(coord) for coord in box_line]

        if with_topology:
            obj.bonds, obj.angles ,obj.top_system_data, vbox,obj.atom_types= build_bonds_angles(obj.atoms)
            
            a = vbox[3] - vbox[0]
            b = vbox[4] - vbox[1]
            c = vbox[5] - vbox[2]
            obj.lmp_box[0][0] = 0.5 * (a - obj.box[0])
            obj.lmp_box[0][1] = 0.5 * (a + obj.box[0])
            obj.lmp_box[1][0] = 0.5 * (b - obj.box[1])
            obj.lmp_box[1][1] = 0.5 * (b + obj.box[1])
            obj.lmp_box[2][0] = 0.5 * (c - obj.box[2])
            obj.lmp_box[2][1] = 0.5 * (c + obj.box[2])

        return obj

    def write(self, filename: str) -> None:
        """
        Write the GroFile data to a .gro file.

        Parameters:
            filename (str): Path where to save the file.
        """
        if self.temp_top:
            temp_file = filename.replace(".gro", "_temp.top")
            with open(temp_file, 'w') as f:
                f.write("list of atom types and their charge\n")
                datat = self.temp_top["atom_types"]
                for data in datat:
                    f.write(f"{data[1]} {data[0]} {data[2]} {data[3]}\n")
                f.write("[system]\n")
                f.write(f"{self.title}\n")
                f.write("[molecules]\n")
                for mol_name, count in self.temp_top["moecule_types"].items():
                    f.write(f"{mol_name} {count}\n")
        with open(filename, 'w') as f:
            # Write the title
            if not self.title:
                raise ValueError("GroFile title is empty or missing")
            f.write(f"{self.title}\n")

            # Write number of atoms
            natoms = len(self.atoms)
            if natoms == 0:
                raise ValueError("GroFile has no atoms to write")
            f.write(f"{natoms}\n")

            # Write atom lines
            for atom in self.atoms:
                line = groline(atom, mode='write')  # assumes groline returns formatted string with newline
                if line is None:
                    raise ValueError(f"Failed to format atom data: {atom}")
                f.write(line)

            # Write box dimensions
            # Handle if box is list/tuple of floats or string
            if isinstance(self.box, (list, tuple)):
                box_str = " ".join(f"{x:.6f}" for x in self.box)
            elif isinstance(self.box, str):
                box_str = self.box.strip()
            else:
                raise TypeError(f"Unsupported type for box attribute: {type(self.box)}")
            f.write(f"{box_str}\n")

        def add_atom(self, atom_dict):
            """
            Add an atom dictionary to the GroFile.

            Parameters:
                atom_dict (dict): Atom data with required keys.
            """
            self.atoms.append(atom_dict)
            self.natoms = len(self.atoms)
