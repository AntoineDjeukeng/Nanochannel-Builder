import numpy as np
from md_simulation.io.gro_line import groline
from md_simulation.io.lmp_line import lmpline

class Surface:
    all_surfaces = {}
    total_charge = 0

    def __init__(self, a, b):
        self.a = a  # Desired x length
        self.b = b  # Desired y length
        self.box = []  # Final box dimensions
        self.sturcture = None  # Will store the structure of the surface
        self.coords = None  # Will store atom coordinates
        self.molecule = None


    def create(self, molecule,m):
        self.molecule = molecule
        self.charge_nultiplicity = m
        name = molecule.name

        # Track how many times a surface of a given type has been used
        if name not in Surface.all_surfaces:
            Surface.all_surfaces[name] = 0

        mbox = molecule.box  # Unit cell box size

        # Determine replication counts
        m = round(self.a / mbox[0])
        n = round(self.b / mbox[1])

        # Update surface box size
        self.box = [m * mbox[0], n * mbox[1], mbox[2]]

        # Build surface by replicating the unit cell
        ai = 0  # Atom index counter
        mi = 0  # Molecule index counter
        self.coords = []

        for i in range(m):
            for j in range(n):
                mi += 1
                offset = np.array([(i + 0.5) * mbox[0], (j + 0.5) * mbox[1], 0.5 * mbox[2]])
                atoms = molecule.build(offset)

                newmol = np.zeros((atoms.shape[0], 5))
                for k in range(atoms.shape[0]):
                    ai += 1
                    newmol[k][0] = ai       # Atom index
                    newmol[k][1] = mi       # Molecule index
                    newmol[k][2:] = atoms[k]

                self.coords.append(newmol)

        # Stack all molecule arrays into one array
        self.coords = np.vstack(self.coords)

    def translate(self, vectors, a=0, b=0, t=0):
        """
        Translate and format the surface for LAMMPS and GRO files.
        Args:
            vectors: 3-element array for translation [x, y, z]
            a: molecule index offset
            b: atom index offset
        Returns:
            lammps_lines, gro_lines, nmol, natoms
        """
        if self.coords is None:
            self.create(self.molecule)

        coords = self.coords
        lammps_lines = []
        gro_lines = []

        last_mol_id = -1
        nmole = 0
        natoms = 0
        if self.charge_nultiplicity != 0:
            self.charge_nultiplicity = int(coords.shape[0]*self.charge_nultiplicity+0.5)/coords.shape[0]
        Surface.total_charge += coords.shape[0]*self.charge_nultiplicity
        for i in range(coords.shape[0]):
            atom_index = int(coords[i, 0]) + b
            mol_index = int(coords[i, 1]) + a
            x, y, z = coords[i, 2:] + vectors

            # Detect new molecule
            if mol_index != last_mol_id :
                Surface.all_surfaces[self.molecule.name] += 1
                nmole += 1
                last_mol_id = mol_index
                atom_local_index = 0
            else:
                atom_local_index += 1
            if t == 0:
                mol_index = atom_index
            m_charge = self.molecule.charges[atom_local_index]
            m_charge = m_charge * self.charge_nultiplicity
            data = {
                'atom_index': atom_index,
                'mol_index': mol_index,
                'mol_name': self.molecule.name,
                'atom_name': self.molecule.atom_names[atom_local_index],
                'atom_type': self.molecule.satom_types[atom_local_index],
                'atom_charge': m_charge,
                'x': 0.1 *x,
                'y': 0.1 *y,
                'z': 0.1 *z
            }
            lammps_lines.append(data)
            gro_lines.append(data)
            natoms += 1
        return lammps_lines, gro_lines, nmole, natoms
