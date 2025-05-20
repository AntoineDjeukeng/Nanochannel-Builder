from collections import defaultdict
from math import inf

def build_bonds_angles(atoms):
    """
    Build bonds and angles from a list of atoms.
    Assumes all atoms have 'mol_index', 'mol_name', 'atom_name', 'atom_index'.

    Returns:
        bonds: list of dicts with keys (index, func, i, j)
        angles: list of dicts with keys (index, func, i, j, k)
    """
    bonds = []
    angles = []
    top_data = {}
    atom_types = {}

    vbox = [inf, inf, inf, -inf, -inf, -inf]  # [xmax, ymax, zmax, xmin, ymin, zmin]

    # Group atoms by molecule
    molecules = defaultdict(list)
    for atom in atoms:
        if atom['atom_name'] not in atom_types:
            atom_types[atom['atom_name']] = len(atom_types) + 1
        vbox[0] = min(vbox[0], atom['x'])
        vbox[3] = max(vbox[3], atom['x'])
        vbox[1] = min(vbox[1], atom['y'])
        vbox[4] = max(vbox[4], atom['y'])
        vbox[2] = min(vbox[2], atom['z'])
        vbox[5] = max(vbox[5], atom['z'])
        molecules[atom['mol_index']].append(atom)
    
    bond_index = 1
    angle_index = 1

    for mol_atoms in molecules.values():
        
        mol_name = mol_atoms[0]['mol_name']
        if mol_name not in top_data:
            top_data[mol_name] = 0
        
        top_data[mol_name] += 1
        if mol_name != 'SOL':
            continue

        atom_map = {atom['atom_name']: atom['atom_index'] for atom in mol_atoms}

        # Check for water atoms
        try:
            n1 = atom_map['OW']
            n2 = atom_map['HW1']
            n3 = atom_map['HW2']

            bonds.append({"index": bond_index, "func": 1, "i": n1, "j": n2})
            bond_index += 1
            bonds.append({"index": bond_index, "func": 1, "i": n1, "j": n3})
            bond_index += 1

            angles.append({"index": angle_index, "func": 1, "i": n2, "j": n1, "k": n3})
            angle_index += 1
        except KeyError:
            # skip malformed water molecules
            continue
    
    return bonds, angles, top_data, vbox, atom_types
