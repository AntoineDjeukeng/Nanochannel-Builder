import os
from collections import defaultdict
from md_simulation.io.gro_file import GroFile
from md_simulation.io.lmp_file import LmpFile

def test_write_lmp_file_from_gro(tmp_path):
    gro_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'react_pos_7_1.gro')
    gro = GroFile.read(gro_path)
    lmp = LmpFile.from_gro(gro)

    output_path = tmp_path / "test_output.data"
    lmp.write(output_path)

    assert output_path.exists()

    # Check that atoms section exists and atoms are grouped by molecule contiguously
    with open(output_path) as f:
        lines = f.readlines()

    # Extract atoms lines after "Atoms" section
    atom_lines = []
    in_atoms_section = False
    for line in lines:
        if line.strip() == "Atoms":
            in_atoms_section = True
            continue
        if in_atoms_section:
            if not line.strip():
                break
            atom_lines.append(line.strip())

    mol_atoms = defaultdict(list)
    for line in atom_lines:
        parts = line.split()
        # Format: atom_id mol_id atom_type charge x y z
        atom_id = int(parts[0])
        mol_id = int(parts[1])
        mol_atoms[mol_id].append(atom_id)

    # Assert atoms within each molecule appear contiguously and in ascending order
    for mol_id, atoms in mol_atoms.items():
        sorted_atoms = sorted(atoms)
        assert atoms == sorted_atoms, f"Atoms of molecule {mol_id} not contiguous or ordered"

    print(f"✔ Atoms for {len(mol_atoms)} molecules correctly grouped and ordered")
