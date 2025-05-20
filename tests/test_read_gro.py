import os
from collections import Counter
from md_simulation.io.gro_file import GroFile

def test_read_gro_file():
    # Path to test GRO file in the data folder
    gro_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'react_pos_7_1.gro')
    assert os.path.exists(gro_path), f"GRO file not found: {gro_path}"

    # Read the GRO file
    gro = GroFile.read(gro_path)

    # Basic assertions
    assert gro.natoms > 0
    assert len(gro.atoms) == gro.natoms
    assert gro.title != ""
    assert gro.box != ""

    # Count molecules by mol_name
    mol_counter = Counter(atom['mol_name'] for atom in gro.atoms)

    print("\n--- GroFile Read Summary ---")
    print(f"Title: {gro.title}")
    print(f"Number of atoms: {gro.natoms}")
    print(f"LMP box: {gro.lmp_box}")
    print(f"Box: {gro.box}")

    print("Molecule counts:")
    for mol_name, number in gro.top_system_data.items():
        print(f"  {mol_name}: {number} atoms")
    print(f"atom_types: {gro.atom_types}")
    print("\nFirst 5 atoms:")
    for atom in gro.atoms[:5]:
        print(atom)

    if hasattr(gro, 'bonds') and gro.bonds:
        print("\nFirst 5 bonds:")
        for bond in gro.bonds[:5]:
            print(f"  {bond.strip()}" if isinstance(bond, str) else bond)

    if hasattr(gro, 'angles') and gro.angles:
        print("\nFirst 5 angles:")
        for angle in gro.angles[:5]:
            print(f"  {angle.strip()}" if isinstance(angle, str) else angle)

    print("\n----------------------------\n")
