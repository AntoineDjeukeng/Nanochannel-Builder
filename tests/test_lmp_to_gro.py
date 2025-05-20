import os
from md_simulation.io.lmp_file import LmpFile

def test_lmp_to_gro_conversion():
    # Path to test LAMMPS data file
    lmp_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'test_output.data')
    assert os.path.exists(lmp_path), f"LAMMPS data file not found: {lmp_path}"

    # Read LAMMPS data
    lmp = LmpFile.read(lmp_path)

    # Mapping of molecule types and atom names to atom types
    mol_atom_map = {
        "SOL": {"OW": 1, "HW1": 2, "HW2": 3},
        "Na": {"NA": 4},
        "Cl": {"Cl": 5}
    }

    # Convert to GRO
    gro = lmp.to_gro(mol_atom_map)

    # Basic assertions
    assert gro.natoms == len(gro.atoms)
    assert gro.title != ""
    assert gro.box is not None
    assert isinstance(gro.atom_types, list)
    assert gro.top_system_data

    # Test writing to file
    test_gro_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'test_write.gro')
    gro.write(test_gro_path)
    assert os.path.exists(test_gro_path), "GRO output file was not created"

    # Optional: Basic content check (file not empty, line counts)
    with open(test_gro_path, 'r') as f:
        lines = f.readlines()
    assert len(lines) > gro.natoms + 2, "GRO file seems incomplete"

    print("\n--- LAMMPS to GRO Conversion and Write Summary ---")
    print(f"Title: {gro.title}")
    print(f"Number of atoms: {gro.natoms}")
    print(f"Box: {gro.box}")
    print(f"Output GRO file: {test_gro_path}")

    # Optional cleanup if you don't want to keep the file
    # os.remove(test_gro_path)
