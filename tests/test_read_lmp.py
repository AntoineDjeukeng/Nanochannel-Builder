import os
from md_simulation.io.lmp_file import LmpFile

def test_read_lmp_file():
    # Path to test LMP file in the data folder
    lmp_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'test_output.data')
    assert os.path.exists(lmp_path), f"LMP file not found: {lmp_path}"

    # Read LMP file
    lmp = LmpFile.read(lmp_path)

    # Basic sanity checks
    assert lmp.natoms > 0
    assert len(lmp.atoms) == lmp.natoms
    assert lmp.title != ""

    print("\n--- LmpFile Read Summary ---")
    print(f"Title: {lmp.title}")
    print(f"Number of atoms: {lmp.natoms}")
    print(f"Number of bonds: {len(lmp.bonds) if hasattr(lmp, 'bonds') else 0}")
    print(f"Number of angles: {len(lmp.angles) if hasattr(lmp, 'angles') else 0}")

    # Print first 5 atoms for a quick look
    print("\nFirst 5 atoms:")
    for atom in lmp.atoms[:5]:
        print(atom)

    print("\n----------------------------\n")
