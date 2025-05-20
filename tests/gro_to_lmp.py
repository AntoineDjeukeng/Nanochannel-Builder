import os
from md_simulation.io.gro_file import GroFile
from md_simulation.io.lmp_file import LmpFile

def test_read_gro_and_write_lmp():
    # Path to test GRO file in the data folder
    gro_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'react_pos_7_1.gro')
    assert os.path.exists(gro_path), f"GRO file not found: {gro_path}"

    # Read GRO file
    gro = GroFile.read(gro_path)

    # Basic sanity checks on GRO file read
    assert gro.natoms > 0
    assert len(gro.atoms) == gro.natoms

    # Convert GRO to LmpFile
    lmp = LmpFile.from_gro(gro)

    # Prepare output path for LMP file
    lmp_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'test_output.data')
    os.makedirs(os.path.dirname(lmp_path), exist_ok=True)

    # Write LMP file and catch any errors
    try:
        lmp.write(lmp_path)
    except Exception as e:
        assert False, f"Failed to write LMP file: {e}"

    # Verify LMP file was created and is non-empty
    assert os.path.exists(lmp_path), "LMP output file was not created"
    assert os.path.getsize(lmp_path) > 0, "LMP output file is empty"

    print(f"LMP file successfully written to {lmp_path}")
