***MD simulation***
# Molecular Dynamics Simulation of solid-liquid interface

This project prepares input files for molecular dynamics (MD) simulations using **GROMACS** and **LAMMPS**. Currently, it includes two main features:

1. ## Conversion between **LAMMPS** and **GROMACS** formats

a. **LAMMPS** to **GROMACS**

An example is provided in the test folder (test_lmp_to_gro.py) that converts a **LAMMPS** data file (test.lmp) to a **GROMACS** .gro file (data/test.gro).

The code uses a dictionary to map **LAMMPS** atom types to **GROMACS** atom and residue names. For example:

```python
import os
from md_simulation.io.lmp_file import LmpFile

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
test_gro_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'test_write.gro')
gro.write(test_gro_path)
```
If **HW1** and **HW2** have the same type, the code will automatically resolve naming issues.

b. **GROMACS** to **LAMMPS**

The process is similar, but the atom types are automatically assigned based on order of appearance. No user mapping is needed:

```python
import os
from md_simulation.io.gro_file import GroFile
from md_simulation.io.lmp_file import LmpFile

gro_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'react_pos_7_1.gro')
assert os.path.exists(gro_path), f"GRO file not found: {gro_path}"

gro = GroFile.read(gro_path)
lmp_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'test_output.data')
os.makedirs(os.path.dirname(lmp_path), exist_ok=True)

lmp.write(lmp_path)
```
2. # Membrane Builder

Currently, the script can build two types of layered membranes:

- Graphitic membrane (**CDI**)

- Hexagonal boron nitride (**hBN**)

The input is parsed using:
```python
parser = argparse.ArgumentParser(description="Build a layered membrane surface.")
...
python scripts/build_membrane.py\   
--system CDI   \
--name_surf positive \  
--output CDI2\   
--Lx 20   \
--Ly 52   \
--LX 200   \
--z_gap 7.0  \ 
--shift_y 0.0   \
--write_type 1  \ 
--structure 2 2 1 1 2 2

```
Key arguments:

- `system`: Membrane type. Currently supports "CDI" and "hBNs".

- `name_surf`: Layer label (useful for charged vs. neutral layers in CDI).

- `output`: Base name for output files.

- `Lx`, `Ly`: approximated dimensions of the simulation box along X and Y (angstrom). These values are used to calculate the number of unit cells in the X and Y directions.

- `LX`: Simulation box length along X (angstrom), adjusting this value will create a buffer zone for a reservoir in the x direction. See the example below.

- `z_gap`: Gap between top and bottom layers (angstrom).

- `shift_y`: Y-axis shift between layers (used to create staggered structures).

- `write_type`: 1 = write each atom with a different residue ID (GRO+LAMMPS), 0 = all atoms in the same uint cell have the same residue ID.

- `structure`: Defines the layer pattern as a list of integers. Each integer represents a unit cell type:

- For CDI: 1 = contact (charged) layer, 2 = neutral graphene layer

- For hBNs: 1 = hBN layer, 2 = neutral graphene layer

***Example***: --structure 2 2 1 1 2 2 means:

- In the case fo CDI: 4 layers of neutral graphene and 2 contact layers to the electrolyte.

**Example CDI unit cell**:
```python
car = UnitCell(coords=base_coords, d=1.42, cz=3.604)
car.name = "CAR"
car.atom_names = ["CW"] * 4
car.atom_types = [5] * 4
car.charges = [0.0] * 4
cells["CAR"] = car
```
- `d`: Distance between atoms in a layer (C–C bond length).

- `cz`: Distance between top and bottom layers.

If write_type is set to 1, atoms in the same layer get different residue IDs.
## Example Membranes

### CDI Membrane:
![CDI](/data/CDI1.png)
###chBN Membrane:
![hBNs](/data/hBNs1.png)
