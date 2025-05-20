import numpy as np
from md_simulation.builder.surface import Surface
from md_simulation.builder.unit_cells import load_unit_cells
from md_simulation.io.gro_file import GroFile
from md_simulation.io.lmp_file import LmpFile
import argparse
def select_unit_cells(system, name_surf, unit_cells):
    if system == "CDI":
        unit_cell_map = {
            'positive': unit_cells["CAP"],
            'negative': unit_cells["CAN"],
            'neutral': unit_cells["CAR"],
        }
        unit_cell1 = unit_cell_map[name_surf]
        unit_cell2 = unit_cells["CAG"]
    elif system == "hBNs":
        unit_cell1 = unit_cells["hBN"]
        unit_cell2 = unit_cells["GRA"]
    else:
        raise ValueError(f"Unknown system: {system}")
    return unit_cell1, unit_cell2


def build_layered_surface(unit_cell1, unit_cell2, Lx, Ly, charges, structure, z_gap=7, LX=700, write_type=0, shift_y=0):
    surf1 = Surface(Lx, Ly)
    surf1.create(unit_cell1, charges[0])

    surf2 = Surface(Lx, Ly)
    surf2.create(unit_cell2, charges[1])

    cz1, cz2 = unit_cell1.box[2], unit_cell2.box[2]
    cy1, cy2 = unit_cell1.box[1] / 3, unit_cell2.box[1] / 3
    zt = 0.5 * (LX - surf1.box[0])

    h1, h2, y1, y2 = [], [], [], []
    z_current = -cz1

    for i, layer in enumerate(structure):
        gap = z_gap - 0.5 * (cz1 + cz2) if i == len(structure) // 2 else 0
        shift = 0.5 * (-1) ** i

        if layer == 1:
            z_current += cz1 + gap
            h1.append(z_current)
            y1.append(shift * cy1 * shift_y)
        elif layer == 2:
            z_current += cz2 + gap
            h2.append(z_current)
            y2.append(shift * cy2 * shift_y)

    Lz = z_current + cz2
    nmole = natoms = 0
    lmp_lines = []
    gro_lines = []

    for z, y in zip(h1, y1):
        v = [zt, y, z]
        ll, gl, n_mol, n_at = surf1.translate(v, nmole, natoms, write_type)
        nmole += n_mol
        natoms += n_at
        lmp_lines += ll
        gro_lines += gl

    if surf1.total_charge and natoms and surf1.charge_nultiplicity != 0:
        surf1.charge_nultiplicity = int(surf1.total_charge + 0.5) / natoms

    for z, y in zip(h2, y2):
        v = [zt, y, z]
        ll, gl, n_mol, n_at = surf2.translate(v, nmole, natoms, write_type)
        nmole += n_mol
        natoms += n_at
        lmp_lines += ll
        gro_lines += gl

    if surf2.total_charge and natoms and surf2.charge_nultiplicity != 0:
        surf2.charge_nultiplicity = int(surf2.total_charge + 0.5) / natoms

    box = surf1.box.copy()
    box[2] = Lz
    box[0] = LX

    return lmp_lines, gro_lines, natoms, box, surf1, surf2


def write_output_files(gro_lines, lammps_lines, natoms, box, surf1, surf2, write_type, title, filename):
    gro = GroFile()
    gro.title = title
    gro.natoms = natoms
    gro.box = [0.1 * x for x in box]
    gro.atoms = gro_lines

    lmp = LmpFile()
    lmp.title = title
    lmp.natoms = natoms
    lmp.box = np.array([[0, 0.1 * box[0]], [0, 0.1 * box[1]], [0, 0.1 * box[2]]])
    lmp.atoms = lammps_lines

    names = {}
    types_list = []
    atom_counter = 0

    for surf in [surf1, surf2]:
        for i, name in enumerate(surf.molecule.atom_names):
            charge = surf.charge_nultiplicity * surf.molecule.charges[i]
            type_id = surf.molecule.satom_types[i]
            entry = [type_id, surf.molecule.name, name, charge]

            if name not in names:
                names[name] = len(names) + 1
                types_list.append(entry)
            elif write_type != 0:
                types_list.append(entry)
            atom_counter += 1

    lmp.atom_types = names
    gro.temp_top = {
        "atom_types": types_list,
        "moecule_types": {
            surf1.molecule.name: Surface.all_surfaces[surf1.molecule.name] * (1 if write_type else len(surf1.molecule.atom_names)),
            surf2.molecule.name: Surface.all_surfaces[surf2.molecule.name] * (1 if write_type else len(surf2.molecule.atom_names))
        }
    }

    gro.write(f"{filename}.gro")
    lmp.write(f"{filename}.data")



 

import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Build a layered membrane surface.")
    
    parser.add_argument("--system", type=str, required=True, choices=["CDI", "hBNs"],
                        help="System type (CDI or hBNs)")
    parser.add_argument("--name_surf", type=str, required=True,
                        help="Name of the surface (e.g., positive)")
    parser.add_argument("--output", type=str, required=True,
                        help="Base name for output files")
    parser.add_argument("--Lx", type=int, required=True,
                        help="Number of unit cells along x")
    parser.add_argument("--Ly", type=int, required=True,
                        help="Number of unit cells along y")
    parser.add_argument("--LX", type=int, required=True,
                        help="Length of the box along X (angstrom)")
    parser.add_argument("--z_gap", type=float, required=True,
                        help="Gap between layers in z (angstrom)")
    parser.add_argument("--shift_y", type=float, required=True,
                        help="Y-axis shift between top and bottom surfaces")
    parser.add_argument("--write_type", type=int, required=True,
                        help="Write type: 1 = GRO+LAMMPS, 0 = GRO only")
    parser.add_argument("--structure", type=int, nargs='+', required=True,
                        help="Layer structure as a list of integers (e.g., --structure 2 2 1 1 2 2)")
    
    return parser.parse_args()

def main():
    args = parse_args()
    unit_cells = load_unit_cells()

    structure = args.structure
    charges = [0.0104 * 2, 0] if args.system == "CDI" else [0, 0]
    title = "Generated by MD Simulation"

    unit_cell1, unit_cell2 = select_unit_cells(args.system, args.name_surf, unit_cells)
    lammps_lines, gro_lines, natoms, box, surf1, surf2 = build_layered_surface(
        unit_cell1, unit_cell2,
        args.Lx, args.Ly, charges,
        structure, args.z_gap, args.LX, args.write_type, args.shift_y
    )

    write_output_files(
        gro_lines, lammps_lines, natoms,
        box, surf1, surf2,
        args.write_type, title, f"data/{args.output}"
    )

if __name__ == "__main__":
    main()

# python scripts/build_membrane.py \
#   --system CDI \
#   --name_surf positive \
#   --output mymembrane \
#   --Lx 20 \
#   --Ly 52 \
#   --LX 200 \
#   --z_gap 7.0 \
#   --shift_y 0.0 \
#   --write_type 1 \
#   --structure 2 2 1 1 2 2