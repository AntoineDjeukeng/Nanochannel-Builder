import numpy as np
from md_simulation.builder.unit_cell import UnitCell


def load_unit_cells():
    cells = {}

    base_coords = np.array([
        [0.000000, 0.000000, 0.000000],
        [1.229756, 0.710000, 0.000000],
        [1.229756, 2.130000, 0.000000],
        [0.000000, 2.840000, 0.000000]
    ])

    # CAP
    cap = UnitCell(coords=base_coords, d=1.42, cz=3.604)
    cap.name = "CAP"
    cap.atom_names = ["CW"] * 4
    cap.atom_types = [6] * 4
    cap.satom_types = [1] * 4
    cap.charges = [1.0] * 4
    cells["CAP"] = cap

    # CAN
    can = UnitCell(coords=base_coords, d=1.42, cz=3.604)
    can.name = "CAN"
    can.atom_names = ["CW"] * 4
    can.atom_types = [7] * 4
    can.satom_types = [2] * 4
    can.charges = [-1.0] * 4
    cells["CAN"] = can

    # CAR
    car = UnitCell(coords=base_coords, d=1.42, cz=3.604)
    car.name = "CAR"
    car.atom_names = ["CW"] * 4
    car.atom_types = [5] * 4
    car.satom_types = [3] * 4
    car.charges = [0.0] * 4
    car.satom_types = [0] * 4
    cells["CAR"] = car

    # CAG
    cag = UnitCell(coords=base_coords, d=1.42, cz=3.604)
    cag.name = "CAG"
    cag.atom_names = ["CA"] * 4
    cag.atom_types = [8] * 4
    cag.satom_types = [4] * 4
    cag.charges = [0.0] * 4
    cells["CAG"] = cag

    # hBN
    hbn_coords = np.array([
        [0.000000, 0.000000, 0.000000],
        [1.259200937, 0.727000000, 0.000000],
        [0.000000000, 2.908000000, 0.000000],
        [1.259200937, 2.181000000, 0.000000]
    ])
    hbn = UnitCell(coords=hbn_coords, d=1.454, cz=3.35)
    hbn.name = "hBN"
    hbn.atom_names = ["B", "N", "N", "B"]
    hbn.atom_types = [1, 2, 2, 1]
    hbn.satom_types = [0] * 4
    hbn.charges = [0.3, -0.3, -0.3, 0.3]
    cells["hBN"] = hbn

    # GRA
    gra = UnitCell(coords=base_coords, d=1.418536585, cz=3.35)
    gra.name = "GRA"
    gra.atom_names = ["C"] * 4
    gra.atom_types = [3] * 4
    gra.satom_types = [0] * 4
    gra.charges = [0.0] * 4
    cells["GRA"] = gra
    
    return cells
