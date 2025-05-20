def lmpline(data, mode):
    """
    Parse or format a single LMP coordinate line.

    Parameters:
        data: dict or str
            If mode=='read', data is the line string to parse.
            If mode=='write', data is a dict with atom info.
        mode: str
            'read' to parse a line, 'write' to format a line.

    Returns:
        dict (if reading) or str (if writing)
    """
    if mode == 'read':
        return read(data)
    elif mode == 'write':
        return write(data)
    else:
        raise ValueError(f"Invalid mode '{mode}', expected 'read' or 'write'.")

def write(data):
    """
    Generate a LAMMPS data line from atom info dictionary.

    Expected keys in `data`:
        - atom_index (int)
        - mol_index (int)
        - atom_type (int)
        - atom_charge (float)
        - x, y, z (float)

    Returns:
        str: Formatted LAMMPS line, or None if input is invalid.
    """
    try:
        # Format each field for consistent spacing
        atom_index = f"{data['atom_index']:>6}  "
        mol_index = f"{data['mol_index']:>4}  "
        atom_type = f"{data['atom_type']:>4}  "
        atom_charge = f"{data['atom_charge']:>9.6f}  "
        x = f"{data['x']:>10.5f}"
        y = f"{data['y']:>10.5f}"
        z = f"{data['z']:>10.5f}"

        # Combine all fields into one line
        return f"{atom_index}{mol_index}{atom_type}{atom_charge}{x}{y}{z}\n"

    except KeyError as e:
        print(f"[ERROR] Missing key in input data: {e}")
        return None
    except (ValueError, TypeError) as e:
        print(f"[ERROR] Invalid value in input data: {e}")
        return None


def read(line):
    try:
        parts = line.strip().split()
        if len(parts) != 7:
            raise ValueError("LAMMPS line must have exactly 7 fields")

        atom_id = int(parts[0])
        mol_id = int(parts[1])
        atom_type = int(parts[2])
        charge = float(parts[3])
        x = float(parts[4]) / 10.0  # Convert from Å to nm
        y = float(parts[5]) / 10.0
        z = float(parts[6]) / 10.0

        return {
            'atom_index': atom_id,
            'mol_index': mol_id,
            'atom_type': atom_type,
            'charge': charge,
            'x': x,
            'y': y,
            'z': z,
            'mol_name': f'mol{mol_id}',      # default placeholder
            'atom_name': f'at{atom_type}'    # default placeholder
        }
    except (ValueError, IndexError) as e:
        print(f"Failed to parse LAMMPS line: {e}")
        return None
