def groline(data, mode):
    """
    Parse or format a single GRO coordinate line.

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
    try:
        # Limit indices to max 99999 to fit GRO format
        mol_index = data['mol_index'] % 100000
        atom_index = data['atom_index'] % 100000

        mol_index_str = f"{mol_index:5d}"
        mol_name_str = f"{data['mol_name'].strip():<5}"
        atom_name_str = f"{data['atom_name'].strip():>5}"
        atom_index_str = f"{atom_index:5d}"

        # GRO files use nm units, convert from Å if necessary
        # Confirm your coordinate unit convention
        x_str = f"{ data['x']:8.3f}"
        y_str = f"{ data['y']:8.3f}"
        z_str = f"{ data['z']:8.3f}"

        line = f"{mol_index_str}{mol_name_str}{atom_name_str}{atom_index_str}{x_str}{y_str}{z_str}\n"
        return line

    except KeyError as e:
        print(f"Missing key in input data: {e}")
        return None
    except (ValueError, TypeError) as e:
        print(f"Error processing data: {e}")
        return None

def read(line):
    try:
        mol_index = int(line[0:5].strip())
        mol_name = line[5:10].strip()
        atom_name = line[10:15].strip()
        atom_index = int(line[15:20].strip())
        x = float(line[20:28].strip())
        y = float(line[28:36].strip())
        z = float(line[36:44].strip())

        return {
            'mol_index': mol_index,
            'mol_name': mol_name,
            'atom_name': atom_name,
            'atom_index': atom_index,
            'x': x,
            'y': y,
            'z': z
        }

    except (IndexError, ValueError) as e:
        print(f"Error processing line: {line.strip()}")
        print(f"Error details: {e}")
        return None
