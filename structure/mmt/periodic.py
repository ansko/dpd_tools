import random




def mmt_sheet_double_diagonals_periodic(**kwargs):
    '''
    Create a sigle periodic sheet of MMT having thickness of 2 atoms.
    Next bonds (if exist) are set to 2*bead_radius:
        top_{i-1}{j}, top_{i}{j}
        top_{i}{j-1}, top_{i}{j}
        bottom_{i-1}{j}, bottom_{i}{j}
        bottom_{i}{j-1}, bottom_{i}{j}
        top_{i-1}{j-1}, bottom{i}{j}   # diagonal
        bottom_{i-1}{j-1}, top{i}{j}   # diagonal
    Required parameters:
        cell_side - size of the cell (i.e. size of sheet)
        bead_radius - bead radius in DPD
            (bond_length equals to 2*bead_radius)
    Optional parameters:
        atom_type=1 - type to assign to all atoms
        bond_type_edge=1 - type to assign to bonds between neighbors
        bond_type_diagonal=2 - type to assign to diagonal bonds
        x=0,y=0,z=0 - coordinates of platelet's geometrical center
        start_atom_id=1 - id of the first atom
        start_bond_id=1 - id of the first bond
        start_angle_id=1 - id of the first angle
        bead_charge - charge of every charged bead
        charged_count - count of beads to be charged
    '''
    bead_radius = kwargs['bead_radius']
    atom_type = kwargs['atom_type'] if 'atom_type' in kwargs.keys() else 1
    try:
        bond_type_edge = kwargs['bond_type_edge']
        bond_type_diagonal = kwargs['bond_type_diagonal']
    except KeyError:
        bond_type_edge = 1
        bond_type_diagonal = 2
    x = kwargs['x'] if 'x' in kwargs.keys() else 0
    y = kwargs['y'] if 'y' in kwargs.keys() else 0
    z = kwargs['z'] if 'z' in kwargs.keys() else 0
    cell_side = kwargs['cell_side']
    atom_id = kwargs['start_atom_id'] if 'start_atom_id' in kwargs.keys() else 1
    bond_id = kwargs['start_bond_id'] if 'start_bond_id' in kwargs.keys() else 1
    angle_id = kwargs['start_angle_id'] if 'start_angle_id' in kwargs.keys() else 1
    structure = {'atoms': {}, 'bonds': {}}
    assert(not (len(structure['atoms'].keys()) != 0 and atom_id == 1))
    assert(not (len(structure['bonds'].keys()) != 0 and bond_id == 1))
    atoms_before = len(structure['atoms'])
    atom_ids = {}
    # Add atoms
    for idxx in range(cell_side):
        atom_ids[idxx] = {}
        dx = idxx * 2*bead_radius
        for idxy in range(cell_side):
            atom_ids[idxx][idxy] = {}
            dy = idxy * 2*bead_radius
            structure['atoms'][atom_id] = {
                'x': x + dx,
                'y': y + dy,
                'z': z - bead_radius,
                'charge': 0,
                'type': atom_type,
                'phase': 'filler'
            }
            atom_ids[idxx][idxy]['bottom'] = atom_id
            atom_id += 1
            structure['atoms'][atom_id] = {
                'x': x + dx,
                'y': y + dy,
                'z': z + bead_radius,
                'charge': 0,
                'type': atom_type,
                'phase': 'filler'
            }
            atom_ids[idxx][idxy]['top'] = atom_id
            atom_id += 1
            structure['bonds'][bond_id] = {
                'atoms': [atom_ids[idxx][idxy]['top'],
                          atom_ids[idxx][idxy]['bottom']],
                'type': bond_type_edge
            }
            bond_id += 1
    # Add other bonds
    for idxx in atom_ids.keys():
        for idxy in atom_ids[idxx].keys():
            if idxx == 0:
                idxx_other_minus = cell_side - 1
                idxx_other_plus = 1
            elif idxx == cell_side - 1:
                idxx_other_minus = idxx - 1
                idxx_other_plus = 0
            else:
                idxx_other_minus = idxx - 1
                idxx_other_plus = idxx + 1
            if idxy == 0:
                idxy_other_minus = cell_side - 1
                idxy_other_plus = 1
            elif idxy == cell_side - 1:
                idxy_other_minus = idxy - 1
                idxy_other_plus = 0
            else:
                idxy_other_minus = idxy - 1
                idxy_other_plus = idxy + 1
            # edge bonds
            structure['bonds'][bond_id] = {
                'atoms': [atom_ids[idxx_other_minus][idxy]['top'],
                          atom_ids[idxx][idxy]['top']],
                'type': bond_type_edge
            }
            bond_id += 1
            structure['bonds'][bond_id] = {
                'atoms': [atom_ids[idxx_other_minus][idxy]['bottom'],
                          atom_ids[idxx][idxy]['bottom']],
                'type': bond_type_edge
            }
            bond_id += 1
            structure['bonds'][bond_id] = {
                'atoms': [atom_ids[idxx][idxy_other_minus]['top'],
                          atom_ids[idxx][idxy]['top']],
                'type': bond_type_edge
            }
            bond_id += 1
            structure['bonds'][bond_id] = {
                'atoms': [atom_ids[idxx][idxy_other_minus]['bottom'],
                          atom_ids[idxx][idxy]['bottom']],
                'type': bond_type_edge
            }
            bond_id += 1
            # diagonal
            structure['bonds'][bond_id] = {
                'atoms': [atom_ids[idxx_other_minus][idxy_other_minus]['bottom'],
                          atom_ids[idxx][idxy]['top']],
                'type': bond_type_diagonal
            }
            bond_id += 1
            structure['bonds'][bond_id] = {
                'atoms': [atom_ids[idxx_other_minus][idxy_other_plus]['bottom'],
                          atom_ids[idxx][idxy]['top']],
                'type': bond_type_diagonal
            }
            bond_id += 1
            structure['bonds'][bond_id] = {
                'atoms': [atom_ids[idxx_other_plus][idxy_other_minus]['bottom'],
                          atom_ids[idxx][idxy]['top']],
                'type': bond_type_diagonal
            }
            bond_id += 1

            structure['bonds'][bond_id] = {
                'atoms': [atom_ids[idxx_other_plus][idxy_other_plus]['bottom'],
                          atom_ids[idxx][idxy]['top']],
                'type': bond_type_diagonal
            }
            bond_id += 1
    # Assign (or not) charges
    if 'bead_charge' not in kwargs.keys():
        return structure
    bead_charge = kwargs['bead_charge']
    charged_count = kwargs['charged_count']
    chosen_atoms = set()
    while len(chosen_atoms) < charged_count:
        chosen_atoms.add(random.choice(
            range(atoms_before + 1, len(structure['atoms']))))
    for k in chosen_atoms:
        structure['atoms'][k]['charge'] = bead_charge
    return structure
