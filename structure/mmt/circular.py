# The function is made to add a circular MMT platelet into system.
# The platelet is parallel to the xy plane.
#
# MMT platelet consists of 2 layers of beads (top and bottom). Beads
# are all of one type, and are connected by bonds of two types: former
# bond type is for neighboring beads while former is for diagonal
# bonds:
#     top_{i-1}{j}, top_{i}{j}       # neighboring
#     top_{i}{j-1}, top_{i}{j}       # neighboring
#     bottom_{i-1}{j}, bottom_{i}{j} # neighboring
#     bottom_{i}{j-1}, bottom_{i}{j} # neighboring
#     top_{i-1}{j-1}, bottom{i}{j}   # diagonal
#     bottom_{i-1}{j-1}, top{i}{j}   # diagonal
#
# Required parameters are:
#     platelet_radius - radius of the circular sheet (in beads). When
#         it is odd, atoms with single bond may appear.
#     bead_radius - bead radius in DPD (they are really unitless, but i
#         treat them like they are expressed in angstroms). Bond 
#         lengths equal to 2*bead_radius for neighboring and
#         2*sqrt(3)*bead_radius for diagonal.
# Optional parameters are:
#     atom_type=1 - type to assign to all atoms
#     bond_type_neig=1 - type to assign to bonds between neighbors
#     bond_type_diag=2 - type to assign to diagonal bonds
#     x=0, y=0, z=0 - coordinates of platelet's geometrical center
#     start_atom_id=1 - id of the first atom
#     start_bond_id=1 - id of the first bond
#     start_angle_id=1 - id of the first angle
#     structure - if defined, MMT is added into the structure, otherwise
#         new structure is created
#     charged_count=0 - number of MMT beads carrying charges
#     bead_charge=0 - charge of every charged bead


import random


def circular(**kwargs):
    # Define required parameters
    platelet_radius = kwargs['platelet_radius']
    bead_radius = kwargs['bead_radius']
    atom_type = kwargs['atom_type'] if 'atom_type' in kwargs.keys() else 1
    # Define optional parameters
    try:
        bond_type_edge = kwargs['bond_type_edge']
    except KeyError:
        bond_type_edge = 1
    try:
        bond_type_diagonal = kwargs['bond_type_diagonal']
    except KeyError:
        bond_type_diagonal = bond_type_edge + 1
    x = kwargs['x'] if 'x' in kwargs.keys() else 0
    y = kwargs['y'] if 'y' in kwargs.keys() else 0
    z = kwargs['z'] if 'z' in kwargs.keys() else 0
    try:
        atom_id = kwargs['start_atom_id']
    except KeyError:
        atom_id = 1
    try:
        bond_id = kwargs['start_bond_id']
    except KeyError:
        bond_id = 1
    try:
        structure = kwargs['structure']
    except KeyError:
        structure = {'atoms': {}, 'bonds': {}}
    try:
        charged_count = kwargs['charged_count']
        bead_charge = kwargs['bead_charge']
    except KeyError:
        charged_count = 0
    try:
        atom_args = kwargs['atom_args']
    except KeyError:
        atom_args = {}

    atom_id_start = atom_id
    atom_ids = {}
    # Add atoms
    for idxx in range(-platelet_radius+1, platelet_radius):
        atom_ids[idxx] = {}
        dx = idxx * 2*bead_radius
        for idxy in range(-platelet_radius+1, platelet_radius):
            if platelet_radius**2 < idxx**2 + idxy**2:
                continue
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
    for k in range(atom_id_start, atom_id):
        for k1, v1 in atom_args.items():
            structure['atoms'][k][k1] = v1
    # Add bonds
    for idxx in atom_ids.keys():
        for idxy in atom_ids[idxx].keys():
            if idxx - 1 in atom_ids.keys() and idxy in atom_ids[idxx - 1].keys():
                structure['bonds'][bond_id] = {
                    'atoms': [atom_ids[idxx - 1][idxy]['top'],
                              atom_ids[idxx][idxy]['top']],
                    'type': bond_type_edge
                }
                bond_id += 1
                structure['bonds'][bond_id] = {
                    'atoms': [atom_ids[idxx - 1][idxy]['bottom'],
                              atom_ids[idxx][idxy]['bottom']],
                    'type': bond_type_edge
                }
                bond_id += 1
            if idxy - 1 in atom_ids[idxx].keys():
                structure['bonds'][bond_id] = {
                    'atoms': [atom_ids[idxx][idxy - 1]['top'],
                              atom_ids[idxx][idxy]['top']],
                    'type': bond_type_edge
                }
                bond_id += 1
                structure['bonds'][bond_id] = {
                    'atoms': [atom_ids[idxx][idxy - 1]['bottom'],
                              atom_ids[idxx][idxy]['bottom']],
                    'type': bond_type_edge
                }
                bond_id += 1
            # digonal
            if (idxx - 1 in atom_ids.keys()
                and idxy - 1 in atom_ids[idxx - 1].keys()):
                structure['bonds'][bond_id] = {
                    'atoms': [atom_ids[idxx - 1][idxy - 1]['top'],
                              atom_ids[idxx][idxy]['bottom']],
                    'type': bond_type_diagonal
                }
                bond_id += 1
            if (idxx - 1 in atom_ids.keys()
                and idxy + 1 in atom_ids[idxx - 1].keys()):
                structure['bonds'][bond_id] = {
                    'atoms': [atom_ids[idxx - 1][idxy + 1]['top'],
                              atom_ids[idxx][idxy]['bottom']],
                    'type': bond_type_diagonal
                }
                bond_id += 1
            if (idxx + 1 in atom_ids.keys()
                and idxy - 1 in atom_ids[idxx + 1].keys()):
                structure['bonds'][bond_id] = {
                    'atoms': [atom_ids[idxx + 1][idxy - 1]['top'],
                              atom_ids[idxx][idxy]['bottom']],
                    'type': bond_type_diagonal
                }
                bond_id += 1
            if (idxx + 1 in atom_ids.keys()
                and idxy + 1 in atom_ids[idxx + 1].keys()):
                structure['bonds'][bond_id] = {
                    'atoms': [atom_ids[idxx + 1][idxy + 1]['top'],
                              atom_ids[idxx][idxy]['bottom']],
                    'type': bond_type_diagonal
                }
                bond_id += 1
    # Assign (or not) charges
    if charged_count == 0:
        return True, structure
    if atom_id - atom_id_start < charged_count:
        print('Error: cannot charge {0} atoms since only {1} are added'.format(
            charged_count, atom_id - atom_id_start))
        return False, {}
    charged_ids = set()
    while len(charged_ids) < charged_count:
        charged_ids.add(random.choice(range(atom_id_start, atom_id)))
    for atom_id in charged_ids:
        structure['atoms'][atom_id]['charge'] = bead_charge
    return True, structure
