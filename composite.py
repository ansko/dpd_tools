import copy
import pprint
pprint = pprint.PrettyPrinter(indent=4).pprint

from basic.mmt import mmt_sheet_double_diagonals_circular
from basic.modifier import modifier_random
from basic.polymer import polymer_chain
from basic.get_structure_params import get_structure_params_L
from basic.write_data import write_data


def make_ternary_composite_no_angles_from_structure_params(system=None):
    if system == 'L':
        parameters = get_structure_params_L()
        bead_radius = parameters['bead_radius']
        platelet_radius = parameters['platelet_radius']
        head_charge = parameters['head_charge']
        tail_length = parameters['tail_length']
        polymerization = parameters['polymerization']
        modifiers_count = parameters['modifiers_count']
        polymers_count = parameters['polymers_count']
        lx = parameters['lx']
        ly = parameters['ly']
        lz = parameters['lz']
    else:
        # Default / test configuration
        bead_radius = 1       # bead radius - same for all beads
        platelet_radius = 10  # mmt platelet's radius (in beads)
        head_charge = 1       # charge of the head of the modifier
        tail_length = 2       # length of modifier non-charged tail
        polymerization = 100  # length of polymer chain (in beads)
        modifiers_count = 10  # number of modifier molecules in the system
        polymers_count = 10   # number of polymer chain in the system
        lx = 2*platelet_radius * 2*bead_radius * 2
        ly = lx
        lz = lx

    atom_types = {
        'mmt': 1,
        'modifier_head': 2
    }
    bond_types = {
        'mmt_edge': 1,
        'mmt_diagonal': 2,
    }
    angle_types = {
        'mmt': 1
    }
    if tail_length > 1:
        atom_types['modifier_tail'] = 3
        atom_types['polymer'] = 4
        bond_types['modifier_head_tail'] = 3
        bond_types['modifier_tail_tail'] = 4
        bond_types['polymer'] = 5
    elif tail_length == 1:
        atom_types['modifier_tail'] = 3
        atom_types['polymer'] = 4
        bond_types['modifier_head_tail'] = 3
        bond_types['polymer'] = 4
    else:
        atom_types['polymer'] = 3
        bond_types['polymer'] = 3
    # Box
    structure = {
        'atoms': {}, 'bonds': {}, 'angles': {},
        'cell': {
            'xlo': -lx/2,
            'xhi': lx/2,
            'ylo': -ly/2,
            'yhi': ly/2,
            'zlo': -lz/2,
            'zhi': lz/2
        }
    }
    # 1st MMT platelet
    structure = mmt_sheet_double_diagonals_circular(
        platelet_radius=platelet_radius, bead_radius=bead_radius,
        x=0, y=0, z=-3*bead_radius,
        start_atom_id=1,
        start_bond_id=1,
        start_angle_id=1,
        atom_type = atom_types['mmt'],
        bond_type_edge = bond_types['mmt_edge'],
        bond_type_diagonal = bond_types['mmt_diagonal'],
        angle_type = angle_types['mmt'],
        charge=modifiers_count*head_charge/2,
        structure=structure)
    atoms_in_platelet = len(structure['atoms'])
    bonds_in_platelet = len(structure['bonds'])
    angles_in_platelet = len(structure['angles'])
    # 2nd MMT platelet
    structure = mmt_sheet_double_diagonals_circular(
        platelet_radius=platelet_radius, bead_radius=bead_radius,
        x=0, y=0, z=3*bead_radius,
        start_atom_id=atoms_in_platelet + 1,
        start_bond_id=bonds_in_platelet + 1,
        start_angle_id=angles_in_platelet + 1,
        atom_type = atom_types['mmt'],
        bond_type_edge = bond_types['mmt_edge'],
        bond_type_diagonal = bond_types['mmt_diagonal'],
        angle_type = angle_types['mmt'],
        charge=modifiers_count*head_charge/2,
        structure=structure)
    assert(atoms_in_platelet * 2 == len(structure['atoms']))
    assert(bonds_in_platelet * 2 == len(structure['bonds']))
    assert(angles_in_platelet * 2 == len(structure['angles']))
    mmt_atoms = atoms_in_platelet * 2
    mmt_bonds = bonds_in_platelet * 2
    mmt_angles = angles_in_platelet * 2
    # modifiers
    modifiers_done = 0
    modifiers_fails_done = 0
    modifiers_fails_allowed = modifiers_count * 10
    while (modifiers_done < modifiers_count and
        modifiers_fails_done < modifiers_fails_allowed):
        old_structure = copy.deepcopy(structure)
        if tail_length > 1:
            status, new_structure = modifier_random(
                head_charge=head_charge,
                tail_length=tail_length,
                bead_radius=bead_radius,
                start_atom_id=len(structure['atoms']) + 1,
                start_bond_id=len(structure['bonds']) + 1,
                head_atom_type=atom_types['modifier_head'],
                tail_atom_type=atom_types['modifier_tail'],
                head_tail_bond_type=bond_types['modifier_head_tail'],
                tail_tail_bond_type=bond_types['modifier_tail_tail'],
                structure=structure)
        elif tail_length == 1:
            status, new_structure = modifier_random(
                head_charge=head_charge,
                tail_length=tail_length,
                bead_radius=bead_radius,
                start_atom_id=len(structure['atoms']) + 1,
                start_bond_id=len(structure['bonds']) + 1,
                head_atom_type=atom_types['modifier_head'],
                tail_atom_type=atom_types['modifier_tail'],
                head_tail_bond_type=bond_types['modifier_head_tail'],
                structure=structure)
        else:
            status, new_structure = modifier_random(
                head_charge=head_charge,
                tail_length=tail_length,
                bead_radius=bead_radius,
                atom_id=len(structure['atoms']) + 1,
                bond_id=len(structure['bonds']) + 1,
                head_atom_type=atom_types['modifier_head'],
                structure=structure)
        if status:
            modifiers_done += 1
            old_structure = new_structure
            print('modifiers:', modifiers_done)
        else:
            structure = copy.deepcopy(old_structure)
            modifiers_fails_done += 1
    modifier_atoms = len(structure['atoms']) - mmt_atoms
    modifier_bonds = len(structure['bonds']) - mmt_bonds
    assert(modifier_atoms == modifiers_count * (1 + tail_length))
    assert(modifier_bonds == modifiers_count * tail_length)
    '''
    # polymer
    polymers_done = 0
    polymers_fails_done = 0
    polymers_fails_allowed = polymers_count * polymerization
    while (polymers_done < polymers_count and
        polymers_fails_done < polymers_fails_allowed):
        old_structure = copy.deepcopy(structure)
        if polymerization > 1:
            status, new_structure = polymer_chain(
                polymerization=polymerization,
                bead_radius=bead_radius,
                atom_type=atom_types['polymer'],
                bond_type=bond_types['polymer'],
                start_atom_id=len(structure['atoms']) + 1,
                start_bond_id=len(structure['bonds']) + 1,
                structure=structure)
        else:
            status, new_structure = polymer_chain(
                polymerization=polymerization,
                bead_radius=bead_radius,
                atom_type=atom_types['polymer'],
                start_atom_id=len(structure['atoms']) + 1,
                structure=structure)
        if status:
            structure = new_structure
            polymers_done += 1
            print('polymers:', polymers_done)
        else:
            structure = copy.deepcopy(old_structure)
            polymers_fails_done += 1
    polymer_atoms = len(structure['atoms']) - mmt_atoms - modifier_atoms
    polymer_bonds = len(structure['bonds']) - mmt_bonds - modifier_bonds
    assert(polymer_atoms == polymers_count * polymerization)
    assert(polymer_bonds == polymers_count * (polymerization - 1))
    structure_description = {
        'atom_types': atom_types,
        'bond_types': bond_types,
        'angle_types': angle_types,
        'atoms_count': {
            'mmt': mmt_atoms,
            'modifier': modifier_atoms,
            'polymer': polymer_atoms,
            'all': mmt_atoms + modifier_atoms + polymer_atoms
        },
        'bonds_count': {
            'mmt': mmt_bonds,
            'modifier': modifier_bonds,
            'polymer': polymer_bonds,
            'all': mmt_bonds + modifier_bonds + polymer_bonds
        },
        'angles_count': {'mmt': mmt_angles},
        'phases': {
            'mmt': {
                'platelets': 2,
                'atoms_in_platelet': atoms_in_platelet
            },
            'modifier': {
                'molecules_count': modifiers_count,
                'head_charge': head_charge,
                'tail_length': tail_length
            },
            'polymer': {
                'molecules_count': polymers_count,
                'polymerization': polymerization
            }
        }
    }
    '''
    print('mmt atoms', atoms_in_platelet)
    write_data('composite.data', structure)


if __name__ == '__main__':
    make_ternary_composite_no_angles_from_structure_params(system='L')
