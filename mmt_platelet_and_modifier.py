from basic.mmt import mmt_sheet_double_diagonals_periodic
from basic.modifier import modifier_near_clays
from basic.fill_by_beads import fill_by_beads
from basic.write_data import write_data
from basic.polymer import polymer_chain


def make_mmt_platelet_and_modifier():
    '''
    Make periodic mmt platelet in parallelepiped box
    surrounded by charged modifiers.
    Vacuum is filled with uncharged beads.
    '''
    bead_radius = 1       # bead radius - same for all beads
    cell_side = 10        # horizontal cell size along x and y (in beads)
    cell_height = 15      # vertical cell side (in beads)
    tail_length = 2       # length of the tail of the modifiers

    charged_count = 10
    bead_charge = 1

    structure = mmt_sheet_double_diagonals_periodic(
        cell_side=cell_side,
        bead_radius=bead_radius,
        atom_type=1,
        bond_type_edge=1,
        bond_type_diagonal=2,
        charged_count=charged_count,
        bead_charge=bead_charge)
    structure['cell'] = {
        'xlo': 0, 'xhi': cell_side * 2*bead_radius,
        'ylo': 0, 'yhi': cell_side * 2*bead_radius,
        'zlo': -cell_height/2*bead_radius, 'zhi': cell_height/2 * bead_radius
    }
    mmt_atoms_count = len(structure['atoms'])
    print('atoms in mmt', mmt_atoms_count)

    modifiers_done = 0
    fails_done = 0
    fails_allowed = charged_count * 10
    while modifiers_done < charged_count and fails_done < fails_allowed:
        status, structure = modifier_near_clays(
            structure=structure,
            head_charge=bead_charge,
            tail_length=tail_length,
            bead_radius=bead_radius,
            head_atom_type=2,
            tail_atom_type=3,
            head_tail_bond_type=3,
            tail_tail_bond_type=4)
        if status:
            modifiers_done += 1
        else:
            fails_done += 1

    print('atoms in charged modifier', len(structure['atoms']) - mmt_atoms_count)
    print('group charged id 1:{0}'.format(len(structure['atoms'])))
    N_is = len(structure['atoms'])
    N_should_be = cell_side**2 * bead_radius**3 * 4 * cell_height * 3
    N_last = int(N_should_be - N_is) / (4*bead_radius**3)
    structure = fill_by_beads(
        structure=structure,
        bead_radius=bead_radius,
        beads_count=N_last)
    print('all atoms', len(structure['atoms']))

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

    write_data('mmt_and_mod.data', structure)


if __name__ == '__main__':
    make_mmt_platelet_and_modifier()
