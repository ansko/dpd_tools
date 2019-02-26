from basic.mmt import mmt_sheet_double_diagonals_periodic
from basic.fill_by_beads import fill_by_beads
from basic.write_data import write_data


def make_mmt_platelet_uncharged():
    '''
    Make periodic mmt platelet in parallelepiped box.
    Vacuum is filled with beads.
    '''
    bead_radius = 1       # bead radius - same for all beads
    cell_side = 25        # horizontal cell size along x and y (in beads)
    cell_height = 10      # vertical cell side (in beads)
    structure = mmt_sheet_double_diagonals_periodic(
        cell_side=cell_side,
        bead_radius=bead_radius,
        atom_type=1,
        bond_type_edge=1,
        bond_type_diagonal=2,
        charge=0)
    structure['cell'] = {
        'xlo': 0, 'xhi': cell_side * 2*bead_radius,
        'ylo': 0, 'yhi': cell_side * 2*bead_radius,
        'zlo': -cell_height/2*bead_radius, 'zhi': cell_height/2 * bead_radius
    }
    structure = fill_by_beads(
        structure=structure,
        bead_radius=bead_radius,
        fraction=0.75)

    write_data('mmt.data', structure)


def make_mmt_platelet_charged():
    '''
    Make periodic mmt platelet in parallelepiped box.
    Vacuum is filled with beads.
    '''
    bead_radius = 1       # bead radius - same for all beads
    cell_side = 15        # horizontal cell size along x and y (in beads)
    cell_height = 15      # vertical cell side (in beads)

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
    structure = fill_by_beads(
        structure=structure,
        bead_radius=bead_radius,
        beads_count=charged_count,
        charged_count=charged_count,
        bead_charge=-bead_charge)
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
    write_data('mmt.data', structure)


if __name__ == '__main__':
    #make_mmt_platelet_uncharged()
    make_mmt_platelet_charged()
