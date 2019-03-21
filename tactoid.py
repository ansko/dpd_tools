import copy
import random


from structure.mmt.circular import circular
from structure.modifier.between_clays import between_clays
from structure.polymer.random_chain import random_chain
from other.write_data import write_data


def make_tactoid():
    '''
    Make circular mmt platelet in the parallelepiped box surrounded by
    the charged modifiers and uncharged polymer in the cubic box.
    '''
    # CEC: 92.6 milliequiv / 100g
    # mass of mmt: 24*0 + 8*Si + 4*H + 4*Al + 1*Na = 720amu + 23(Na)amu [* 9]
    # 100 * 1/720 * 2/3 = 0.0926
    # -> CEC = 92.6 milliequiv / 100g
    # mmt surface denstiy: 720amu / (90*80 A2) ~ 0.1 amu/A2

    # From MD:
    lx = 90; ly = 80; lz = 65
    mmt_thickness = 8
    md_soft_atoms = 7560 + 34380

    # General parameters:
    platelet_radius = 6   # radius of mmt lamella (in beads)
    tail_length = 2       # length of the tail of the modifiers
    bead_charge = 1       # charge of every modifier and substitution in mmt
    stacking = 2          # mmt lamellas in tactoid
    polymerization = 100  # length of every polymer
    real_interlayer = 8  # interlayer emptiness thickness

    # L system in MD:
    #     Size: 90A x 80A x 65A
    #     Atoms: 48420
    #         6480  MMT
    #         7560  modifier
    #         34380 polymer
    #     Soft bead size: 90*80*(65-8) / (7560 + 34380) ~ 9.8 A3
    # DPD:
    #     density = 3
    #     DPD bead radius [lj]: sqrt3(density * soft_bead_volume) ~ 3.1A

    # r_c = (rho * 4/3 pi r**3)**(1/3)
    # r_c (1 / rho * 3/4/pi)**(1/3) = r
    # r = 0.43 * r_c
    dpd_rho = 3
    real_bead_radius = (lx * ly * (lz - mmt_thickness) / md_soft_atoms /
                        (4/3*3.14))**(1/3)
    real_cutoff = (dpd_rho * 4/3 * 3.14)**(1/3) * real_bead_radius
    lj_bead_radius = 1
    lj_interlayer = real_interlayer * lj_bead_radius/real_bead_radius

    # So, in the L system there were 108 substitutions per 90*80 A2,
    # i.e. they had densty 0.015 substitution/A2
    # Now MMT as size: platelet_radius**2 * 3.1**2 A2
    # Get modifiers count per one lamella:
    real_mmt_area = 3.14 * platelet_radius**2 * real_bead_radius**2
    print(real_mmt_area, real_mmt_area/7200)
    charged_count = 0.015 * real_mmt_area
    charged_count = int(round(charged_count, 0)) * stacking
    print('charged count is', charged_count)
    # Compute box size
    min_height = real_interlayer*stacking + 4*real_bead_radius * stacking
    cube_edge = 2 * platelet_radius * real_bead_radius
    cube_edge = max(cube_edge, min_height) * lj_bead_radius/real_bead_radius

    # Now adjust polymers count:
    free_volume = ((cube_edge*real_bead_radius - 8)**3
                   -charged_count * (1 + tail_length))
    polymer_volume = polymerization * real_bead_radius**3
    polymers_count = int(round(dpd_rho * free_volume / polymer_volume, 0))
    print('should be {0} polymers'.format(polymers_count))

    # From now all units are only lj
    # Create box:
    structure = {
        'cell': {
            'xlo': -cube_edge, 'xhi': cube_edge,
            'ylo': -cube_edge, 'yhi': cube_edge,
            'zlo': -cube_edge, 'zhi': cube_edge
        },
        'atoms': {}, 'bonds': {}
    }
    # Add MMT:
    mmt_common_args = {
        'platelet_radius': platelet_radius,
        'bead_radius': lj_bead_radius,
        'atom_type': 1,
        'bond_type_edge': 1,
        'bond_type_diagonal': 2,
        'charged_count': charged_count,
    }
    half_stacking = stacking // 2
    atom_id = 1
    bond_id = 1
    if stacking % 2 == 0:
        dz = lj_interlayer/2 + 2*lj_bead_radius
    else:
        status, new_structure = circular(**mmt_common_args,
            bead_charge=-bead_charge,
            structure=structure)
        if status:
            structure = new_structure
        dz = lj_interlayer + 4*lj_bead_radius
        atom_id = len(structure['atoms']) + 1
        bond_id = len(structure['bonds']) + 1
    for idx in range(half_stacking):
        # To the top of the tactoid
        status, new_structure = circular(**mmt_common_args,
            start_atom_id=atom_id,
            start_bond_id=bond_id,
            z=dz,
            bead_charge=-bead_charge,
            structure=structure)
        if status:
            structure = new_structure
            atom_id = len(structure['atoms']) + 1
            bond_id = len(structure['bonds']) + 1
        # To the bottom of the tactoid
        status, new_structure = circular(**mmt_common_args,
            start_atom_id=atom_id,
            start_bond_id=bond_id,
            z=-dz,
            bead_charge=-bead_charge,
            structure=structure)
        if status:
            structure = new_structure
            atom_id = len(structure['atoms']) + 1
            bond_id = len(structure['bonds']) + 1
            dz += lj_interlayer + 4*lj_bead_radius
    mmt_atoms_count = len(structure['atoms'])

    top = lj_interlayer/2
    bottom = -lj_interlayer/2
    #print('atoms in mmt', mmt_atoms_count)
    modifiers_done = 0
    fails_done = 0
    fails_allowed = charged_count * 10
    while modifiers_done < charged_count and fails_done < fails_allowed:
        status, new_structure = between_clays(
            top=top, bottom=bottom,
            structure=structure,
            head_charge=bead_charge,
            tail_length=tail_length,
            bead_radius=lj_bead_radius,
            head_atom_type=2,
            tail_atom_type=3,
            head_tail_bond_type=3,
            tail_tail_bond_type=4)
        if status:
            modifiers_done += 1
            structure = new_structure
        else:
            fails_done += 1
    print('charged', len(structure['atoms']))

    polymers_done = 0
    polymers_fails_done = 0
    polymers_fails_allowed = 100 * polymers_count
    while (polymers_done < polymers_count and
        polymers_fails_done < polymers_fails_allowed):
        old_structure = copy.deepcopy(structure)
        if polymerization > 1:
            status, new_structure = random_chain(
                polymerization=polymerization,
                bead_radius=lj_bead_radius,
                atom_type=4,
                bond_type=5,
                start_atom_id=len(structure['atoms']) + 1,
                start_bond_id=len(structure['bonds']) + 1,
                structure=structure)
        else:
            status, new_structure = polymer_chain(
                polymerization=polymerization,
                bead_radius=bead_radius,
                atom_type=4,
                start_atom_id=len(structure['atoms']) + 1,
                structure=structure)
        if status:
            structure = new_structure
            polymers_done += 1
            print('polymers:', polymers_done)
        else:
            structure = copy.deepcopy(old_structure)
            polymers_fails_done += 1

    write_data('tactoid.data', structure)


if __name__ == '__main__':
    make_tactoid()
