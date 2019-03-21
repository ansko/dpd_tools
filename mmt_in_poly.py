#!/usr/bin/env python3


# Make MMT platelet comprised from circular lamellae in the
# parallelepiped box saturated by the charged ions and surrounded
# by uncharged polymer in the cubic box.


import copy
import random


from structure.mmt.circular import circular
from structure.modifier.between_clays import between_clays
from structure.polymer.random_chain import random_chain
from other.write_data import write_data


def make_tactoid():
    # General system parameters:
    platelet_radius = 4   # radius of mmt lamella [in beads]
    tail_length = 0       # OMMT has Na cations in interlayer
    bead_charge = 0.01       # charge of every cation and substitution in MMT
    stacking = 2          # MMT lamellas in tactoid
    polymerization = 100  # length of every polymer [in monomers/beads]
    real_interlayer = 4   # MMT interlayer emptiness
    lj_bead_radius = 1
    coeff = 3             # ratio between box xy size and platelet radius

    # MMT CEC == 92.6 milliequiv / 100g
    # 2 substituaion per 3 * 720 amu: 100 * 1/720 * 2/3 = 0.0926 equiv / 100g
    # mmt surface denstiy: 18*9*720amu / (90*80 A2) ~ 16.2 amu/A2
    mmt_surf_density = 16.2  # single MMT layer density [amu/A2]
    mmt_cec = 20#93  # cation exchange capacity of MMT [mequiv/100g]

    # Cell parameters from MD:
    lx = 90; ly = 80; lz = 65
    mmt_thickness = 8  # space where soft phase atoms are absent
    soft_beads = 3*70 + 20*90

    # Compute DPD bead radius:
    # (i do not completely understand the concept)
    # real_cutoff is further used as lj distance unit
    dpd_rho = 3  # used in literature
    real_bead_radius = (lx * ly * (lz - mmt_thickness) / soft_beads /
                        (4/3*3.14))**(1/3)
    #real_cutoff = (dpd_rho * 4/3 * 3.14)**(1/3) * real_bead_radius
    real_cutoff = (dpd_rho * lx * ly * (lz - mmt_thickness) / soft_beads)**(1/3)
    lj_interlayer = real_interlayer * lj_bead_radius/real_bead_radius
    print('Real cutoff and bead radius:', real_cutoff, real_bead_radius)
    print(real_bead_radius / real_cutoff)  # 0.43

    # Na cation size:
    real_Na_radius = 2.27
    lj_Na_radius = real_Na_radius / real_cutoff
    print('vdW Na radius:', lj_Na_radius)

    # Get cations count:
    real_mmt_area = 3.14 * platelet_radius**2 * real_bead_radius**2
    print('area', real_mmt_area)
    charged_count = mmt_cec / 1e5 * (real_mmt_area*mmt_surf_density) * stacking
    charged_count = int(round(charged_count, 0))
    print('Modifers/substitions:', charged_count)

    # Compute box size
    min_height = real_interlayer*(stacking + 1) + 4*real_bead_radius * stacking
    cube_edge = coeff * platelet_radius * 2*real_bead_radius
    cube_edge = max(cube_edge, min_height) * lj_bead_radius/real_bead_radius

    # Adjust polymers count according to space left in the box:
    free_volume = ((cube_edge*real_bead_radius - 8)**3
                   -charged_count * (1 + tail_length))
    polymer_volume = polymerization * 4/3*3.14 * real_bead_radius**3
    #polymers_count = int(round(dpd_rho * free_volume / polymer_volume, 0))
    polymers_count = int(round(free_volume / polymer_volume, 0))
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
    # 1 - Add lamellae:
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
    print('mmt_mass', real_mmt_area/40 * 720 / mmt_atoms_count)

    # 2 - Add interlayer cations:
    cation_places = stacking + 1  # count of places for cations: interlayers + 2
    cations_in_place = int(round(charged_count/cation_places, 0))
    for idx in range(half_stacking):
        # top
        modifiers_done = 0
        fails_done = 0
        fails_allowed = cations_in_place * 100
        while fails_done < fails_allowed and modifiers_done < cations_in_place:
            if stacking % 2 == 0:
                top = lj_interlayer/2 + idx * (lj_interlayer + 4*lj_bead_radius)
            else:
                top = lj_interlayer * idx + 4*lj_bead_radius * (idx + 0.5)
            bottom = top -lj_interlayer
            status, new_structure = between_clays(
                top=top,
                bottom=bottom,
                structure=structure,
                head_charge=bead_charge,
                tail_length=tail_length,
                bead_radius=lj_bead_radius,
                head_atom_type=2,
                tail_atom_type=None)
            if status:
                modifiers_done += 1
                structure = new_structure
            else:
                fails_done += 1
        # bottom
        modifiers_done = 0
        fails_done = 0
        if stacking % 2 == 0 and idx == 0:
            continue
        while fails_done < fails_allowed and modifiers_done < cations_in_place:
            if stacking % 2 == 0:
                top = lj_interlayer/2 - idx * (lj_interlayer + 4*lj_bead_radius)
            else:
                top = lj_interlayer * (-idx + 1) + 4*lj_bead_radius * (-idx + 0.5)
            bottom = top -lj_interlayer
            status, new_structure = between_clays(
                top=top,
                bottom=bottom,
                structure=structure,
                head_charge=bead_charge,
                tail_length=tail_length,
                bead_radius=lj_bead_radius,
                head_atom_type=2,
                tail_atom_type=None)
            if status:
                modifiers_done += 1
                structure = new_structure
            else:
                fails_done += 1
    # from top and bottom
    if stacking % 2 == 0:
        top = (lj_interlayer/2
               + half_stacking * lj_interlayer
               + half_stacking * 4*lj_bead_radius)
    else:
        top = (2*lj_bead_radius + lj_interlayer 
               + half_stacking * lj_interlayer
               + half_stacking * 4*lj_bead_radius)
    bottom = top - lj_interlayer/4
    modifiers_done = 0
    fails_done = 0
    while fails_done < fails_allowed and modifiers_done < cations_in_place:
        status, new_structure = between_clays(
            top=top,
            bottom=bottom,
            structure=structure,
            head_charge=bead_charge,
            tail_length=tail_length,
            bead_radius=lj_bead_radius,
            head_atom_type=2,
            tail_atom_type=None)
        if status:
            modifiers_done += 1
            structure = new_structure
        else:
            fails_done += 1
    if stacking % 2 == 0:
        top = (-lj_interlayer/2
               - half_stacking * lj_interlayer
               - half_stacking * 4*lj_bead_radius)
    else:
        top = (-2*lj_bead_radius + lj_interlayer 
               - half_stacking * lj_interlayer
               - half_stacking * 4*lj_bead_radius)
    bottom = top - lj_interlayer/4
    modifiers_done = 0
    fails_done = 0
    while fails_done < fails_allowed and modifiers_done < cations_in_place:
        status, new_structure = between_clays(
            top=top,
            bottom=bottom,
            structure=structure,
            head_charge=bead_charge,
            tail_length=tail_length,
            bead_radius=lj_bead_radius,
            head_atom_type=2,
            tail_atom_type=None)
        if status:
            modifiers_done += 1
            structure = new_structure
        else:
            fails_done += 1
    print('charged', len(structure['atoms']))

    # Add polymers:
    polymers_done = 0
    polymers_fails_done = 0
    polymers_fails_allowed = 10# * polymers_count
    while (polymers_done < polymers_count and
        polymers_fails_done < polymers_fails_allowed):
        old_structure = copy.deepcopy(structure)
        if polymerization > 1:
            status, new_structure = random_chain(
                polymerization=polymerization,
                bead_radius=lj_bead_radius,
                atom_type=3,
                bond_type=3,
                start_atom_id=len(structure['atoms']) + 1,
                start_bond_id=len(structure['bonds']) + 1,
                structure=structure)
        else:
            status, new_structure = polymer_chain(
                polymerization=polymerization,
                bead_radius=bead_radius,
                atom_type=3,
                start_atom_id=len(structure['atoms']) + 1,
                structure=structure)
        if status:
            structure = new_structure
            polymers_done += 1
            print('polymers:', polymers_done)
        else:
            structure = copy.deepcopy(old_structure)
            polymers_fails_done += 1

    write_data('mmt_in_poly.data', structure)


if __name__ == '__main__':
    make_tactoid()
