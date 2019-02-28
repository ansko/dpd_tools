import pprint
pprint = pprint.PrettyPrinter(indent=4).pprint


def get_structure_params_L():
    '''
    This function returns input parameters for creation of composite
    originating from the L system:
        MMT       82 A x 90 A (6480 atoms == 108 subtitions)
        modifier  108 molecules
        polymer   90 molecules x 20 monomers (382 atoms each)
        ---
        6480 + 108*70 + 90*382 = 48420 atoms
        63 A x 90 A x 82 A ~~ 450000 A**3 volume

    Resulting system:
    2 MMT platelets in box with every side increased 2 times:
        lx = 2 * lx
        ly = 2 * ly
        lz = 2 * lz
    '''
    lx = 90
    ly = 82
    lz = 63
    mmt_thickness = 7
    atoms_count = 48420
    mmt_atoms_count = 6480
    polymerization = 20
    polymers_count = 90
    modifiers_count = 108
    tail_length = 1  # ??? why ???
    head_charge = 1

    polymer_beads = polymerization * polymers_count
    modifier_beads = modifiers_count * (1 + tail_length)
    soft_volume = lx * ly * (lz - mmt_thickness)
    soft_atom_volume = soft_volume / (polymer_beads + modifier_beads)
    soft_bead_radius = (soft_atom_volume / 4)**(1/3)
    platelet_radius = int((lx * ly / soft_bead_radius / soft_bead_radius)**0.5)

    new_volume = 8 * lx*ly*lz
    new_mmt_volume = 3.14 * platelet_radius**2 * mmt_thickness
    new_soft_volume = new_volume - new_mmt_volume
    polymer_volume = new_soft_volume - modifier_beads * 4*soft_bead_radius**2
    polymer_volume /= 4*soft_bead_radius**2
    polymers_count = int(polymer_volume / polymerization)

    parameters = {
        'bead_radius': soft_bead_radius,
        'platelet_radius': platelet_radius,
        'head_charge': head_charge,
        'tail_length': tail_length,
        'polymerization': polymerization,
        'modifiers_count': modifiers_count,
        'polymers_count': polymers_count,
        'lx': lx * 2,
        'ly': ly * 2,
        'lz': lz * 2,
    }
    pprint(parameters)

    return parameters
