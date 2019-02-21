import random


def fill_by_beads(**kwargs):
    '''
    DPD does not treat vacuum properly, this function fills vacuum by beads.
    Required parameters:
        structure - structure where new beads are added
        bead_radius - size of DPD bead
        beads_count - number of beads to add
    Optional parameters
        bead_type=1 - type to assing to all added beads
        charge=0 - charge to assign to all added beads
        bead_charge=0 - charged of every charged bead
        charged_count=0 - count of beads that should be charged
    '''
    structure = kwargs['structure']
    bead_radius = kwargs['bead_radius']
    beads_count = kwargs['beads_count']
    try:
        bead_type = kwargs['bead_type']
    except KeyError:
        try:
            bead_type = max([v['type'] for v in structure['atoms'].values()]) + 1
        except ValueError:
            bead_type = 1
    try:
        bead_charge = kwargs['bead_charge']
        charged_count = kwargs['charged_count']
    except KeyError:
        bead_charge = 0
        charged_count = 0
    xlo = structure['cell']['xlo']
    xhi = structure['cell']['xhi']
    lx = xhi - xlo
    ylo = structure['cell']['ylo']
    yhi = structure['cell']['yhi']
    ly = yhi - ylo
    zlo = structure['cell']['zlo']
    zhi = structure['cell']['zhi']
    lz = zhi - zlo
    N_already_is = len(structure['atoms'])
    new_beads = {}
    atom_id = N_already_is + 1
    N_last = beads_count
    fails_allowed = beads_count * 10
    fails_done = 0
    while N_last > 0 and fails_done < fails_allowed:
        x = xlo + random.random() * lx
        y = ylo + random.random() * ly
        z = zlo + random.random() * lz
        is_close = False
        for atom in structure['atoms'].values():
            dx = atom['x'] - x
            dy = atom['y'] - y
            dz = atom['z'] - z
            if dx**2 + dy**2 + dz**2 < bead_radius**2 * 2:
                is_close = True
                break
        if is_close:
            fails_done += 1
            continue
        new_beads[atom_id] = {
            'type': bead_type, 'charge': 0,
            'x': x, 'y': y, 'z': z
        }
        atom_id += 1
        N_last -= 1
    if N_last != 0:
        print('{0} not appended from {1}'.format(N_last, fails_allowed//10))
    beads_to_charge = set()
    while len(beads_to_charge) < charged_count:
        beads_to_charge.add(random.choice(list(new_beads.keys())))
    for k, v in new_beads.items():
        structure['atoms'][k] = v
        if k in beads_to_charge:
            structure['atoms'][k]['charge'] = bead_charge
    return structure
