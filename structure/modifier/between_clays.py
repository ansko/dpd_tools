# The function is made to add modifier into the interlayer
# (having value of z coordintate between top and bottom).
#
# Required parameters:
#     head_charge - charge of the head of the modifier
#     tail_length - length (in beads) of uncharged tail (>=0)
#     bead_raduis - bead radius in DPD
#     structure - structure where charged filler already exists
#     top - top z value of space for modifiers
#     bottom - bottom z value of space for modifiers
# Optional parameters:
#     head_atom_type=1 - type to assign to head atom
#     tail_atom_type=1 - type to assign to all tail atoms
#     start_atom_id=1 - id of the first atom
#     start_bond_id=1 - id of the first bond
#     head_tail_bond_type=1 - type of bond between head and tail
#     tail_tail_bond_type=2 - type of bond between tail atoms1


import math
import random


def between_clays(**kwargs):
    # Define required parameters
    head_charge = kwargs['head_charge']
    tail_length = kwargs['tail_length']
    bead_radius = kwargs['bead_radius']
    top = kwargs['top']
    bottom = kwargs['bottom']
    xlo = kwargs['structure']['cell']['xlo']
    xhi = kwargs['structure']['cell']['xhi']
    ylo = kwargs['structure']['cell']['ylo']
    yhi = kwargs['structure']['cell']['yhi']
    zlo = kwargs['structure']['cell']['zlo']
    zhi = kwargs['structure']['cell']['zhi']
    lx = xhi - xlo
    ly = yhi - ylo
    lz = zhi - zlo
    # Define optional parameters
    try:
        structure = kwargs['structure']
    except KeyError:
        structure = {'atoms': {}, 'bonds': {}}
    try:
        head_type = kwargs['head_atom_type']
        tail_type = kwargs['tail_atom_type']
    except KeyError:
        try:
            head_type = max([v['type'] for v in structure['atoms'].values()])
        except ValueError:
            head_type = 1
        tail_type = head_type + 1
    try:
        atom_id = kwargs['start_atom_id']
    except KeyError:
        try:
            atom_id = max(structure['atoms'].keys()) + 1
        except ValueError:
            atom_id = 1
    try:
        bond_id = kwargs['start_bond_id']
    except KeyError:
        try:
            bond_id = max(structure['bonds'].keys()) + 1
        except ValueError:
            bond_id = 1
    if tail_length > 0:
        try:
            ht_bond_type = kwargs['head_tail_bond_type']
            tt_bond_type = kwargs['tail_tail_bond_type']
        except KeyError:
            try:
                ht_bond_type = max([
                    v['type'] for v in structure['atoms'].values()])
                tt_bond_type = ht_bond_type + 1
            except ValueError:
                ht_bond_type = 1
                tt_bond_type = 2

    # Start algorithm
    charged_clay_ids = set()
    for k, v in structure['atoms'].items():
        if v['phase'] == 'filler' and v['charge'] != 0:
            charged_clay_ids.add(k)
    charged_clay_ids = list(charged_clay_ids)
    in_chain_id = 1
    modifier = {}
    fails_allowed = 100
    fails_done = 0
    while fails_done < fails_allowed:
        if not modifier:
        #    clay_bead = structure['atoms'][random.choice(charged_clay_ids)]
        #    x = clay_bead['x']
        #    y = clay_bead['y']
        #    z = clay_bead['z']
            x = xlo + lx/4 + random.random() * lx/2
            y = ylo + ly/4 + random.random() * ly/2
            z = bottom + (top - bottom) * random.random()
        else:
            x = modifier[max(modifier.keys())]['x']
            y = modifier[max(modifier.keys())]['y']
            z = modifier[max(modifier.keys())]['z']
        theta = random.random() * math.pi
        phi = random.random() * 2*math.pi
        x += 2.25*bead_radius * math.sin(theta) * math.cos(phi)
        y += 2.25*bead_radius * math.sin(theta) * math.sin(phi)
        z += 2.25*bead_radius * math.cos(theta)
        all_is_ok = True
        for atom in structure['atoms'].values():
            dx = x - atom['x']
            dy = y - atom['y']
            dz = z - atom['z']
            dr2 = dx**2 + dy**2 + dz**2
            if dr2 < 4*bead_radius**2:
                all_is_ok = False
                break
        if not all_is_ok:
            fails_done += 1
            continue
        for atom in modifier.values():
            dx = x - atom['x']
            dy = y - atom['y']
            dz = z - atom['z']
            dr2 = dx**2 + dy**2 + dz**2
            if dr2 < 4*bead_radius**2:
                all_is_ok = False
                break
        if not all_is_ok:
            fails_done += 1
            continue
        if not modifier:  # now appending head atom
            charge = head_charge
        else:
            charge = 0
        modifier[in_chain_id] = {
            'type': head_type if not modifier else tail_type,
            'charge': charge,
            'x': x if xlo < x < xhi else x - lx if x > xhi else x + lx,
            'y': y if ylo < y < yhi else y - ly if y > yhi else y + ly,
            'z': z if zlo < z < zhi else z - lz if z > zhi else z + lz,
            'phase': 'modifier'
        }
        in_chain_id += 1
        if len(modifier) == tail_length + 1:
            break
    if len(modifier) != tail_length + 1:
        return False, {}
    for k, v in modifier.items():
        structure['atoms'][atom_id + k - 1] = v
        if k == 2:
            structure['bonds'][bond_id + k - 2] = {
                'type': ht_bond_type,
                'atoms': [atom_id + k - 2, atom_id + k - 1]
            }
        elif k > 2:
            structure['bonds'][bond_id + k - 2] = {
                'type': tt_bond_type,
                'atoms': [atom_id + k - 2, atom_id + k - 1]
            }
    return True, structure
