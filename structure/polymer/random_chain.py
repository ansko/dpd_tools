import math
import random

import pprint
pprint=pprint.PrettyPrinter(indent=4).pprint


def random_chain(**kwargs):
    '''
    Create/Add a single polymer chain composed of same beads
    having length defined by polymerization.
    Required parameters:
        polymerization - length of the polymer chain (>=1)
        bead_radius - bead radius in DPD
            (bond_length equals to 2*bead_radius)
    Optional parameters:
        atom_type=1 - type to assign to all atoms
        bond_type=1 - type to assign to all bonds
        start_atom_id=1 - id of the first atom
        start_bond_id=1 - id of the first bond
        structure - if set, polymer is added into the structure, otherwise
            new structure is created
    '''
    polymerization = kwargs['polymerization']
    bead_radius = kwargs['bead_radius']
    try:
        structure = kwargs['structure']
    except KeyError:
        structure = {'atoms': {}, 'bonds': {}}
    try:
        atom_type = kwargs['atom_type']
    except KeyError:
        try:
            atom_type = max([v['type'] for v in structure['atoms'].values()]) + 1
        except ValueError:
            atom_type = 1
    try:
        bond_type = kwargs['bond_type']
    except KeyError:
        try:
            bond_type = max([v['type'] for v in structure['bonds'].values()]) + 1
        except ValueError:
            bond_type = 1
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
    cell_was_set = False
    try:
        xlo = structure['cell']['xlo']
        xhi = structure['cell']['xhi']
        ylo = structure['cell']['ylo']
        yhi = structure['cell']['yhi']
        zlo = structure['cell']['zlo']
        zhi = structure['cell']['zhi']
        cell_was_set = True
    except KeyError:
        try:
            inflation = 1.25  # if atoms exist but no borders, system increases
            xlo = min(v['x'] for v in structure['atoms'].values()) * inflation
            xhi = max(v['x'] for v in structure['atoms'].values()) * inflation
            ylo = min(v['y'] for v in structure['atoms'].values()) * inflation
            yhi = max(v['y'] for v in structure['atoms'].values()) * inflation
            zlo = min(v['z'] for v in structure['atoms'].values()) * inflation
            zhi = max(v['z'] for v in structure['atoms'].values()) * inflation
        except KeyError:
            chain_length_approximate = 2*bead_radius * polymerization**0.5
            xlo = -chain_length_approximate
            xhi = chain_length_approximate
            ylo = -chain_length_approximate
            yhi = chain_length_approximate
            zlo = -chain_length_approximate
            zhi = chain_length_approximate
    lx = xhi - xlo
    ly = yhi - ylo
    lz = zhi - zlo
    chain = {}
    in_chain_id = 1
    fails_done = 0
    fails_allowed = polymerization * 10

    while fails_done < fails_allowed:
        if not chain:
            x = xlo + random.random() * lx
            y = ylo + random.random() * ly
            z = zlo + random.random() * lz
        else:
            last_atom = chain[max(chain.keys())]
            theta = math.pi * random.random()
            phi = 2*math.pi * random.random()
            x = last_atom['x'] + 2*bead_radius * math.sin(theta) * math.cos(phi)
            y = last_atom['y'] + 2*bead_radius * math.sin(theta) * math.sin(phi)
            z = last_atom['z'] + 2*bead_radius * math.cos(theta)
        all_is_ok = True
        for atom in structure['atoms'].values():
            dx = x - atom['x']
            dy = y - atom['y']
            dz = z - atom['z']
            dr2 = dx**2 + dy**2 + dz**2
            if atom['phase'] == 'filler' and dr2 < 9*bead_radius**2:
                all_is_ok = False
                break
            elif dr2 < 4*bead_radius**2:
                all_is_ok = False
                break
        if not all_is_ok:
            fails_done += 1
            continue
        for atom in chain.values():
            dx = x - atom['x']
            dy = y - atom['y']
            dz = z - atom['z']
            dr2 = dx**2 + dy**2 + dz**2
            if dr2 < 3.5*bead_radius**2:
                all_is_ok = False
                break
        if not all_is_ok:
            continue
        chain[in_chain_id] = {
            'type': atom_type, 'charge': 0,
            'x': x if xlo < x < xhi else x - lx if x > xhi else x + lx,
            'y': y if ylo < y < yhi else y - ly if y > yhi else y + ly,
            'z': z if zlo < z < zhi else z - lz if z > zhi else z + lz,
            'phase': 'matrix'
        }
        in_chain_id += 1
        if len(chain) == polymerization:
            break
    if len(chain) != polymerization:
        return False, {}
    for k, v in chain.items():
        structure['atoms'][atom_id + k - 1] = v
        if k > 1:
            structure['bonds'][bond_id + k - 2] = {
                'type': bond_type,
                'atoms': [atom_id + k - 2, atom_id + k - 1]
            }
    if not cell_was_set:
        xlo = min(v['x'] for v in structure['atoms'].values())
        xhi = max(v['x'] for v in structure['atoms'].values())
        ylo = min(v['y'] for v in structure['atoms'].values())
        yhi = max(v['y'] for v in structure['atoms'].values())
        zlo = min(v['z'] for v in structure['atoms'].values())
        zhi = max(v['z'] for v in structure['atoms'].values())
        structure['cell'] = {
            'xlo': xlo, 'xhi': xhi,
            'ylo': ylo, 'yhi': yhi,
            'zlo': zlo, 'zhi': zhi
        }
    return True, structure
