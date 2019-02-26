import math
import random


def modifier_random(**kwargs):
    '''
    Required parameters:
        head_charge - charge of the head of the modifier
        tail_length - length (in beads) of uncharged tail (>=0)
        bead_raduis - bead radius in DPD
    Optional parameters:
        head_type=1 - type to assign to head atom
        tail_type=1 - type to assign to all tail atoms
        start_atom_id=1 - id of the first atom
        start_bond_id=1 - id of the first bond
        head_tail_bond_type=1 - type of bond between head and tail
        tail_tail_bond_type=2 - type of bond between tail atoms1
        structure - if set, modifier is added into the structure, otherwise
            new structure is created
    '''
    head_charge = kwargs['head_charge']
    tail_length = kwargs['tail_length']
    bead_radius = kwargs['bead_radius']
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
            tail_type = head_type + 1
        except ValueError:
            head_type = 1
            tail_type = 2
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
            rough_molecule_size = (1 + tail_length)**0.5 * 2*bead_radius
            xlo = -rough_molecule_size
            xhi = rough_molecule_size
            ylo = -rough_molecule_size
            yhi = rough_molecule_size
            zlo = -rough_molecule_size
            zhi = rough_molecule_size
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
    lx = xhi - xlo
    ly = yhi - ylo
    lz = zhi - zlo
    in_chain_id = 1
    modifier = {}
    fails_allowed = 100
    fails_done = 0

    while fails_done < fails_allowed:
        if not modifier:
            x = xlo + lx * random.random()
            y = ylo + ly * random.random()
            z = zlo + lz * random.random()
        else:
            x = modifier[max(modifier.keys())]['x']
            y = modifier[max(modifier.keys())]['y']
            z = modifier[max(modifier.keys())]['z']
            theta = math.pi * random.random()
            phi = 2*math.pi * random.random()
            x += 2*bead_radius * math.sin(theta) * math.cos(phi)
            y += 2*bead_radius * math.sin(theta) * math.sin(phi)
            z += 2*bead_radius * math.cos(theta)
        all_is_ok = True
        for atom in structure['atoms'].values():
            dx = x - atom['x']
            dy = y - atom['y']
            dz = z - atom['z']
            dr2 = dx**2 + dy**2 + dz**2
            if dr2 < 3.5*bead_radius**2:
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
            if dr2 < 3.5*bead_radius**2:
                all_is_ok = False
                break
        if not all_is_ok:
            fails_done += 1
            continue
        if not modifier:  # appending head atom
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
    if not cell_was_set:
        xlo = min(v['x'] for v in structure['atoms'].values())
        xhi = max(v['x'] for v in structure['atoms'].values())
        ylo = min(v['y'] for v in structure['atoms'].values())
        yhi = max(v['y'] for v in structure['atoms'].values())
        zlo = min(v['z'] for v in structure['atoms'].values())
        zhi = max(v['z'] for v in structure['atoms'].values())
        structure['cell'] = {
            'xlo': xlo - 2*bead_radius, 'xhi': xhi + 2*bead_radius,
            'ylo': ylo - 2*bead_radius, 'yhi': yhi + 2*bead_radius,
            'zlo': zlo - 2*bead_radius, 'zhi': zhi + 2*bead_radius
        }
    return True, structure


def modifier_near_clays(**kwargs):
    '''
    Required parameters:
        head_charge - charge of the head of the modifier
        tail_length - length (in beads) of uncharged tail (>=0)
        bead_raduis - bead radius in DPD
    Optional parameters:
        head_type=1 - type to assign to head atom
        tail_type=1 - type to assign to all tail atoms
        start_atom_id=1 - id of the first atom
        start_bond_id=1 - id of the first bond
        head_tail_bond_type=1 - type of bond between head and tail
        tail_tail_bond_type=2 - type of bond between tail atoms1
        structure - if set, modifier is added into the structure, otherwise
            new structure is created
    '''
    head_charge = kwargs['head_charge']
    tail_length = kwargs['tail_length']
    bead_radius = kwargs['bead_radius']
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
            tail_type = head_type + 1
        except ValueError:
            head_type = 1
            tail_type = 2
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
            rough_molecule_size = (1 + tail_length)**0.5 * 2*bead_radius
            xlo = -rough_molecule_size
            xhi = rough_molecule_size
            ylo = -rough_molecule_size
            yhi = rough_molecule_size
            zlo = -rough_molecule_size
            zhi = rough_molecule_size
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
    lx = xhi - xlo
    ly = yhi - ylo
    lz = zhi - zlo
    in_chain_id = 1
    modifier = {}
    fails_allowed = 100
    fails_done = 0

    while fails_done < fails_allowed:
        if not modifier:
            x = xlo + lx * random.random()
            y = ylo + ly * random.random()
            z = zlo + lz * random.random()
        else:
            x = modifier[max(modifier.keys())]['x']
            y = modifier[max(modifier.keys())]['y']
            z = modifier[max(modifier.keys())]['z']
            theta = math.pi * random.random()
            phi = 2*math.pi * random.random()
            x += 2*bead_radius * math.sin(theta) * math.cos(phi)
            y += 2*bead_radius * math.sin(theta) * math.sin(phi)
            z += 2*bead_radius * math.cos(theta)
        all_is_ok = True
        near_clay = False
        for atom in structure['atoms'].values():
            dx = x - atom['x']
            dy = y - atom['y']
            dz = z - atom['z']
            dr2 = dx**2 + dy**2 + dz**2
            if dr2 < 9*bead_radius**2 and atom['phase'] == 'filler':
                near_clay = True
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
        if not all_is_ok or not near_clay:
            fails_done += 1
            continue
        if not modifier:  # appending head atom
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
    if not cell_was_set:
        xlo = min(v['x'] for v in structure['atoms'].values())
        xhi = max(v['x'] for v in structure['atoms'].values())
        ylo = min(v['y'] for v in structure['atoms'].values())
        yhi = max(v['y'] for v in structure['atoms'].values())
        zlo = min(v['z'] for v in structure['atoms'].values())
        zhi = max(v['z'] for v in structure['atoms'].values())
        structure['cell'] = {
            'xlo': xlo - 2*bead_radius, 'xhi': xhi + 2*bead_radius,
            'ylo': ylo - 2*bead_radius, 'yhi': yhi + 2*bead_radius,
            'zlo': zlo - 2*bead_radius, 'zhi': zhi + 2*bead_radius
        }
    return True, structure
