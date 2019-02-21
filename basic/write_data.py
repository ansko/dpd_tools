def write_data(out_fname, structure):
    with open(out_fname, 'w') as f:
        f.write('LAMMPS data file by Anton\n\n')
        if 'atoms' in structure.keys() and structure['atoms']:
            f.write('{0} atoms\n'.format(len(structure['atoms'])))
            f.write('{0} atom types\n'.format(len(set([v['type'] for v in
                structure['atoms'].values()]))))
        if 'bonds' in structure.keys() and structure['bonds']:
            f.write('{0} bonds\n'.format(len(structure['bonds'])))
            f.write('{0} bond types\n'.format(len(set([v['type'] for v in
                structure['bonds'].values()]))))
        if 'angles' in structure.keys() and structure['angles']:
            f.write('{0} angles\n'.format(len(structure['angles'])))
            f.write('{0} angle types\n'.format(len(set([v['type'] for v in
                structure['angles'].values()]))))
        if 'dihedrals' in structure.keys() and structure['dihedrals']:
            f.write('{0} dihedrals\n'.format(len(structure['dihedrals'])))
            f.write('{0} dihedral types\n'.format(len(set([v['type'] for v in
                structure['dihedrals'].values()]))))
        if 'impropers' in structure.keys() and structure['impropers']:
            f.write('{0} impropers\n'.format(len(structure['impropers'])))
            f.write('{0} improper types\n'.format(len(set([v['type'] for v in
                structure['impropers'].values()]))))
        f.write('\n')
        for axis in 'xyz':
            lo = structure['cell']['{0}lo'.format(axis)]
            hi = structure['cell']['{0}hi'.format(axis)]
            f.write('{0} {1} {2}lo {2}hi\n'.format(lo, hi, axis))
        f.write('\n')
        if 'atoms' in structure.keys() and structure['atoms']:
            f.write('Atoms\n\n')
            for k, v in structure['atoms'].items():
                q = v['charge'] if 'charge' in v.keys() else 0
                f.write('{0} 1 {1} {2} {3} {4} {5} 0 0 0\n'.format(
                    k, v['type'], q, v['x'], v['y'], v['z']))
        if 'bonds' in structure.keys() and structure['bonds']:
            f.write('\nBonds\n\n')
            for k, v in structure['bonds'].items():
                f.write('{0} {1} {2} {3}\n'.format(k, v['type'], *v['atoms']))
        if 'angles' in structure.keys() and structure['angles']:
            f.write('\nAngles\n\n')
            for k, v in structure['angles'].items():
                f.write('{0} {1} {2} {3} {4}\n'.format(k, v['type'], *v['atoms']))
        if 'dihedrals' in structure.keys() and structure['dihedrals']:
            f.write('\nDihedrals\n\n')
            for k, v in structure['dihedrals'].items():
                f.write('{0} {1} {2} {3} {4} {5}\n'.format(k, v['type'],
                    *v['atoms']))
        if 'impropers' in structure.keys() and structure['impropers']:
            f.write('\nImpropers\n\n')
            for k, v in structure['impropers'].items():
                f.write('{0} {1} {2} {3} {4} {5}\n'.format(k, v['type'],
                    *v['atoms']))
