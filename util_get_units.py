#!/usr/bin/env python3


# This script helps to compute units from real to lj.


def composite_L():
    rho_dpd = 3
    lx = 90; ly = 82; lz = 64;  # box dimensions
    th_mmt = 8  # mmt thickness (excluded volume for soft phase)
    soft_atoms_count = 70*108 + 382*90
    mmt_atoms_count = 6480
    soft_beads_count = 70*3 + 20*90
    soft_bead_volume = lx * ly * (lz - th_mmt) / soft_beads_count
    dpd_cutoff = (soft_bead_volume * rho_dpd)**(1/3)
    print('dpd cutoff:', dpd_cutoff)
    soft_bead_radius = (soft_bead_volume * 3/4/3.14)**(1/3)
    print('soft bead radius:', soft_bead_radius)


def main():
    composite_L()


if __name__ == '__main__':
    main()
