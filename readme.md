# DPD tools

A reposistory to develop simple scripts that are for DPD modeling of a ternary 
nanocomposite MMT - modifier - polymer.


## Physical structures

* MMT: 3-atoms-thick platelet (with bonds between neighboring atoms and angles
equal to pi/2 between (middle_{i-1}, middle_{i}, top_{i}) and 
(middle_{i-1}, middle_{i}, bottom_{i}) that may be charged
* Modifier: 1 charged head and >= 0 not charged tails
* Polymer: >= 1 not charged beads


## Data structures (all are dicts)


### Atom:

{<br>
    'x':      float<br>
    'y':      float<br>
    'z':      float<br>
    'charge': float<br>
    'type':   int<br>
    'phase':  str<br>
}


### Bond / Angle / Dihedral / Improper:

{<br>
    'type':  int<br>
    'atoms': list/set/tuple<br>
}
