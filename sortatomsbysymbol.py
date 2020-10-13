#!/usr/bin/env python

from ase.build import sort
from ase.io import read
import sys

def sortatomsbysym(atom, syms):
	tags = atom.get_chemical_symbols()
	#print([(tag, i) for i, tag in enumerate(tags)])
	decs = sorted([(tag, i) for i, tag in enumerate(tags)], key = lambda x: syms.index(x[0]))
	indices = [i for tag, i in decs]
	return atom[indices]

atom = read(sys.argv[1])
chemicals = list(set(atom.get_chemical_symbols()))
print("symbols in POSCAR:" + " " + str(chemicals))
nsyms = input("input chemical symbols, seprate using whitespace:").split()
atom = sortatomsbysym(atom, nsyms)
atom.write('POSCAR_SORTED', format = "vasp", vasp5 = True, direct = True)


