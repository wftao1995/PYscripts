#!/usr/bin/env python
#import poscartoolkit
import sys
from ase.io import read
#a = poscartoolkit.cstruc()
#a.readposcar(sys.argv[1])
#a.setmode('C')
#a.write_poscar('POSCAR_CARSCOORD.vasp')
a = read('sys.argv[1]')
a.write('POSCAR_CAR', vasp5 = True, direct = True)
