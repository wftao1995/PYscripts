#!/usr/bin/env python
#import poscartoolkit
import sys
from ase.io import read
#a = poscartoolkit.cstruc()
#a.readposcar(sys.argv[1])
#a.setmode('D')
#a.write_poscar('POSCAR_DIRECT.vasp')
a = read(sys.argv[1])
if len(sys.argv) ==2 :
    a.write('POSCAR_DIR', vasp5=True, direct = True)
else:
    a.write(sys.argv[2], vasp5=True, direct = True)
