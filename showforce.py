#!/usr/bin/env python
from ase.io import read
import numpy as np
import sys

if len(sys.argv) == 2:
    a0 = read('vasprun.xml', index=':')
    for i in a0:
        forces = i.get_forces()
        lst = []
        for j in forces:
            lst.append(np.linalg.norm(j))
        print(max(lst))
    exit()


a0 = read('OUTCAR')
# print(a0.get_forces())
forces = a0.get_forces()
lst = []
for i in forces:
	lst.append(np.linalg.norm(i))

print(lst)
print('maxforce: ', max(lst))
