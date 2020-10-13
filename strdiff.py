#!/usr/bin/env python
from ase.io import read
import sys
import numpy as np

def wrapcenter(atom):
    posarr = atom.get_scaled_positions(wrap = False)
    for i in range(len(posarr)):
        for j in range(3):
            if posarr[i][j] > 0.5:
                posarr[i][j] = posarr[i][j] - 1
            elif posarr[i][j] < -0.5:
                posarr[i][j] = posarr[i][j] + 1
    atom.set_scaled_positions(posarr)

if len(sys.argv) == 1:
    pos1 = read('POSCAR')
    pos2 = read('CONTCAR')
elif sys.argv[1] == 'XDATCAR':
    a = read('XDATCAR',index=':')
    pos1 = a[-2]
    pos2 = a[-1]
else:
    pos1 = read(sys.argv[1])
    pos2 = read(sys.argv[2])

posarr1 = pos1.get_positions()
posarr2 = pos2.get_positions()

cell1 = pos1.get_cell()
cell2 = pos2.get_cell()
distc1 = np.linalg.norm(cell1[0] - cell2[0])
distc2 = np.linalg.norm(cell1[1] - cell2[1])
distc3 = np.linalg.norm(cell1[2] - cell2[2])
distmax = max(distc1, distc2, distc3)
if distmax > 10e-8:
    print('Warning: base vectors are not the same\n')
    vol1 = np.linalg.det(cell1)
    vol2 = np.linalg.det(cell2)
    print('vol1, vol2: ' + str(vol1) + '    ' + str(vol2))

posdiff = posarr1 - posarr2
pos1.set_positions(posdiff)
wrapcenter(pos1)
magtit = 0
totnormdiff = 0
for i in pos1.positions:
    norm = np.linalg.norm(i)
    totnormdiff = totnormdiff + norm
    magtit = magtit + norm**2
magtit = np.sqrt(magtit)
print('Total: ' + str(totnormdiff))
print('Average: ' + str(totnormdiff/len(posdiff)))
print('TOT: ' + str(magtit))
print('AVG: ' + str(magtit/len(posdiff)))



pos1.write('POSCAR.diff', format="vasp", vasp5 = True, direct = False)
# pos1.write('POSCAR.direct.diff', format="vasp", vasp5 = True, direct = True)

