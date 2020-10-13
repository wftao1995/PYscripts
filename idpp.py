#!/usr/bin/env python
#accepted arguments: ini fin images, accepted formats: POSCAR CONTCAR ase.traj
import os
import sys 
from ase.neb import NEB
from ase.io import read
from ase.io.trajectory import Trajectory

def i2str(num):
    if num < 10:
        return '0' + str(num)
    else : 
        return str(num)
#accept POSCAR OUTCAR or xml
def do_vasp(pos1, pos2, snimgs):
    nimgs = int(snimgs)
    ini = read(pos1)
    fin = read(pos2)
    ini.wrap()
    fin.wrap()
    traj = [ini]
    traj += [ini.copy() for i in range(nimgs)]
    traj += [fin]
    neb = NEB(traj)
    neb.interpolate('idpp')
    images = neb.images
    if not os.path.exists('00'):
        os.mkdir('00')
    if not os.path.exists(i2str(nimgs+1)):
        os.mkdir(i2str(nimgs+1))
    images[0].write('00/POSCAR', vasp5 = True, direct = True)
    images[-1].write(i2str(nimgs+1)+'/POSCAR', vasp5 = True, direct = True)
    for i in range(nimgs):
        if not os.path.exists(i2str(i+1)):
            os.mkdir(i2str(i+1))
        images[i+1].write(i2str(i+1) + '/POSCAR', vasp5 = True, direct = True)
    traj = Trajectory('images.tarj','w')
    for i in images:
        traj.write(i)

do_vasp(sys.argv[1], sys.argv[2], sys.argv[3])
