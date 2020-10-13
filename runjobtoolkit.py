#!/usr/bin/env python
import os
import sys
from numpy import pi
import ase.io.vasp
from ase.io.vasp import read_vasp,read_vasp_xml
import ase.dft.kpoints
import ase.calculators.vasp
from ase import Atoms
from ase.data import chemical_symbols, atomic_numbers
import ase.build
from ase.build import cut,bulk
import numpy as np
import ase.geometry
from ase.geometry import get_layers
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d,interp2d,interpnd
import matplotlib as mpl
# one-, two-, or n-level run is provided
# not suitalbe for calculations like band structure
# sequenial, an ionic or lattice relaxation will affect all steps afterwards
# if seq = True, *args should be arranged as follows: calc_n, dirname_n ..., I won't check these in the code
# so check the args carefully before submit to the job server
def runcalc(atom, calc, dirname, do_static = True, calc_s = 0, seq = False, *args):
    rootdir = os.getcwd()
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    os.chdir(dirname)
    atom.set_calculator(calc)
    eng = 0
    try: 
        eng = atom.get_potential_energy()
    except:
        fp = open('_TAG_FAILED','w')
        fp.close()
        return 0
    if do_static:
        if not os.path.exists('static'):
            os.mkdir('static')
        os.chdir('static')
        try :
            atom.set_calculator(calc_s)
            eng = atom.get_potential_energy()
        except:
            fp = open('_TAG_STATIC_FAILED','w')
            fp.close()
            return 0
    os.chdir(rootdir)
    if seq:
        #True: dir, False: calc
        t = True
        rootdir = os.getcwd()
        for i in args:
            if t:
                t = not t
                os.chdir(rootdir)
                atom.set_calculator(i)
            else:
                t = not t
                if not os.path.exists(i):
                    os.mkdir(i)
                os.chdir(i)
                try :
                    atom.get_potential_energy()
                except :
                    fp = open('_TAG_STATIC_FAILED','w')
                    fp.close()
    return eng
#for test purpose only
def runcalc_dymmy(atom, calc, dirname, do_static = True, calc_s = 0):
    rootdir = os.getcwd()
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    os.chdir(dirname)
    atom.set_calculator(calc)
    eng = 0
    try: 
        if not os.path.exists('DUMMY_RUN'):
            fp = open('DUMMY_RUN','w')
            fp.close()
    except:
        fp = open('_TAG_FAILED','w')
        fp.close()
        return 0
    if do_static:
        if not os.path.exists('static'):
            os.mkdir('static')
        os.chdir('static')
        try :
            atom.set_calculator(calc_s)
            if not os.path.exists('DUMMY_RUN'):
                fp = open('DUMMY_RUN','w')
                fp.close()
        except:
            fp = open('_TAG_STATIC_FAILED','w')
            fp.close()
            return 0
    os.chdir(rootdir)
    return eng
#wrap of ase.geometry.get_layers()
def set_layer_info(atoms, miller = (0,0,1), tolerance = 0.001):
    tags = get_layers(atoms, miller, tolerance)
    atoms.set_tags(tags[0])
#wrap of ase.build.sort
#note that the ase.build.sort method retruns an ordered atom object instead of doing an on-place run
def sort_by_layer_info(atoms):
    atoms = ase.build.sort(atoms, atoms.get_tags())
#Here a classical Universal Binding Energy Curves (UBEC) method is used (only dummy now, @TODO: finish this).
#interplotation is done via cubic spineline
#atoms : Atoms object; n1: layers of phase 1, dists: a array of distances; axis: 0-x,1-y,2-z
#no realxation here, @TODO: relaxation
#note in ase, all coordinates are sorted in carstain mode
# for some surfaces, the optimial distance might be negative
# one may condiser a distance range contains negative values in such circumstances
def get_opt_interface_seperation_distance(atoms, n1, dists, calc, axis = 2):
    if sum(atoms.get_tags()) == 0:
        raise Exception("Please provide an atom object with tag info!") 
    st = ase.build.sort(atoms, atoms.get_tags())
    natoms = len(st.positions)
    layers = st.get_tags()
    revlayers = layers.tolist().reverse()
    posid1 = revlayers.index(n1)
    posid2 = posid1 + 1
    offset = st.positions[posid2][axis] - st.positions[posid1][axis]
    cnt = 0
    rootdir = os.getcwd()
    engs = []
    for i in dists:
        os.chdir(rootdir)
        cst = st.copy()
        if not os.path.exists(str(cnt)):
            os.mkdir(str(cnt))  
        cnt = cnt + 1  
        os.chdir(str(cnt))
        for j in range(posid2, natoms):
            cst.positions[j][axis] = cst.positions[j][axis] - offset + i
        cst.set_calculator(calc)
        try:
            eng = cst.get_potential_energy()
            engs.append(eng)
        except:
            fp = open('_TAG_FAILID','w')
            fp.close()
            print('---ERROR---')
            raise Exception('CALCULATION FAILED, SEE vasp.out FOR DETIALS')
    #post-procesing
    os.chdir(rootdir)
    f = interp1d(dists,engs,kind='cubic')
    xnew = np.linspace(dists[0], dists[-1], num=np.floor((dists[-1]-dists[0])/0.01), endpoint=True)
    plt.plot(dists,engs,'o',xnew,f(xnew),'-')
    plt.savefig('energies.svg')
    fp = open('energies.dat', 'w')
    fp.write('dist: '.join(f'{j:16f}' for j in dists) + '\n')
    fp.write('energy: '.join(f'{j:16f}' for j in engs) + '\n')
    fp.close()
    fp = open('plotdata.dat','w')
    fp.write('xnew: '.join(f'{j:16f}' for j in xnew) + '\n')
    fp.write('ynew: '.join(f'{j:16f}' for j in f(xnew)) + '\n')
    fp.close()

def shift_by_layers_upper(atoms, interpos, a1 = 0, a2 = 0, a3 = 0, wp = True):
    natoms = len(atoms.positions)
    dt = [a1, a2, a3]
    for i in range(interpos+1, natoms):
        atoms.positions[i] = atoms.positions[i] + dt
    if wp:
        atoms.wrap()
# note that for a 10x10 grid, 100 calculations are needed, one may want to use a minimal POTCAR, and a loose k-mesh
# and check symmtery to reduce the mesh size
# do not relax the structure or relax in $axis direction only
def get_opt_layer_mishift(atoms, n1, maxx, maxy, interx, intery, calc, axis = 2):
    if sum(atoms.get_tags()) == 0:
        raise Exception("Please provide an atom object with tag info!") 
    st = ase.build.sort(atoms, atoms.get_tags())
    layers = st.get_tags()
    stepx = maxx/interx
    stepy = maxy/intery
    revlayers = layers.tolist().reverse()
    posid1 = revlayers.index(n1)
    engs = np.zeros((interx,intery),dtype=float)
    rootdir = os.getcwd()
    meshx, meshy = np.mgrid[0.0:maxx:interx*1j, 0.0:maxy:intery*1j]
    for i in range(interx):
        for j in range(intery):
            atc = st.copy()
            atc.set_calculator(calc)
            shift_by_layers_upper(atc, posid1, i*stepx, j*stepy)
            try:
                if os.path.exists(str(i) + '-' + str(j)):
                    os.mkdir(str(i) + '-' + str(j))
                os.chdir(str(i) + '-' + str(j))
                engs[i][j] = atc.get_potential_energy()
            except:
                fp = open('_TAG_FAIL','w')
                fp.close()
            os.chdir(rootdir)
    fnc = interp2d(meshx,meshy,engs,kind='cubic')
    nx = np.linspace(0,maxx,maxx/0.02,endpoint=True)
    ny = np.linspace(0,maxy,maxy/0.02,endpoint=True)
    nz = fnc(nx,ny)
    fig, ax = plt.subplots(1,1)
    #jet dark blue -> green -> red
    cp = ax.contourf(nx,ny,nz,80,cmap = mpl.cm.get_cmap(name = 'jet'))
    fig.colorbar(cp)
    plt.savefig('energies.svg')
    fp = open('energies.dat','w')
    fp.write('Paras:\n')
    fp.write(str(maxx) + ' ' + str(maxy) + ' ' + str(interx) + ' ' + str(intery) + '\n')
    # fp.write('Mesh x:\n')
    # for i in meshx:
    #     fp.write(' '.join(f'{j:16f}' for j in meshx[i]) + '\n')
    # fp.write('Mesh y:\n')
    # for i in meshy:
    #     fp.write(' '.join(f'{j:16f}' for j in meshy[i]) + '\n') 
    fp.write('eng:\n')
    for i in range(interx):
        fp.write(' '.join(f'{j:16f}' for j in engs[i]) + '\n')
    fp.close()
    fp = open('plot.dat')
    fp.write('x:\n')
    fp.write(' '.join(f'{j:16f}' for j in nx) + '\n')
    fp.write('y:\n')
    fp.write(' '.join(f'{j:16f}' for j in ny) + '\n')
    fp.write('z:\n')
    for i in nz:
        fp.write(' '.join(f'{j:16f}' for j in i) + '\n')
    fp.close()
    

#WARNNING: COMPUTATION VERY DEMANDING, AND THE RESULTS ARE USUALLY USELESS
#not tested yet
#since scipy only support linear and nearest interpolation, the density countour-surface might be not so smooth
def get_opt_layer_real(atoms, n1, dists, maxx, maxy, interx, intery, calc, axis = 2):
    interz = len(dists)
    natoms = len(atoms.positions)
    if sum(atoms.get_tags()) == 0:
        raise Exception("Please provide an atom object with tag info!") 
    st = ase.build.sort(atoms, atoms.get_tags())
    layers = st.get_tags()
    stepx = maxx/interx
    stepy = maxy/intery
    revlayers = layers.tolist().reverse()
    posid1 = revlayers.index(n1)
    posid2 = posid1 + 1
    offset = st.positions[posid2][axis] - st.positions[posid1][axis]
    engs = np.zeros((interx,intery),dtype=float)
    rootdir = os.getcwd()
    meshx, meshy = np.mgrid[0.0:maxx:interx*1j, 0.0:maxy:intery*1j]
    cnt = 0
    for z in range(interz):
        os.chdir(rootdir)
        atc = st.copy()
        if not os.path.exists(str(cnt)):
            os.mkdir(str(cnt))    
        os.chdir(str(cnt))
        cnt = cnt + 1
        for j in range(posid2, natoms):
            atc.positions[j][axis] = atc.positions[j][axis] - offset + i
        crootdir = os.getcwd()
        for i in range(interx):
            for j in range(intery):
                #atc = st.copy()
                atc.set_calculator(calc)
                shift_by_layers_upper(atc, posid1, i*stepx, j*stepy)
                try:
                    if os.path.exists(str(i) + '-' + str(j)):
                        os.mkdir(str(i) + '-' + str(j))
                    os.chdir(str(i) + '-' + str(j))
                    engs[z][i][j] = atc.get_potential_energy()
                except:
                    fp = open('_TAG_FAIL','w')
                    fp.close()
                os.chdir(crootdir)
    #@TODO: finish logs


#arragments: calc_n, kpt_n,
def relaxstr(atom, steps, *args):
    if not os.path.exist('relax'):
        os.mkdir('relax')
    os.chdir('relax')
