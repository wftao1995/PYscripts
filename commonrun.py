import os
from ase.io.trajectory import Trajectory

MAGMOM = {
    'Fe' : 3,
    'Mn' : 4,
    'La' : 0.8,
    'H'  : 0.6,
    'Y'  : 0.8,
    'Ce' : 1.2,
    'Gd' : 8,
    'C'  : 0.6,
    'Cr' : 0.6,
    'V'  : 0.7,
    'Mo' : 0.5,
    'Sc' : 3,
    'Pr' : 5,
    'Nd' : 5,
    'Pm' : 5,
    'Sm' : 5,
    'Eu' : 5,
    'Tb' : 5,
    'Dy' : 5,
    'Ho' : 5,
    'Er' : 5,
    'Tm' : 5,
    'Yb' : 5,
    'Lu' : 5,
    'Ni' : 4,
    'Cu' : 1,
    'Co' : 4,
    'W' : 0.6
}

def set_ini_mag(atom):
    inimagmom = []
    for i in atom:
        inimagmom += [MAGMOM[i.symbol]]
    atom.set_initial_magnetic_moments(inimagmom)

# args: calc_1, calc_2, ...
def relaxstr(atom, steps, *args):
    rootdir = os.getcwd()
    traj = Trajectory('realx.traj','w')
    if not os.path.exists('relax'):
        os.mkdir('relax')
    os.chdir('relax')
    for i in range(steps):
        try:
            atom.set_calculator(args[i])
            atom.get_potential_energy()
            traj.write(atom)
        except:
            print('Error occurs at step' + str(i+1))
            fp = open('_TAG_FAILED','w')
            fp.close()
            os.chdir(rootdir)
            return 0
    os.chdir(rootdir)
    return 1
 
def getenergy(atom, calc, dirname = "static"):
    rootdir = os.getcwd()
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    os.chdir(dirname)
    try:
        atom.set_calculator(calc)
        atom.get_potential_energy()
    except:
        print("Error while get energy")
        fp = open('_TAG_FAILED','w')
        fp.close()
        os.chdir(rootdir)
        return 0
    os.chdir(rootdir)
    return 1

def enccoverge(atoms, calc, enclist, dirname = 'enccoverge'):
    import matplotlib.pyplot as plt
    rootdir = os.getcwd()
    traj = Trajectory('encoverge.traj','w')
    engs = []
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    os.chdir(dirname)
    for i in enclist:
        calc.set(encut = i)
        try:
            atom = atoms.copy()
            atom.set_calculator(calc)
            eng = atom.get_potential_energy()
            engs.append(eng)
            traj.write(atom)
        except:
            print("Error while encut = ", str(i))
    plt.plot(enclist, engs, '-')
    plt.xlabel('ENCUT/eV')
    plt.ylabel('Energe/eV')
    plt.savefig('encoverge.svg')
    os.chdir(rootdir)

def kptconverge(atom, calc, kptlist, dirname = 'kptcoverge'):
    import matplotlib.pyplot as plt
    rootdir = os.getcwd()
    traj = Trajectory('kptcoverge.traj','w')
    engs = []
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    os.chdir(dirname)
    for i in kptlist:
        calc.set(kpts = i)
        try:
            atom.set_calculator(calc)
            eng = atom.get_potential_energy()
            engs.append(eng)
            traj.write(atom)
        except:
            print("Error while kpts = ", str(i))
    plt.plot(engs, '-')
    plt.xlabel('KPT')
    plt.ylabel('Energe/eV')
    plt.savefig('kptcoverge.svg')
    os.chdir(rootdir)

# see 'Equation of States' section of ase documentation, many methods provided here
# supported eos: 
# sj, taylor, murnaghan, birch, birchmurnaghan, pouriertarantola, vient, antonschmidt and p3.
# for sj, E(V) = c_0 + c_1*t + c_2*t^2 + c_3*t^3, while t = cubicroot(V)^-1
# for p3, E(V) = c_0 + c_1*V + c_2*V^2 + c_3*V^3
# for parabola, 
def eosfit(atoms, calc, begin = 0.95, end = 1.05, steps = 5, method = "birchmurnaghan"):
    from ase.eos import EquationOfState
    from numpy import linspace
    rootdir = os.getcwd()
    if not os.path.exists('eosfit'):
        os.mkdir('eosfit')
    os.chdir('eosfit')
    cell = atoms.get_cell()
    atoms.set_calculator(calc)
    traj = Trajectory('eosfit.traj','w')
    for x in linspace(begin, end, steps):
        print(str(x))
        atoms.set_cell(cell*x, scale_atoms=True)
        atoms.get_potential_energy()
        traj.write(atoms)
    configs = Trajectory('eosfit.traj','r')
    volumes = [at.get_volume() for at in configs]
    energies = [at.get_potential_energy() for at in configs]
    eos = EquationOfState(volumes, energies, eos = method)
    v0, e0, B = eos.fit()
    eos.plot('eosfit.svg')
    fp = open('eosdata.dat')
    fp.write('Volume\t\tEnergy', 'w')
    for i in range(len(configs)):
        fp.write(str(volumes[i]) + '\t' + str(energies[i]) + '\n')
    fp.write('V0, E0, B: ' + str(v0) + '\t' + str(e0) + '\t' + str(B) + '\n')
    fp.close()
    os.chdir(rootdir)

# specify a range, list or array
def eosfit_spec(atoms, calc, rg, method = "birchmurnaghan"):
    from ase.eos import EquationOfState
    from numpy import linspace
    rootdir = os.getcwd()
    if not os.path.exists('eosfit'):
        os.mkdir('eosfit')
    os.chdir('eosfit')
    cell = atoms.get_cell()
    atoms.set_calculator(calc)
    traj = Trajectory('eosfit.traj','w')
    for x in rg:
        print(str(x))
        atoms.set_cell(cell*x, scale_atoms=True)
        atoms.get_potential_energy()
        traj.write(atoms)
    configs = Trajectory('eosfit.traj','r')
    volumes = [at.get_volume() for at in configs]
    energies = [at.get_potential_energy() for at in configs]
    eos = EquationOfState(volumes, energies, eos = method)
    v0, e0, B = eos.fit()
    eos.plot('eosfit.svg')
    fp = open('eosdata.dat', 'w')
    fp.write('Volume\t\tEnergy')
    for i in range(len(configs)):
        fp.write(str(volumes[i]) + '\t' + str(energies[i]) + '\n')
    os.chdir(rootdir)

