import numpy as np
import re

class CFG:
    def __init__(self):
        self.pos = []
        self.size = 0
        self.force = []
        self.types = []
        self.energy = 0
        self.cell = np.zeros((3, 3))
        self.stress = np.zeros(6)
        self.fenergy = 0
        self.have_energy = False
        self.have_stress = False
        self.have_forces = False
        self.have_features = False
        self.features = {}
    # Returns average magnitude of forces
    def avg_force(self):
        total = 0
        for i in range(self.size):
            total += np.linalg.norm(self.force[i])
        return total / self.size
    # Returns magnitude of forces
    def force_list(self, negflag=False):
        flist = []
        if negflag:
            istart = -1
        else:
            istart = 1
        sgn = 1
        for i in range(self.size):
            sgn *= istart
            flist.append(sgn*np.linalg.norm(self.force[i]))
        return flist

    def setfenergy(self, edict):
        self.fenergy = self.energy
        for i in self.types:
            self.fenergy -= edict[i]

    def read_poscar(self, filename):
        fp = open(filename, 'r')
        cfp = fp.readlines()
        fp.close()
        self.have_energy = False
        self.have_forces = False
        self.have_stress = False
        self.have_features = False
        scale = float(cfp[1])
        for i in range(3):
            line = cfp[2+i].strip()
            vt = np.fromstring(line, sep=' ')
            self.cell[i] = vt * scale
        ncomp = np.fromstring(cfp[6].strip(), sep=' ', dtype=int)
        for i in range(len(ncomp)):
            self.types += [i for j in range(ncomp[i])]
        self.size = sum(ncomp)
        line = cfp[7].strip()
        if line[0] == 'S' or line[0] == 's':
            ntp = 8
        else:
            ntp = 7
        line = cfp[ntp].strip()
        if line[0] == 'D' or line[0] == 'd':
            direct = True
        for i in range(self.size):
            line = cfp[ntp + i + 1].strip()
            vt = np.fromstring(line, sep=' ')
            if direct:
                # strip
                for i in range(3):
                    if vt[i] < 0:
                        vt[i] += 1
                    if vt[i] >= 1:
                        vt[i] -= 1
                self.pos.append(np.dot(vt, self.cell))
            else:
                self.pos.append(vt)

    def wrap(self):
        invcell = np.linalg.inv(self.cell)
        for i in range(len(self.pos)):
            dirpos = np.dot(self.pos[i], invcell)
            for j in range(3):
                if dirpos[j] < 0:
                    dirpos[j] += 1
                if dirpos[j] >= 1:
                    dirpos[j] -= 1
            self.pos[i] = np.dot(dirpos, self.cell)

    def unwrap(self, idealpos, cut=0.5):
        spos = self.get_scaled_positions()
        for i in range(self.size):
            for j in range(3):
                if abs(idealpos[i][j]) < 10e-3:
                    if spos[i][j] > -cut and spos[i][j] < cut:
                        pass
                    elif spos[i][j] > cut:
                        spos[i][j] -= 1
        for i in range(self.size):
            self.pos[i] = np.dot(spos[i], self.cell)

    def get_scaled_positions(self):
        spos = self.pos.copy()
        invcell = np.linalg.inv(self.cell)
        for i in range(len(spos)):
            spos[i] = np.dot(spos[i], invcell)
        return spos

    def get_toteng_ha(self):
        return 0.0367502*self.energy

    def get_force_ha_bohr(self):
        flist = self.force.copy()
        for i in range(len(flist)):
            flist[i] = flist[i] * 0.0367502 / 1.8897259886
        return flist

    def get_toteng(self, unit):
        # Unit: ha, eV, kJ/mol, kcal/mol, cm-1 K J Hz
        if unit == 'ha':
            return 0.0367502*self.energy
        if unit == 'eV':
            return self.energy
        if unit == 'kJ/mol':
            return 96.4869*self.energy
        if unit == 'kcal/mol':
            return 23.0609*self.energy
        if unit == 'cm-1':
            return 8065.73*self.energy
        if unit == 'K':
            return 11604.9*self.energy
        if unit == 'J':
            return 1.60210e-19*self.energy
        if unit == 'Hz':
            return 2.41804e14*self.energy

    def write(self, fp, fenergy=False):
        fp.write("BEGIN_CFG\n")
        fp.write(" Size\n")
        fp.write("    ")
        fp.write(str(self.size) + '\n')
        fp.write(" Supercell\n")
        for i in range(3):
            fp.write("        " +
                     " ".join(f'{j:.10f}' for j in self.cell[i]) + '\n')
        if self.have_forces:
            fp.write(
                " AtomData:  id type       cartes_x      cartes_y      cartes_z           fx          fy          fz\n")
            for i in range(self.size):
                fp.write("            " + str(i+1) +
                         "    " + str(self.types[i]))
                fp.write("   " + "   ".join(f'{j:.10f}' for j in self.pos[i]))
                fp.write(
                    "   " + "   ".join(f'{j:.10f}' for j in self.force[i]))
                fp.write('\n')
        else:
            fp.write(
                " AtomData:  id type       cartes_x      cartes_y      cartes_z\n")
            for i in range(self.size):
                fp.write("            " + str(i+1) +
                         "    " + str(self.types[i]))
                fp.write("   " + "   ".join(f'{j:.10f}' for j in self.pos[i]))
                fp.write('\n')
        if self.have_energy:
            fp.write(" Energy\n")
            if fenergy:
                fp.write("        " + str(self.fenergy) + "\n")
            else:
                fp.write("        " + str(self.energy) + "\n")
        if self.have_stress:
            fp.write(
                " PlusStress:  xx          yy          zz          yz          xz          xy\n")
            fp.write(
                "    " + " ".join(f'{j:.10f}' for j in self.stress) + '\n')
        if self.have_features:
            fkeys = self.features.keys()
            for i in fkeys:
                fp.write(" Feature    " + i + "    " + self.features[i] + '\n')
        fp.write("END_CFG\n\n")


class CFGS:
    def __init__(self):
        self.cfgs = []
        self.ncfgs = 0

    def readcfg(self, filename):
        fp = open(filename)
        contents = fp.readlines()
        i = 0
        while i < len(contents):
            line = contents[i].strip()
            if re.match("BEGIN", line):
                self.ncfgs += 1
                cfg = CFG()
                i += 1
                continue
            if re.match("Size", line):
                i += 1
                cfgsize = int(contents[i].strip())
                cfg.size = cfgsize
                i += 1
                continue
            if re.match("Supercell", line):
                i += 1
                cell = np.zeros((3, 3))
                for j in range(3):
                    vt = contents[i+j]
                    cell[j] = np.fromstring(vt, dtype=float, sep=' ')
                cfg.cell = cell.copy()
                i += 3
                continue
            if re.match("AtomData", line):
                cfg.have_forces = True
                types = []
                pos = []
                force = []
                i += 1
                for j in range(cfgsize):
                    vt = contents[i+j]
                    vtf = np.fromstring(vt, sep=' ')
                    types.append(int(vtf[1]))
                    pos.append(vtf[2:5])
                    force.append(vtf[5:8])
                cfg.pos = pos.copy()
                cfg.types = types.copy()
                cfg.force = force.copy()
                i += cfgsize
                continue
            if re.match("Energy", line):
                cfg.have_energy = True
                i += 1
                cfg.energy = float(contents[i].strip())
                i += 1
                continue
            if re.match("PlusStress", line):
                cfg.have_stress = True
                i += 1
                vt = contents[i]
                vtf = np.fromstring(vt, sep=' ')
                cfg.stress = vtf.copy()
                i += 1
                continue
            if re.match("Feature", line):
                cfg.have_features = True
                features = line.split()
                fkey = features[1]
                fvaule = features[2]
                cfg.features[fkey] = fvaule
                i += 1
                continue
            if re.match("END", line):
                self.cfgs.append(cfg)
                i += 1
                continue
            else:
                i += 1
                continue

    def strip_unconverge(self):
        cfgs = []
        for i in self.cfgs:
            if not i.have_features:
                raise("Error: no features!")
            if i.features["EFS_by"] == "VASP_not_converged":
                pass
            else:
                cfgs.append(i)
        self.cfgs = cfgs

    def force_list(self, negflag=False):
        flist = []
        for i in range(self.ncfgs):
            flist += self.cfgs[i].force_list(negflag)
        return flist

    def energy_list(self):
        elist = []
        for i in range(self.ncfgs):
            elist.append(self.cfgs[i].energy)
        return elist

    def energy_pa_list(self):
        elist = []
        for i in range(self.ncfgs):
            elist.append(self.cfgs[i].energy / self.cfgs[i].size)
        return elist

    def fenergy(self, edict):
        for i in self.cfgs:
            i.setfenergy(edict)

    def write_cfgs(self, filename, fenergy=False):
        fp = open(filename, 'w')
        for i in self.cfgs:
            i.write(fp, fenergy)

    def append(self, cfg):
        self.ncfgs += 1
        self.cfgs.append(cfg)

    def wrap_all(self):
        for i in self.cfgs:
            i.wrap()

    def unwrap_allequal(self, idealpos, cut=0.5):
        for i in self.cfgs:
            i.unwrap(idealpos, cut)

def compare(cfgref, cfgefs, format="svg"):
    import matplotlib.pyplot as plt
    fref = cfgref.force_list(True)
    fefs = cfgefs.force_list(True)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(fref, fefs, s=1)
    xmax = 1.02*max(max(fref), max(fefs))
    xmin = 0.98*min(min(fref), min(fefs))
    x = np.linspace(xmin, xmax, 10)
    ax.plot(x, x, '-', color='red', linewidth=0.8)
    plt.xlim(xmin, xmax)
    plt.ylim(xmin, xmax)
    plt.xlabel("Reference data")
    plt.ylabel("MTP predict data")
    plt.savefig('figf.' + format, dpi=300)

    eref = cfgref.energy_pa_list()
    eefs = cfgefs.energy_pa_list()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(eref, eefs, s=1)
    xmax = 1.0001*max(max(eref), max(eefs))
    xmin = 0.9999*min(min(eref), min(eefs))
    x = np.linspace(xmin, xmax, 10)
    ax.plot(x, x, '-', color='red', linewidth=0.8)
    plt.xlim(xmin, xmax)
    plt.ylim(xmin, xmax)
    plt.xlabel("Reference data")
    plt.ylabel("MTP predict data")
    plt.savefig('fige.' + format, dpi=300)

def comparetf(cfgreft, cfgreff, cfgefst, cfgefsf, format="svg"):
    import matplotlib.pyplot as plt
    fref = cfgreft.force_list(True)
    fefs = cfgefst.force_list(True)
    freff = cfgreff.force_list(True)
    fefsf = cfgefsf.force_list(True)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(fref, fefs, marker='^', label='Train', s=1.4)
    ax.scatter(freff, fefsf, marker='o', label='Test', s=1.4)
    ax.legend()
    xmax = 1.02*max(max(fref), max(fefs))
    xmin = 0.98*min(min(fref), min(fefs))
    x = np.linspace(xmin, xmax, 10)
    ax.plot(x, x, '-', color='red', linewidth=0.8)
    plt.xlim(xmin, xmax)
    plt.ylim(xmin, xmax)
    plt.xlabel("Reference data")
    plt.ylabel("MTP predict data")
    plt.savefig('figf.' + format, dpi=300)

    eref = cfgreft.energy_pa_list()
    eefs = cfgefst.energy_pa_list()
    ereff = cfgreff.energy_pa_list()
    eefsf = cfgefsf.energy_pa_list()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(eref, eefs, marker='^', label='Train', s=5)
    ax.scatter(ereff, eefsf, marker='o', label='Test', s=5)
    ax.legend()
    xmax = 1.0001*max(max(eref), max(eefs))
    xmin = 0.9999*min(min(eref), min(eefs))
    x = np.linspace(xmin, xmax, 10)
    ax.plot(x, x, '-', color='red', linewidth=0.8)
    plt.xlim(xmin, xmax)
    plt.ylim(xmin, xmax)
    plt.xlabel("Reference data")
    plt.ylabel("MTP predict data")
    plt.savefig('fige.' + format, dpi=300)

def compare_gendata(cfgref, cfgefs):
    fp1 = open('plotdata-f.dat', 'w')
    fp2 = open('plotdata-e.dat', 'w')
    fref = cfgref.force_list(True)
    fefs = cfgefs.force_list(True)
    fp1.write("DFT MLP\n")
    for i in range(len(fref)):
        fp1.write(f'{fref[i]:.10f}')
        fp1.write("  ")
        fp1.write(f'{fefs[i]:.10f}')
        fp1.write('\n')

    eref = cfgref.energy_pa_list()
    eefs = cfgefs.energy_pa_list()
    fp2.write("DFT MLP\n")
    for i in range(len(eref)):
        fp2.write(f'{eref[i]:.10f}')
        fp2.write("  ")
        fp2.write(f'{eefs[i]:.10f}')
        fp2.write('\n')

def force_diff(cfgref, cfgefs):
    fdiff = []
    for i in range(cfgref.ncfgs):
        for j in range(cfgref.cfgs[i].size):
            fdiff += [np.linalg.norm(cfgref.cfgs[i].force[j] - cfgefs.cfgs[i].force[j])]
    return fdiff

def compareforce(cfgref, cfgefs, format="svg"):
    import matplotlib.pyplot as plt
    fdiff = force_diff(cfgref, cfgefs)
    fref = cfgref.force_list(False)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(fref, fdiff, s=1)
    xmax = 1.02*max(max(fref), max(fref))
    xmin = 0.98*min(min(fref), min(fref))
    ymax = 1.1*max(fdiff)
    plt.xlim(xmin, xmax)
    plt.ylim(0, ymax)
    plt.xlabel("Reference data")
    plt.ylabel("MTP predict error")
    plt.savefig('figferr.' + format, dpi=300)
