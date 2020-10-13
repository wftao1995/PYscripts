import sys
import numpy as np
import spglib as spg

#@TODO: it seems that we need a chemical map
class cstruc:
    #positions in Direct
    pos = []
    #d or c, only affects writing POSCAR
    mode = 'C'
    #scale factor
    scale = 1.0
    #if select dynamics
    select = False
    #contains volocity
    volc = False
    #discription
    dispstr = ''
    #base vector
    basevec = np.zeros((3,3), dtype=float, order='C')
    has_atom_type = False
    #@TODO: merge specices and neachatoms into a dictionary
    specices = []
    #@TODO: remove this
    nspecices = 0
    neachatoms = []
    chmmap = []
    tatoms = 0
    #selective dynamics, 1 for Fix
    fix_mat = []
    #space group, using spglib
    spgroup = 0
    symmdataset = []
    symprec = 1e-5
    cell = () #for spglib
    layer_tags = [] # label the layer which each atom lies
    symm_tags = [] #symmtery tag for each atom, to see if two atoms are equavelent

    def readposcar(self, filename = 'POSCAR', volc = False):
        i = 1
        self.select = False
        self.volc = False
        self.has_atom_type = False
        with open(filename) as poscarf:
            for line in poscarf:
                if  line.strip() == '' or line.strip()[0] == '#':
                    continue
                if i == 1:
                    self.dispstr = line
                elif i == 2:
                    self.scale = float(line.strip())
                elif i <= 5:
                    #strip comments
                    istr = line.split()
                    self.basevec[i-3] = np.array(list(map(float, istr[0:3])))
                #now test if POSCAR has any elements name
                elif i == 6 or i == 7:
                    cl = line.strip().split()
                    nsp = 0
                    for j in cl:
                        if j[0] == '#':
                            break
                        nsp = nsp + 1
                    self.nspecices = nsp
                    try:
                        ncl = list(map(int,cl[0:nsp]))
                    except:
                        self.specices = cl[0:nsp]
                        self.has_atom_type = True
                    else:
                        self.neachatoms = ncl
                        self.tatoms = sum(ncl)
                        #self.pos = np.zeros((self.tatoms,3))
                        self.pos = []
                        #self.fix_mat = np.zeros((self.tatoms,3),dtype=int)
                        self.fix_mat = []
                        i = 7
                #if there any "Selective dynamics"
                elif i == 8 or i == 9:
                    if line.strip()[0].upper() == 'S':
                        self.select = True
                    elif line.strip()[0].upper() == 'C' or line.strip()[0].upper() == 'K':
                        self.mode = 'C'
                        i = 9
                    elif line.strip()[0].upper() == 'D':
                        self.mode = 'D'
                        i = 9
                    else:
                        raise Exception("无法确定坐标类型")
                #read atom positions
                elif i <= 9 + self.tatoms:
                    #self.pos[i-10] = np.array(list(map(float,line.strip().split()[0:3])))
                    self.pos.append(np.array(list(map(float,line.strip().split()[0:3]))))
                    if self.select:
                        cl = [0 if j == 'T' else 1 for j in line.strip().split()[3:6]]
                        #self.fix_mat[i-10] = np.array(cl)
                        self.fix_mat.append(cl)
                else:
                    #not to read extra data currently
                    break    
                i = i + 1
        #convert c to d
        if self.mode == 'C' or self.mode == 'K':
            for j in range(self.tatoms):
                self.pos[j] = np.dot(self.pos[j],np.linalg.inv(self.basevec))
        #build chmmap
        for i in range(self.nspecices):
            for j in range(self.neachatoms[i]):
                self.chmmap.append(self.specices[i])
        #build a dummy fix-mat
        if not self.select:
            self.fix_mat = np.zeros((self.tatoms,3),dtype=int).tolist()
    #@TODO improvement
    def buildsupercell(self, a):
        mp = a[0]*a[1]*a[2]
        #npos = np.zeros((self.tatoms*mp,3))
        #nfix_mat = np.zeros((self.tatoms*mp,3),dtype=int)
        npos = []
        nfix_mat = []
        nchmmap = []
        vc = np.array([0,0,1])
        vb = np.array([0,1,0])
        va = np.array([1,0,0])
        for i in range(3):
            self.basevec[i] = self.basevec[i] * a[i]
        #t = 0
        for s in range(self.tatoms):
            bpos = self.pos[s]
            for i in range(a[0]):
                for j in range(a[1]):
                    for k in range(a[2]):
                        #npos[t] = (bpos + i * va + j * vb + k * vc) / np.array(a)
                        npos.append((bpos + i * va + j * vb + k * vc) / np.array(a))
                        #nfix_mat[t] = self.fix_mat[s]
                        nfix_mat.append(self.fix_mat[s])
                        nchmmap.append(self.chmmap[s])
                        #t = t + 1
        self.neachatoms = [ i * mp for i in self.neachatoms ]
        self.tatoms = self.tatoms * mp
        self.pos = npos
        self.fix_mat = nfix_mat

    def write_poscar(self,filename = "POSCAR_NEW"):
        fp = open(filename,'w')
        fp.write(self.dispstr)
        fp.write(' ' + f'{self.scale:.10f}' + '\n')
        for i in range(3):
            ws = '    ' + ' '.join(f'{j:.10f}' for j in self.basevec[i]) + '\n'
            fp.write(ws)
        if self.has_atom_type:
            fp.write('      ' + ' '.join(i for i in self.specices) + '\n')
        fp.write('      ' + ' '.join(str(i) for i in self.neachatoms) + '\n')
        if self.select:
            fp.write('Selective Dynamics' + '\n')
            if self.mode == 'D':
                fp.write('Direct' + '\n')
                for i in range(self.tatoms):
                    dyn = [ 'T' if j == 0 else 'F' for j in self.fix_mat[i] ]
                    fp.write('    ' + ' '.join(f'{j:.10f}' for j in self.pos[i]) + ' ' + ' '.join(j for j in dyn) + '\n')
            else:
                fp.write('Carstain' + '\n')
                for i in range(self.tatoms):
                    dyn = [ 'T' if j == 0 else 'F' for j in self.fix_mat[i] ]
                    npos = np.dot(self.pos[i], self.basevec)
                    fp.write('    ' + ' '.join(f'{j:.10f}' for j in npos) + ' ' + ' '.join(j for j in dyn) + '\n')
        else:
            if self.mode == 'D':
                fp.write('Direct' + '\n')
                for i in range(self.tatoms):
                    fp.write('    ' + ' '.join(f'{j:.10f}' for j in self.pos[i]) + '\n')
            else:
                fp.write('Carstain' + '\n')
                for i in range(self.tatoms):
                    npos = np.dot(self.pos[i], self.basevec)
                    fp.write('    ' + ' '.join(f'{j:.10f}' for j in npos) + '\n')
        
    def setmode(self, mode):
        if mode.upper() == 'C' or mode.upper() == 'K':
            self.mode = 'C'
        else:
            self.mode = 'D'
        
    def set_atom_type(self,atomtypes):
        if self.has_atom_type:
            print("Warning: cell has a type info already : " + self.specices)
        self.has_atom_type = True
        self.specices = atomtypes

    def cell_for_sym(self):
        lattice = self.basevec
        #positions = self.pos.tolist()
        positions = list(map(lambda x : x.tolist(), self.pos))
        numbers = []
        for i in range(self.nspecices):
            for j in range(self.neachatoms[i]):
                numbers.append(i)
        self.cell = (lattice, positions, numbers)

    def find_sym(self):
        self.cell_for_sym()
        self.spgroup = spg.get_spacegroup(self.cell, symprec=self.symprec)
        self.symmdataset = spg.get_symmetry_dataset(self.cell, symprec=self.symprec)
    #@TODO: build primitive cell
    def find_primitive(self):
        if self.select:
            raise Exception("无法构建原胞，请解除限制")
        return spg.find_primitive(self.cell)
    #read a 3d array
    def redefine_lattice(self, newbase, need_right_hand = True):
        if need_right_hand:
            #test if cross product of A and B meets the direction of C
            if np.dot(np.cross(newbase[0],newbase[1]),newbase[2]) < 0:
                raise Exception("该取向非右手系，请将某个向量取相反方向!")
        #note that this is not a copy but just acts like "pointers"
        #npos = self.pos
        for i in range(self.tatoms):
            self.pos[i] = np.dot(self.pos[i],np.dot(self.basevec,np.linalg.inv(newbase)))
        self.basevec = newbase
    #surface is a tuple (a, b, c)
    #for hex, howto??
    #@TODO: finish this
    def find_surface_info(self,surface, tol = 1e-4):
        pass
    #平移，给定平移矢量，direct is a list， rearrange表示是否将超出边界的原子用晶胞内的镜像位置表示，rearrange_tol为容许超出边界的最大值
    def translate(self, direct, mode = 'D', rearrange = True, rearrange_tol = 0.01):
        dec = np.array(direct)
        if(mode !='D'):
            dec = np.dot(dec, np.linalg.inv(self.basevec))
        for i in range(self.tatoms):
            self.pos[i] = self.pos[i] + dec
            for j in range(3):
                while self.pos[i][j] > 1 + rearrange_tol:
                    self.pos[i][j] = self.pos[i][j] - 1
                while self.pos[i][j] < -rearrange_tol:
                    self.pos[i][j] = self.pos[i][j] + 1
    def group_by_element(self, sortby = []):
        nduplist = list(set(self.specices))
        dups = [x for x in self.specices if self.specices.count(x) > 1]
        nneachatoms = []
        for i in nduplist:
            if i in dups:
                dupids = [x for x in range(self.specices) if self.specices[x] == i]
                nt = sum([self.neachatoms[x] for x in dupids])
                nneachatoms.append(nt)
            else:
                nneachatoms.append(self.neachatoms[self.specices.index(i)])
        if len(set(sortby)) != len(sortby):
            raise Exception("Wrong list")
        if sortby == []:
            sortby = nduplist
        #@TODO: more efficient implementaion
        npos = self.pos.copy()
        nfix_mat = self.fix_mat.copy()
        nchmm = self.chmmap.count()
        index = list(range(self.tatoms))
        cmindex = []
        for i in sortby:
            for j in range(len(sortby)):
                if nduplist[j] == i:
                    cmindex.append(j)
            for j in range(self.tatoms):
                if self.chmmap[j] == i:
                    index[i] = j
        for i in range(self.tatoms):
            self.pos[i] = npos[index[i]]
            self.fix_mat[i] = nfix_mat[index[i]]
            self.chmmap[i] = nchmm[index[i]]
        self.specices = sortby
        self.neachatoms = []
        for i in range(len(sortby)):
            self.neachatoms.append(nneachatoms[cmindex[i]])

    #for intersitials
    def add_atom_by_position(self, position, symbol, mode = 'D', use_cell_symm = False):
        if mode == 'D':
            iposition = position
        else:
            iposition = np.dot(np.array(position), self.basevec)
        self.chmmap.append(symbol)
        if use_cell_symm:
            #@TODO: finish this
            pass
        else:
            self.tatoms = self.tatoms + 1
            self.pos.append(iposition)
            if symbol in self.specices:
                #@TODO: finish this
                findid = self.specices.index(symbol)
                self.neachatoms[findid] = self.neachatoms[findid] + 1
                self.group_by_element()
            else:
                self.specices.append(symbol)
                self.nspecices = self.nspecices + 1
                self.neachatoms.append(1)

    def remove_atom_by_id(self, atom_id, use_cell_symm = False):
        chm = self.chmmap[atom_id]
        cindex = self.specices.index(chm)
        self.neachatoms[cindex] = self.neachatoms[cindex] - 1
        if self.neachatoms[cindex] == 0:
            del self.neachatoms[cindex]
            del self.specices[cindex]
            self.nspecices = self.nspecices - 1
        del self.pos[atom_id]
        del self.fix_mat[atom_id]
        del self.chmmap[atom_id]
        self.tatoms = self.tatoms - 1

    def get_atom_id_by_position(self, position, prec = 1e-4):
        for i in range(self.tatoms):
            if np.linalg.norm(self.pos[i] - position) < prec:
                return i

    def remove_atom_by_position(self, position, use_cell_symm = False, prec = 1e-4):
        try:
            atomid = self.get_atom_id_by_position(position)
        except:
            raise Exception("No atom at given position")
        del self.pos[atomid]
        atomtype = self.chmmap[atomid]
        del self.chmmap[atomid]
        del self.fix_mat[atomid]
        self.tatoms = self.tatoms - 1
        typeid = self.specices.index(atomtype)
        self.neachatoms[typeid] = self.neachatoms[typeid] - 1
        if self.neachatoms[typeid] == 0:
            del self.neachatoms[typeid]
            del self.specices[typeid]
            self.nspecices = self.nspecices - 1
    
    def sub_atom_by_id(self, id, use_cell_symm = False):
        pass
    def sub_atom_by_position(self, position, use_cell_symm = False, prec = 1e-4):
        pass
    def remove_atoms_by_type(self, type):
        pass
    #@TODO: how about layer id?
    #@TODO: mode
    def apply_constrain_by_axis_value_range(self, z1 = 0, z2 = 1, constrain = [0,0,0], axis = 2, mode = 'D', default_constrain=[1,1,1]):
        zmax = max(z1,z2)
        zmin = min(z1,z2)
        if len(default_constrain) != 3 or axis > 2 or axis < 0:
            raise Exception("Invalid value!")
        if not self.select:
            self.select = True
            #1 for fix in fix_mat
            #self.fix_mat = np.zeros((self.tatoms,3),dtype=int)
            self.fix_mat = []
            for i in range(self.tatoms):
                #self.fix_mat[i] = default_constrain
                self.fix_mat.append(default_constrain)
        for i in range(self.tatoms):
            if self.pos[i][axis] > zmin and self.pos[i][axis] < zmax:
                self.fix_mat[i] = np.array(constrain)
    #axis : 0 for x, 1 for y and 2 for z
    #flag : false for ascen, true for descen 
    #algo : only 1 is avaliable yet
    def sort_atoms_by_axis(self, axis, flag = False, algo = 1):
        if algo != 1:
            raise Exception("Not implemented yet!")
        #需要排序指标
        index = list(range(self.tatoms))
        index.sort(key = lambda y: self.pos[y][axis], reverse = flag)
        posnew = self.pos.copy()
        fixmat = self.fix_mat.copy()
        chemm = self.chmmap.copy()
        for i in range(self.tatoms):
            self.pos[i] = posnew[index[i]]
            self.fix_mat[i] = fixmat[index[i]]
            self.chmmap[i] = chemm[index[i]]
        self.group_by_element()

    #position: leave, bottom, upper, center
    def build_vacuum_along_axis(self, axis = 2, size = 12, bulk_position = 'leave'):
        for i in range(self.tatoms):
            self.pos[i] = np.dot(self.pos[i], self.basevec)
        len0 = np.linalg.norm(self.basevec[axis])
        t = 1.0 + size/len0
        self.basevec[axis] = self.basevec[axis] * t
        #@TODO: finish this
        if bulk_position == 'bottom':
            pass
        #turn back into direct
        for i in range(self.tatoms):
            self.pos[i] = np.dot(self.pos[i],np.linalg.inv(self.basevec))