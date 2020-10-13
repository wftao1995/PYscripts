#!/bin/env python

from analysecfg import CFG
import os

a = CFG()
a.read_poscar('POSCAR')
fp = open('POSCAR.cfg','w')
a.write(fp)
