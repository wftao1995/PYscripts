#!/usr/bin/env python3

import sys
import analysecfg as cfg

a = cfg.CFGS()
a.readcfg(sys.argv[1])
a.strip_unconverge()
a.write_cfgs(sys.argv[2])


