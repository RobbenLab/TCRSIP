#!/usr/bin/env python
# -*- coding: UTF8 -*-
"""
author: Guillaume Bouvier
email: guillaume.bouvier@ens-cachan.org
creation date: 2014 05 20
license: GNU GPL
Please feel free to use and modify this, but keep the above information.
Thanks!
"""

import sys
sys.path.append('/c5/shared/pymol/1.7.0.0-python-2.7.5-shared/lib/python2.7/site-packages/')
 
import __main__
__main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
import pymol
pymol.finish_launching()

pdb_file =sys.argv[1]
outdir =sys.argv[2]
pdb_name =pdb_file.split('.')[0]
pymol.cmd.load(pdb_file, pdb_name)
pymol.cmd.remove('hydrogens')
pymol.cmd.save(outdir)