#!/usr/bin/env python

"""
Coarse-grained protein in MARTINI.
"""

from amx import *

init()
make_step(settings.step)
if not os.path.isfile(state.start_structure): 
	raise Exception('cannot locate state.start_structure: %s'%state.start_structure)
shutil.copy(state.start_structure,state.here+'protein-start.pdb')
martinize()
