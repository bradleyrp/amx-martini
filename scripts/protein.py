#!/usr/bin/env python

"""
Coarse-grained protein in MARTINI.
"""

from amx import *

init()
make_step(settings.step)
shutil.copy(state.start_structure,state.here+'protein-start.pdb')
martinize()
