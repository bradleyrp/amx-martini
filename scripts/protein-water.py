#!/usr/bin/env python

"""
PROTEIN SIMULATION
Atomistic protein in water.
"""

from amx import *

init()
make_step(settings.step)
write_mdp()
#---we require protein_prepared from the coarse-graining step
shutil.copyfile(state.protein_prepared['top'],state.here+'vacuum.top')
shutil.copyfile(state.protein_prepared['gro'],state.here+'vacuum-alone.gro')
state.composition = state.protein_prepared['composition']
if state.itp==None: state.itp = []
for itp in state.martinize_itps:
	state.itp.append(os.path.basename(itp))
	shutil.copyfile(itp,state.here+state.itp[-1])
#---minimize, solvate, and equilibrate
write_top('vacuum.top')
gmx('editconf',
	structure='vacuum-alone',
	gro='vacuum',
	c=True,d='%.2f'%settings.water_buffer,
	log='editconf-vacuum-room')
minimize('vacuum',method='steep')
solvate(
	structure='vacuum-minimized',
	gro='solvate')
write_top('solvate.top')
minimize('solvate')
counterions(
	structure='solvate-minimized',
	top='solvate',
	ff_includes='ions')
move_file('counterions.gro','counterions_wrong_names.gro')
fix_martini_ions(structure='counterions_wrong_names.gro',gro='counterions.gro')
minimize('counterions')
copy_file('counterions-minimized.gro','system.gro')
grouper(ndx='system-groups',protein=True,lipids=False)
write_top('system.top')
equilibrate(groups='system-groups')
