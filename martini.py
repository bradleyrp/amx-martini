#!/usr/bin/env python

import os,glob

def martinize():
	"""
	Use martinize to generate a coarse-grained protein.
	"""
	name = 'protein'
	martinize_fn = state.martinize_path
	#---this function is run from the step but martinize_path is relative to root
	if not os.path.isfile(martinize_fn): raise Exception('cannot find martinize at %s'%martinize_fn)
	cmd = 'python '+os.path.abspath(martinize_fn)+' -v -p backbone '
	cmd += ' -f protein-start.pdb -o %s.top -x %s.pdb'%(name,name)
	if state.dssp_path: cmd += ' -dssp %s'%os.path.abspath(os.path.expanduser(state.dssp_path))
	if state.martinize_ff_version: cmd += ' -ff %s'%state.martinize_ff_version
	if state.martinize_flags: cmd += ' '+state.martinize_flags
	bash(cmd,cwd=state.here,log='martinize')
	if not os.path.isfile(state.here+'protein.pdb'): raise Exception('martinize failed')
	gmx_run(state.gmxpaths['editconf']+' -f %s.pdb -o %s.gro'%(name,name),log='editconf-convert-pdb')
	#---only allow Z-restraints because this is probably for a bilayer
	#---! WHY?
	bash("sed -i 's/POSRES_FC    POSRES_FC    POSRES_FC/0 0 POSRES_FC/g' Protein.itp",cwd=state.here)
	#---! only works for no-simulation because we pick up all ITP files in the folder for later
	#---! assume that all ITP files are relevant. you could also get itps directly from martinize.py?
	itps = [i for i in glob.glob(state.here+'/*.itp') if len(GMXTopology(i).molecules)>0]
	#---save these for later, particularly for the composition detector in the adhere_protein_bilayer
	state.martinize_itps = itps
