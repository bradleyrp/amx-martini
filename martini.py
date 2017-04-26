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
	#---! write a separate position restrained version if desired by dropping -p None
	cmd += ' -f protein-start.pdb -o %s.top -x %s.pdb '%(name,name)
	#---dssp use should be standard for MARTINI proteins
	if state.dssp_path: 
		dssp_fn = os.path.abspath(os.path.expanduser(state.dssp_path))
		if not os.path.isfile(dssp_fn):
			raise Exception('cannot find %s'%dssp_fn)
		cmd += ' -dssp %s'%dssp_fn
	if state.martinize_ff_version: cmd += ' -ff %s'%state.martinize_ff_version
	if state.martinize_flags: cmd += ' '+state.martinize_flags
	bash(cmd,cwd=state.here,log='martinize')
	if not os.path.isfile(state.here+'protein.pdb'): raise Exception('martinize failed')
	gmx_run(state.gmxpaths['editconf']+' -f %s.pdb -o %s.gro'%(name,name),log='editconf-convert-pdb')
	#---we infer that the only ITP files are generated by martinize
	itps = glob.glob(state.here+'*.itp')
	#---only allow Z-restraints because this is probably for a bilayer
	#---we wish to carefully control the position restraints for fixing the bilayer so we cleave them from the ITP
	#---...previously this was done with bash and "sed -i 's/POSRES_FC    POSRES_FC    POSRES_FC/0 0 POSRES_FC/g"
	#---...however we used the "-p None" flag to accomplish this within martinize.py
	#---! only works for no-simulation because we pick up all ITP files in the folder for later
	#---! assume that all ITP files are relevant. you could also get itps directly from martinize.py?
	itps = [i for i in glob.glob(state.here+'/*.itp') if len(GMXTopology(i).molecules)>0]
	#---save these for later, particularly for the composition detector in the adhere_protein_bilayer
	state.martinize_itps = itps
	#---register this step with the state for subsequent steps
	state.protein_prepared = {'gro':state.here+'protein.gro','top':state.here+'protein.top'}
