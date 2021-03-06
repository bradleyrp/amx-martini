#!/usr/bin/env python

import os,glob

def martinize():
	"""
	Use martinize to generate a coarse-grained protein.
	Note that this is the original implementation, but it has been refined in Multimer.py for more
	general use-cases.
	"""
	name = 'protein'
	martinize_fn = state.martinize_path
	#---this function is run from the step but martinize_path is relative to root
	if not os.path.isfile(martinize_fn): raise Exception('cannot find martinize at %s'%martinize_fn)
	#---we retain backbone position restraints which are activated by the posre flag from the MDP flag
	cmd = 'python '+os.path.abspath(martinize_fn)+' -v -p backbone '
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
	#---moved any control over position restraints to the "flat" procedures
	#---! assume that all ITP files are relevant. you could also get itps directly from martinize.py?
	itps = [i for i in glob.glob(state.here+'/*.itp') if len(GMXTopology(i).molecules)>0]
	#---save these for later, particularly for the composition detector in the adhere_protein_bilayer
	state.martinize_itps = itps
	#---register this step with the state for subsequent steps
	state.protein_prepared = {'gro':state.here+'protein.gro','top':state.here+'protein.top'}
	#---! how is martinize() used in the bilayer adhesion example, and how does this compare to the water?
	#---! hardcoded composition for now, for use in protein-water.py
	state.protein_prepared.update(composition=[('Protein_A',1)])

def fix_martini_ions(structure,gro):
	"""
	MARTINI likes to use NA+ and ION but the counterion functions in GROMACS are belligerent.
	"""
	struct = GMXStructure(state.here+structure)
	import numpy as np
	ions = np.where(struct.residue_names=='NA+')[0]
	struct.residue_names[ions] = 'ION'
	struct.atom_names[ions] = 'NA+'
	ions = np.where(struct.residue_names=='CL-')[0]
	struct.residue_names[ions] = 'ION'
	struct.atom_names[ions] = 'CL-'
	struct.write(state.here+gro)
