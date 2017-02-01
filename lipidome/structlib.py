#!/usr/bin/env python

import json
import os,sys,shutil
import tempfile
import numpy as np
import amx

#---! note that you can explicitly use amx in these libraries if you want

###---ORIGINAL, ENERGY MINIMIZE METHOD HERE

"""
Collect the MARTINI lipidome for use in AUTOMACS.

peudocode:
	import all ITP files
	group ITP
	for each ITP
		make a temporary directory for each ITP
		get the number of points and add them to a GRO
		run a short minimization in vacuum
"""

#---! make sure the error induced by np.zeros() is obvious. weird that it's only in ast

#---! started developing here but that was dump because it couldn't get the sidewash
#---! then moved to structlib.py, but those functions assumed amx was in amx and not in globals
#---! then figured it out and removed "amx." from the extension module

def straighten_lipid(structure,gro='protein-flat.gro',direction='z',cwd='./'):
	"""
	Assume the first atom should be positive.
	Heavily adapted from the lipid-handling routine in lib_place_proteins.py.
	"""
	struct = amx.GMXStructure(structure)
	lpts = struct.points
	#---define the reference axis
	refaxis = np.zeros(3)
	refaxis['xyz'.index(direction)] = 1
	#---rotate the lipid to the reference axis
	lpts -= np.mean(lpts,axis=0)
	#---take the first principal component as the axis
	eigs = np.linalg.eig(np.dot(lpts.T,lpts))
	principal_axis_index = np.argsort(eigs[0])[-1]
	axis = vecnorm(eigs[1][:,principal_axis_index])
	#---sometimes lipids are flipped upside down so we choose the option that brings the lipid
	#---...closest to the "lipid_top" atom
	xyzs = np.dot(lpts,rotation_matrix(vecnorm(np.cross(refaxis,axis)),
		np.arccos(np.dot(refaxis,axis)))) 
	xyzs_b = np.dot(lpts,rotation_matrix(vecnorm(np.cross(refaxis,-1*axis)),
		np.arccos(np.dot(refaxis,-1*axis)))) 
	#---we assume that the "top" of the lipid is the first bead in the structure
	lipid_top = 0
	xyzs_b = np.dot(lpts,rotation_matrix(vecnorm(np.cross(refaxis,-1*axis)),
		np.arccos(np.dot(refaxis,-1*axis)))) 
	candidates = [xyzs,xyzs_b]
	xyzs = candidates[np.argmin([np.linalg.norm(i) 
		for i in [j[lipid_top]-j.mean(axis=0)-refaxis 
		for j in candidates]])]
	struct.points = xyzs-xyzs.mean(axis=0)
	struct.write(gro)

def generate_simple_structure_via_em(name,itps):
	"""
	INCORRECT METHOD
	Use energy minimization to generate a simple structure from an ITP.
	"""
	me = itps.molecules[name]
	amx.make_step(name)
	amx.write_mdp()
	n_atoms = len(me['atoms'])
	#######coords[:,0] = 0.4*np.arange(n_atoms)
	coords = np.random.random_sample((n_atoms,3))*0.4
	repack = {
		'pts':coords,
		'atom_names':np.array([i['atom'] for i in me['atoms']]),
		'residue_names':np.array([i['resname'] for i in me['atoms']]),
		'residue_indices':np.ones((n_atoms)),
		'box':np.array([10.,10.,10.]),
		}
	struct = amx.GMXStructure(**repack)
	struct.write(amx.state.here+'init.gro')
	amx.component(name,count=1)
	amx.write_top('vacuum.top')
	amx.minimize('init',top='vacuum')
	copy_file('em-init-steep.gro','vacuum.gro')
	amx.equilibrate(structure='vacuum',top='vacuum')
	amx.minimize('md.part0001',top='vacuum')
	straighten_lipid(structure=amx.state.here+'init.gro',gro=amx.state.here+name+'.gro')

###---EXPLICIT METHOD

#---magic for local import when you run from elsewhere
sys.path.insert(0,os.path.dirname(os.path.relpath(os.path.abspath(__file__),os.getcwd())))
from structlib_defs import defs_lipids,defs_ptdins,defs_pts_lipids,defs_pts_ptdins

def guess_lipid_coords(name,itps):
	"""
	"""
	me = itps.molecules[name]
	amx.make_step(name)
	amx.write_mdp()
	n_atoms = len(me['atoms'])
	#---parse the definitions for lipids only 
	defs = dict([(group_name,
		dict([(j[0],j[1:]) for j in [i.split() for i in d.splitlines()] if len(j)>0]))
		for group_name,d in [('lipids',defs_lipids),('ptdins',defs_ptdins)]])
	defs_pts = dict([('lipids',defs_pts_lipids),('ptdins',defs_pts_ptdins)])
	#---select the group for this lipid
	group_has = [key for key in defs if name in defs[key]]
	if not group_has: raise Exception('cannot find lipid %s in the definitions'%name)
	elif len(group_has)>1: raise Exception('lipid %s in more than one definition set'%name) 
	else: group_name = group_has[0]
	#---select the correct group
	defs = defs[group_name]
	defs_pts = defs_pts[group_name]
	if name not in defs: raise Exception('lipid %s not in definitions for group %s'%(name,group_name))
	#---get coordinates from standard definitions
	#---note that this information was published in the "insane" paper
	base_pts = np.array(defs_pts).T
	#---subselect the base points by atom names
	subsel = []
	for i in itps.molecules[name]['atoms']:
		if i['atom'] not in defs[name]: 
			raise Exception('cannot find %s from molecule %s in name list %s:\n%s'%(
				i['atom'],name,defs[name],str(itps.molecules[name])))
		subsel.append(defs[name].index(i['atom']))
	coords = base_pts[subsel]*0.25
	#---center the points at zero (the bilayer maker handles the box)
	coords = coords-coords.mean(axis=0)
	repack = {
		'pts':coords,
		'atom_names':np.array([i['atom'] for i in me['atoms']]),
		'residue_names':np.array([i['resname'] for i in me['atoms']]),
		'residue_indices':np.ones((n_atoms)),
		'box':np.array([10.,10.,10.]),
		}
	struct = amx.GMXStructure(**repack)
	out_fn = amx.state.here+name+'.gro'
	struct.write(out_fn)
	return out_fn

def martini_lipidome():
	"""
	Generate structures for all martini Lipids.
	"""
	#---prepare a directory in the inputs/martini repo to deposit these
	if os.path.isdir(amx.state.deposit_at): raise Exception('refusing to overwrite %s'%amx.state.deposit_at)
	else: os.mkdir(amx.state.deposit_at)
	itp_source = amx.state.q('itp_source',None)
	if not itp_source: raise Exception('settings must include itp_source, the path to a lipidome ITP file')
	itp_source = os.path.abspath(os.path.expanduser(itp_source))
	if not os.path.isfile(itp_source): raise Exception('cannot find ITP %s'%itp_source)
	itps = amx.GMXTopology(itp_source)
	mols = itps.molecules.keys()
	mols = amx.state.q('molecules',itps.molecules.keys())
	missing_mols = [i for i in mols if i not in itps.molecules]
	if any(missing_mols): raise Exception('missing molecules from the ITP: %s'%missing_mols)
	#---loop over each molecule and generate a structure
	fns = [guess_lipid_coords(mol,itps) for mol in mols]
	#---copy paths from assumed step folder numbering (better to be more explicit but this should work)
	for fn in fns: shutil.copyfile(fn,os.path.join(amx.state.deposit_at,os.path.basename(fn)))

def write_martini_landscape():
	"""
	"""
	#---note that the following is hardcoded and lipids are autodetected
	land = {'objects':{},'alias':{}}
	land['alias']['protein'] = ['GLH','ILE','ALAD','GLUH','GLN','HISH','ASN1','HYP','GLY','HIP',
		'ARGN','MSE','CYS1','GLU','CYS2','CYS','HISE','ASP','SER','HSD','HSE','PRO','CYX','ASPH',
		'ORN','HSP','HID','HIE','LYN','DAB','ASN','CYM','HISD','VAL','THR','HISB','HIS','HIS1',
		'HIS2','TRP','HISA','ACE','ASH','CYSH','PGLU','LYS','PHE','ALA','QLN','MET','LYSH','NME',
		'LEU','ARG','TYR']
	#---note that MARTINI generates ION,NA+ (this is handled in lib_place_proteins.detect_composition)
	land['objects'] = {
		'NA+':{'charge':1,'is':'ion','parts':['resname'],'ffs':['martini']},
		'CL-':{'charge':-1,'is':'ion','parts':['resname'],'ffs':['martini']},
		'W':{'charge':0,'is':'water','parts':['resname'],'ffs':['martini']},}
	for mol in state.molecules:
		if mol in land['objects']: raise Exception('molecule %s is already in the landscape'%mol)
		land['objects'][mol] = {'is':'lipid','parts':['resname'],'ffs':['martini']}
	with open(state.landscape_at,'w') as fp: fp.write(json.dumps(land))
