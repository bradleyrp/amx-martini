{

'martinize':{
#####
####
###
##
#
'tags':['cgmd','protein','tag_support','tested_2017.09.14'],
'script':'scripts/protein-martini.py',
'params':None,
'extensions':['@martini/martini.py'],
'settings':"""

USE NOTES:|
	This is a support function that typically precedes a protein-water or protein-bilayer simulation.
	It requires a local copy of DSSP find in this module.
	It effectively wraps the "martinize.py" script and receives any extra settings fia "martinize flags".
	Hence this interaface to MARTINI should accomodate almost all coarse-graining settings.

step: protein                                       # name of this step
start structure: inputs/exo70-monomer-body.pdb      # starting structure
martinize path: @martini/bin/martinize.py           # location of the "martinize.py" script we are wrapping
dssp path: @martini/bin/dssp-2.0.4-linux-amd64      # get secondary structure from DSSP binary located here
martinize ff: martini22                             # which MARTINI force field to use (via martinize.py -ff)
martinize flags: -ed                                # see martinize.py help/docs for all settings here

"""},

'protein_water_martini':{
#####
####
###
##
#
'tags':['cgmd','tag_support','tested_2017.10.24.1700_dev'],
'script':'scripts/protein-water.py',
'params':'@bilayers/parameters.py',
'extensions':['@martini/martini.py'],
'settings':"""

USE NOTES:|
	This is a support function that should always follow a topology-generate step e.g. "martinize".
	Note that we get MDP parameters from the bilayer module.
	It receives information from e.g. martinize which loads the state with e.g. protein_prepared.
	Needs: interchangeable water model, standardized wrapper for a coarse-grained model (a long-term goal)

step: protein                       # name of the folder is s01-protein
equilibration: nvt-short,nvt,npt    # which equilibration steps to use (requires MDP files from mdp_specs)
water buffer: 1.2                   # distance to the wall (nm) which sets the volume of water
ionic strength: 0.150               # desired molar ionic strength
cation: NA+                         # name of the cation for neutralizing the system
anion: CL-                          # name of the anion for neutralizing the system
sol: W                              # resname for MARTINI water

files:       ['@structure-repo/martini/martini-water.gro']   # copy the standard water box 
sources:     ['@martini/martini-sources.ff']    # collect the standard martini force field
force field:  martini-sources                   # specify the name of the force field (minus ".ff" suffix)
solvent:      martini-water                     # solvent box for MARTINI (omit the suffix, copy via files)

#---INTEGRATOR PARAMETERS generated via parameters.py
mdp_specs:| {
	'group':'cgmd',
	'mdps':{
		'input-em-steep-in.mdp':['minimize'],
		'input-em-cg-in.mdp':['minimize',{'integrator':'cg'}],
		'input-md-nvt-eq-in.mdp':['nvt-protein','nvt-protein',{'nsteps':10000}],
		'input-md-nvt-short-eq-in.mdp':['nvt-protein',{'dt':0.001,'nsteps':10000}],
		'input-md-npt-eq-in.mdp':['npt-protein',{'nsteps':10000}],
		'input-md-in.mdp':['npt-protein-production',{'nsteps':100000}],},}

"""},

'demo_helix0_water':{
#####
####
###
##
#
'tags':['cgmd','note_structure_repo','tested_2017.10.24.1700_dev'],
'metarun':[
{'step':'martinize','do':'martinize','settings':"""
# note that this test uses a very small starting structure so if you get domain problems, make it bigger
start structure: @structure-repo/proteins/helix0.pdb
"""},
{'step':'protein','do':'protein_water_martini','settings':"""
water buffer: 1.2
"""},]},

}