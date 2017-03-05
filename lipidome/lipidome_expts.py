{

'lipidome':{
#####
####
###
##
#
'metarun':[
{'quick':'clear_lipidome'},
{'quick':'generate_lipidome_structures'},
{'quick':'generate_lipidome_restraints'}
]},

'clear_lipidome':{
#####
####
###
##
#
'quick':"""

from amx import *
import os,shutil,glob
for fn in (glob.glob(settings.auto_ff+os.sep+'*')+
	glob.glob(settings.structure_drop)):
	if os.path.isdir(fn): shutil.rmtree(fn)
	else: os.remove(fn)
if not os.path.isdir(settings.auto_ff): os.mkdir(settings.auto_ff)

""",'settings':"""

auto ff: inputs/martini/auto_ff
structure drop: inputs/martini/library-lipidome-structs

""",
},

'generate_lipidome_structures':{
#####
####
###
##
#
#
'tags':['cgmd','lipidome'],
'params':'@bilayers/parameters.py',
'extensions':[
	'structlib.py',
	'topology_maker.py',
	'@extras/geometry_tools/*.py'],
'quick':"""

from amx import *
init()
martini_lipidome()
write_martini_landscape()

""",
'settings':"""

USAGE NOTES:|
	clear the "deposit_at" directory above
	run "make clean sure && make prep generate_lipidome_structures && make run"
	this regenerates "crude" starting structures in the deposit at directory
	!!! need to update the landscape.yaml file
	!!! need to clear and regenerate the library whenever you want to add new lipids
	note that the `molecules` setting below controls the available lipids

#---source ITP with all MARTINI lipids
itp_source: inputs/martini/martini-sources.ff/martini_v2.0_lipids_all_201506.itp
molecules: "DOPC DOPS DOPE POP2 POP3 POPI POPC DPPC".split()

equilibration: []
mdp specs:| {'group':'cgmd','mdps':{
	'input-em-steep-in.mdp':['minimize'],
	'input-md-in.mdp':['invisible'],
	}}

force_field: martini_cv2
force_field_defs: martini-v2.2.itp
sources: ['inputs/martini/martini-sources.ff']

#---path to deposit gro and itp folders from root
deposit at:    inputs/martini/library-lipidome-structs
landscape at:  inputs/martini/auto_ff/landscape.json

"""},

'generate_lipidome_restraints':{
#####
####
###
##
#
'tags':['cgmd','lipidome'],
'params':'@bilayers/parameters.py',
'extensions':[
	'structlib.py',
	'topology_maker.py',
	'@extras/geometry_tools/*.py'],
'quick':"""

from amx import *
init()
topology_maker_martini_restraints()

""",
'settings':"""

USAGE NOTES:|
	this method is designed to pre-make any lipid restraints you might want
	its product can be used by 
		(1) the bilayer maker's vacuum packing
		(2) the flat bilayer maker's leaflet-specific restraints
	the "wants" dictionary specifies the outputs and describes their restraints
	the deposit site holds the automatically generated force fields
	this method completely avoids using the "define posre" flags in GROMACS
	all restraints are explicit, but this means you should avoid "define posre" which will restrain water
	!!! note: missing DPP2 and CHL1

base force field: inputs/martini/martini-sources.ff  # source force field to modify
deposit site: inputs/martini/auto_ff                 # where to write new force fields (keys in wants)

#---specify transformed force field copies
wants:|{
	'martini_upright.ff':{
		'restraints':{'martini_glycerol':{'z':100},'martini_tails':{'z':100}},
		'naming':'same','which':'lipids'
		},
	'martini_upright_alt.ff':{
		'restraints':{'martini_glycerol':{'z':1000}},
		'naming':'alternate','which':'lipids'
		},
	'martini_prison.ff':{
		'restraints':{'martini_glycerol':{'z':1000},'martini_tails':{'z':1000}},
		'naming':'alternate_restrain_both','which':'lipids'
		},
	}

"""},

}