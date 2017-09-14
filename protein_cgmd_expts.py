{

'martinize':{
'tags':['cgmd','protein','tag_support','tested_2017.09.14'],
'script':'scripts/protein.py',
'params':None,
'extensions':['@martini/martini.py'],
'settings':"""

start structure: inputs/exo70-monomer-body.pdb
martinize path: @martini/bin/martinize.py
dssp path: @martini/bin/dssp-2.0.4-linux-amd64

step: protein
martinize flags: -ed
martinize ff: martini22

"""},

}