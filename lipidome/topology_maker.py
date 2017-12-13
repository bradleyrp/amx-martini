#!/usr/bin/env python

_not_reported = ['deepcopy']

import os,sys,json,glob,shutil,re
from copy import deepcopy

#---! put this in a central place
argsort = lambda seq : [x for x,y in sorted(enumerate(seq), key = lambda x: x[1])]

def transform_itp(itp_fn,specs,**kwargs):
	"""
	Modify an ITP file.
	"""
	defs = kwargs.pop('defs',{})
	if kwargs: raise Exception('remaining kwargs %s'%kwargs)
	#---restraints options listing
	specs_add_restraints = 'which naming restraints'.split()
	#---each restraints option listing has a different method for changing restraints
	if set(specs.keys())>=set(specs_add_restraints):
		itp = GMXTopology(itp_fn,defs=defs)
		for mol in list(itp.molecules.keys()):
			mol_spec = itp.molecules[mol]

			#---note that the "which":"lipids" checks the tags for each ITP in the meta.json to see if there 
			#---...are any lipids in that ITP (but it might be the combined one, so there are other things)
			#---...this is why we perform a secondary check here

			#---if which lipids, we check whether this molecule is a lipid
			if specs['which']=='lipids':
				#---almost all MARTINI lipids have a GL1
				index_gl1 = [ii for ii,i in enumerate(mol_spec['atoms']) if i['atom']=='GL1']
				atoms_tail_1 = [ii for ii,i in enumerate(mol_spec['atoms']) if re.match('^.+A$',i['atom'])]
				atoms_tail_2 = [ii for ii,i in enumerate(mol_spec['atoms']) if re.match('^.+B$',i['atom'])]
				#---ensure this discriminates sterols properly
				is_sterol = mol=='CHOL'
				is_lipid = any(index_gl1) and any(atoms_tail_1) and any(atoms_tail_2)
			else: raise Exception('under development')
			#---if we are  transforming a lipid, apply martini_glycerol and martini_tails restraints
			if is_lipid or is_sterol:
				#---generate blank restraints
				posres_custom = {'funct': '1','fcy':'0','ai':'1','fcx':'0','fcz':'0'}
				#---hard-coding martini naming rules for head and tail atoms
				posres_all = [dict(posres_custom) for i in range(len(mol_spec['atoms']))]
				if not is_sterol:
					#---we already have tail names, so we get the last one
					tail_names = [sorted([mol_spec['atoms'][j]['atom'] 
						for j in k])[-1] for k in [atoms_tail_1,atoms_tail_2]]
					#---apply tail restraints if necessary
					if 'martini_glycerol' in specs['restraints']:
						atom_name = 'GL1'
						atoms = [i['atom'] for i in mol_spec['atoms']]
						for key,value in specs['restraints']['martini_glycerol'].items():
							if key not in 'xyz': 
								raise Exception('restraints, martini_glycerol keys must be in xyz')
							posres_all[atoms.index(atom_name)]['fc%s'%key] = value
					if 'martini_tails' in specs['restraints']:
						for atom_name in tail_names: 
							atoms = [i['atom'] for i in mol_spec['atoms']]
							for key,value in specs['restraints']['martini_glycerol'].items():
								if key not in 'xyz': 
									raise Exception('restraints, martini_glycerol keys must be in xyz')
								posres_all[atoms.index(atom_name)]['fc%s'%key] = value
				#---off for now
				elif False:
					#---! hard-coded MARTINI cholesterol restraints
					#---! CHECK THIS! IS IT CORRECT?
					if 'martini_sterol_top' in specs['restraints']:
						for atom_name in ['ROH','C2']: 
							atoms = [i['atom'] for i in mol_spec['atoms']]
							for key,value in specs['restraints']['martini_sterol_top'].items():
								if key not in 'xyz': 
									raise Exception('restraints, martini_sterol_top keys must be in xyz')
								posres_all[atoms.index(atom_name)]['fc%s'%key] = value
					#---! ONE IF LOOP? OR MANY?
					if 'martini_sterol' in specs['restraints']:
						for atom_name in ['ROH','C2']: 
							atoms = [i['atom'] for i in mol_spec['atoms']]
							for key,value in specs['restraints']['martini_sterol'].items():
								if key not in 'xyz': 
									raise Exception('restraints, martini_sterol_top keys must be in xyz')
								posres_all[atoms.index(atom_name)]['fc%s'%key] = value
				#---in some cases we want to restrain both the alternate group and the original
				#---...so we apply position restraints to the original here, and later copy it to 
				#---...the alternate
				if specs['naming']=='alternate_restrain_both':
					itp.molecules[mol]['position_restraints'] = posres_all	
				if specs['naming']=='same': mol_name = mol
				elif specs['naming'] in ['alternate','alternate_restrain_both']: 
					if len(mol)>4: raise Exception('name %s is too long. cannot make an alternate'%mol)
					mol_name = mol+'R'
					if mol_name in itp.molecules:
						raise Exception(
							'alternate for %s is %s but that is already in the molecules list'%(mol,mol_name))
					#---in the alternate scheme we make two copies of the lipid, 
					#---...and the "R"-suffixed has the restraints
					itp.molecules[mol_name] = deepcopy(itp.molecules[mol])
					#---apply the name in the correct entries
					#---this is important since the molname and resname are different items
					itp.molecules[mol_name]['moleculetype']['molname'] = mol_name
					for a in itp.molecules[mol_name]['atoms']:
						a['resname'] = mol_name
				else: raise Exception('unclear naming scheme')
				#---apply the position restraints we generated above
				itp.molecules[mol_name]['position_restraints'] = posres_all
	else: raise Exception('transform_itp cannot processes specs keys: "%s"'%specs.keys())
	return itp

def topology_maker_martini_restraints():
	"""
	Change something about a topology.
	"""
	#---! temporary setup for martini restraints only, but worthwhile to generalize

	#---settings
	base_ff = state.base_force_field
	if not os.path.isdir(base_ff): raise Exception('base force field must be a directory: %s'%base_ff)
	if not os.path.isfile(os.path.join(base_ff,'meta.json')): 
		raise Exception('base force field requires a meta.json file: %s'%base_ff)

	#---get the metadata for the incoming force field
	with open(os.path.join(base_ff,'meta.json')) as fp: meta = json.load(fp)

	#---interpret instructions from the settings
	for out_ff,specs in state.wants.items():

		#---we deposit the new force field in the deposit site
		out_ff_dn = os.path.join(state.deposit_site,out_ff)
		if os.path.isdir(out_ff_dn): raise Exception('refusing to overwrite %s'%out_ff_dn)
		else: os.mkdir(out_ff_dn)

		#---make sure the meta.json file describes all ITP files in the directory
		if set(meta.keys())!=set([os.path.basename(i) for i in glob.glob(os.path.join(base_ff,'*.itp'))]):
			import ipdb;ipdb.set_trace()
			raise Exception('keys in meta.json from base topology %s do not match the itp files there'%
				base_ff)

		#---apply rules to all itps
		for itp_name in meta.keys():
			#---select itps from meta.json tags and the which keyword
			if specs['which'] in meta[itp_name]:
				#---modify itp according to rules
				kwargs = dict()
				if 'defs' in specs: kwargs.update(defs=specs['defs'])
				itp_new = transform_itp(os.path.join(base_ff,itp_name),specs,**kwargs)
				#---the naming scheme only modifies the file internally. it is always rewritten
				itp_new.write(os.path.join(out_ff_dn,itp_name))
			#---ignore itps with tags in meta.json that are not equal to specs['which']
			else: shutil.copyfile(os.path.join(base_ff,itp_name),os.path.join(out_ff_dn,itp_name))
		#---retain the meta.json
		#---! without changes?
		shutil.copyfile(os.path.join(base_ff,'meta.json'),os.path.join(out_ff_dn,'meta.json'))
	print('[WARNING] just made a new topology but make sure that your atoms are correctly modified!')