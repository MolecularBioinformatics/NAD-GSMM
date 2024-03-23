#!usr/bin/python
import pandas as pd
import cobra
from os import listdir
import xml.etree.ElementTree as et
import re
from pathlib import Path


def read_sabiork(query_files=None, species=None):
	"""
	Reads a number of SabioRK files and adds them all into a single dataframe.
	:param query_files: List of paths to the SabioRK query files. Defaults to using the 'sabiork_queries' folder.
	:param species: List of strings for species to be used. Defaults to NAD(P) and NAD(P)H.
	:return: Dataframe containing all relevant data: UniProt ID, EC, species, KM value, unit
	"""
	if not query_files:
		query_path = Path('./sabiork_queries')
		query_files = [Path(f) for f in listdir(query_path)]
		query_files = ['./sabiork_queries/' / f for f in query_files]
	if not species:
		species = ['NAD+', 'NADH', 'NADP+', 'NADPH']
	dfs = []
	# create dataframes from all files and add them to a list
	for qf in query_files:
		df = pd.read_csv(qf, delimiter='\t', header=0, index_col=1)
		df = df[df['parameter.type'] == 'Km']
		df = df[df['parameter.associatedSpecies'].isin(species)]
		# set species names into lowercase:
		df.index = [x.lower() for x in list(df.index)]
		dfs.append(df)
	# vertically concatenate dataframes and rename:
	df = pd.concat(dfs)
	df.columns = ['id', 'up', 'ec', 'type', 'species', 'value', 'end_value', 'standard_dev', 'unit']
	# drop useless columns and empty rows:
	df.drop(['id', 'type', 'end_value', 'standard_dev'], axis=1, inplace=True)
	df = df[~df.value.isnull()]
	# convert units from M to mM:
	df['value'] *= 1000
	df['unit'] = 'mM'
	df.index.name = 'animal'
	
	return df
		

def read_brenda(query_files=None):
	"""
	Reads a number of Brenda files and adds them into a single dataframe.
	:param query_files: List of paths to the Brenda query files. Defaults to using the 'brenda_queries' folder.
	:return: Dataframe containing all relevant data: EC, species, KM value, unit
	"""
	if not query_files:
		query_files = ['./brenda_queries/' + f for f in listdir('./brenda_queries')]
	dfs = []
	# create dataframes from all files and add them to a list
	for qf in query_files:
		df = pd.read_csv(qf, delimiter='\t', index_col=4, header=None)
		df.columns = ['ec', 'value', 'species', 'commentary', 'ligand', 'literature']
		df.drop(['commentary', 'ligand', 'literature'], axis=1, inplace=True)
		# set species names into lowercase:
		df.index = [x.lower() for x in list(df.index)]
		dfs.append(df)
	# vertically concatenate dataframes and rename:
	df = pd.concat(dfs)
	# add indicator that units are in mM, no conversion necessary:
	df['unit'] = 'mM'
	df.index.name = 'animal'
	
	return df


def create_mappings_ec_mitocore(model, species_tag='nad', save = False):
	"""
	Creates mapping of mitocore reactions onto ECs. Saves into file and returns dictionary.
	:param model: Path to model file.
	:param species_tag: Tag to identify species by. Defaults to "NAD".
	:return: Dictionary of associations. Reaction IDs as keys.
	"""
	# todo: current version of COBRApy completely killed this by not parsing EC in notes
	model = cobra.io.read_sbml_model(model)
	rea_nad = [r for r in model.reactions if species_tag in r.build_reaction_string()]
	associations = {}
	err = 0
	for r in rea_nad:
		try:
			ec = r.notes['EC NUMBER']
			for s in [',', '+', 'and', 'or']:
				ec = [e.split(s) for e in ec]  # filter out wrong entries
				ec = [e.strip() for ex in ec for e in ex]  # flattens list
			ec = [e for e in ec if '.' in e]
			associations[r.id] = ec
		except KeyError:
			err += 1
	
	print('ECs found for {} reactions.'.format(len(associations)))
	print('No ECs found for {} reactions.'.format(err))
	
	if save:
		assoc_str = []
		for rxn in associations:
			assoc_str.append('{}\t{}'.format(rxn, '\t'.join(associations[rxn])))
		assoc_str = '\n'.join(assoc_str)
		with open('./mappings_mitocore_ec.txt', 'w') as f:
			f.write(assoc_str)
		
	return associations


def create_mappings_ec_mitocore_from_xml(modelFile, species_tag='nad', save_path = None, verbose = False):
	re_ec = re.compile('[0-9]\.[0-9]+\.[0-9]+\.[0-9]+')
	root = et.parse(modelFile).getroot() # parses XML
	model = root[0]
	rxns = model[4]
	
	mapping = {}
	# iterate over reactions in ListOfReactions tag
	nad_found = []
	nad_not_found = []
	for x in rxns:
		rxn_id = x.attrib['id'][2:]
		# find the notes section
		note = x[3][0]
		
		# if reactions don't involve NAD we want to skip them
		nad = False
		for p in note:
			if 'Description' in p.text:
				if (species_tag.upper() in p.text) or (species_tag.lower() in p.text):
					nad = True
					nad_found.append(rxn_id)
		# skip the reaction if NAD is not involved
		if not nad:
			nad_not_found.append(rxn_id)
			continue
		# finally find ECs
		error_count = 0
		
		ec_found = False
		for p in note:
			if 'EC Number' in p.text:
				ec_found = True
				try:
					ecs = re_ec.findall(p.text)
					if not ecs: 
						print(f'No mapping found for {rxn_id}')
						continue
					for ec in ecs:
						try:
							mapping[rxn_id].append(ec)
						except KeyError:
							mapping[rxn_id] = [ec]
				except IndexError: # no ECs found
					print(f'Invalid {p.text} in {rxn_id}')
					error_count += 1
		if verbose and not ec_found:
			print(f'No EC annotation in reaction {rxn_id}')
	if verbose:
		print(f'NAD found in {len(nad_found)} reactions:')
		print(nad_found)
		print()
		print(f'NAD not found in {len(nad_not_found)} reactions:')
		print(nad_not_found)
		print()

	print('ECs found for {} reactions.'.format(len(mapping)))
	print('No ECs found for {} reactions.'.format(error_count))
	
	return mapping


def load_mappings(file_name):
	"""
	Loads associations from a tsv file.
	:param file_name: Path to file with associations.
	:return: Dictionary containing associations.
	"""
	associations = {}
	with open(file_name, 'r') as f:
		for line in f:
			line = line.split()
			associations[line[0]] = line[1:]
			
	return associations
