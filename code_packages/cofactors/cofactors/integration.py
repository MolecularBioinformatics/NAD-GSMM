import pandas as pd
import cobra
from copy import deepcopy
from statistics import median, mean

TOLERANCE = 10**(-5)


def _get_km(kms, reaction, mappings, id_type, decision, preference, verbose=False):
	"""
	Pulls the KM for a given reaction.
	:param kms: Dataframe containing SabioRK or Brenda data.
	:param reaction: Reaction in question.
	:param mappings: Dictionary associating reactions with their UniProt or EC identifiers.
	:param id_type: Type of ID used. Either 'ec' or 'up'
	:param decision: How to choose one of the KMs. 'min', 'avg' or 'med'. Defaults to 'min'.
	:param preference: List of scientific species names. How to rank the different species. First one is most preferred.
	:return: A KM. Float.
	"""
	try:
		ids = mappings[reaction]
	except KeyError:
		return None
	if verbose:
		print(f'Getting Kms for {ids}')
	kms = kms[kms[id_type].isin(ids)]
	organisms = list(kms.index.unique())
	organisms = [o.lower() for o in organisms]
	preference = [o.lower() for o in preference]
	if kms.empty: # empty dataframes not being caught would lead to weird errors
		return None
	if verbose:
		print('All Kms mapped:')
		print(kms)
	for org in preference:
		if org in organisms:
			kms = kms[kms.index == org]
			kms = list(kms['value'])
			kms = [float(x) for x in kms]
			if verbose:
				print(f'Organism {org} found')
				print(f'Picking from Kms {sorted(kms)}')
			break
	#todo: Implement alternative decisions
	return decision(kms)


def create_models(model, mappings, kms, c_old_dict, c_new_dict, id_type='ec', decision=min, preference=None):
	"""
	Creates a model for the given NAD concentrations. Uses pFBA to update fluxes.
	:param model: Model object to be used
	:param mappings: Mappings of reactions to ECs or Uniprot IDs.
	:param kms: Dataframe containing SabioRK or Brenda data.
	:param c_old_dict: Dictionary of compartment to [NAD].
	:param c_new_dict: Dictionary of compartment to [NAD].
	:param id_type: Identifier to use. 'ec' or 'up'.
	:param decision: How to decide between multiple KMs. Defaults to min (picking minimal Km).
	:param preference: List of scientific species names. How to rank the different species. First one is most preferred.
	:return: Model with adapted boundaries.
	"""
	model_new = deepcopy(model)
	if not preference:
		preference = ['Homo sapiens', 'Sus scrofa', 'Bos taurus', 'Rattus norvegicus', 'Mus musculus']
	sol = cobra.flux_analysis.parsimonious.pfba(model_new)
	
	for rxn in model_new.reactions:
		km = _get_km(kms, rxn.id, mappings, id_type, decision, preference)
		if km is None:
			continue
		comp = rxn.get_compartments()[0]
		c_old = c_old_dict[comp]
		c_new = c_new_dict[comp]
		flux = sol.fluxes[rxn.id]
		try:
			f_goal = flux * (c_new / (km + c_new)) * ((km + c_old) / c_old)  # set an ideal new flux using the old flux
		except ZeroDivisionError:
			raise ZeroDivisionError(f"ZeroDivisionError in {rxn.id} with KM {km}, c' {c_new} and c {c_old}")
		if f_goal > 0:
			rxn.upper_bound = f_goal
			rxn.lower_bound = 0
		else:
			rxn.upper_bound = 0
			rxn.lower_bound = f_goal
	
	return model_new


def create_models_fva(model, mappings, kms, c_old_dict, c_new_dict, id_type='ec', decision=min, preference=None, obj_frac=.9, verbose=False):
	"""
	Creates a model for the given NAD concentrations. Uses flux variance to update fluxes.
	:param model: Model object to be used
	:param mappings: Mappings of reactions to ECs or Uniprot IDs.
	:param kms: Dataframe containing SabioRK or Brenda data.
	:param c_old_dict: Dictionary of compartment to [NAD].
	:param c_new_dict: Dictionary of compartment to [NAD].
	:param id_type: Identifier to use. 'ec' or 'up'.
	:param decision: How to decide between multiple KMs. Defaults to min (picking minimal Km).
	:param preference: List of scientific species names. How to rank the different species. First one is most preferred.
	:param obj_frac: Fraction of optimal objective flux that needs to be retained. Float, defaults to 0.9 .
	:param verbose: Print additional output for information.
	:return: Model with adapted boundaries.
	"""
	model_new = deepcopy(model)
	if not preference:
		preference = ['Homo sapiens', 'Sus scrofa', 'Bos taurus', 'Rattus norvegicus', 'Mus musculus']
	sol = cobra.flux_analysis.flux_variability_analysis(model_new, fraction_of_optimum=obj_frac)
	not_adjusted_rxns = []
	adjusted_rxns = []
	for rxn in model_new.reactions:
		if verbose:
			print(f'Adjusting reaction {rxn.id}')
		km = _get_km(kms, rxn.id, mappings, id_type, decision, preference, verbose=verbose) # mM in the Brenda and SabioRK DBs
		if km is None:
			not_adjusted_rxns.append(rxn.id)
			if verbose:
				print('No Km found. Skipping {rxn.id}\n---\n')
			continue
		if verbose:
			print(f'Picked Km {km}')
		comp = rxn.get_compartments()[0]
		if verbose:
			print(f'Compartment used is {comp}')
		c_old = c_old_dict[comp]
		c_new = c_new_dict[comp]
		if verbose:
			print(f'Old and new concentrations are: {c_old} -- {c_new}')
		flux_low = sol.loc[rxn.id, 'minimum']
		flux_high = sol.loc[rxn.id, 'maximum']
		if verbose:
			print(f'FVA bounds are {flux_low:.2f} -- {flux_high:.2f}')
		try:
			f_goal_high = flux_high * (c_new / (km + c_new)) * ((km + c_old) / c_old)
			f_goal_low = flux_low * (c_new / (km + c_new)) * ((km + c_old) / c_old)
		except ZeroDivisionError:
			raise ZeroDivisionError(f"ZeroDivisionError in {rxn.id} with KM {km}, c' {c_new} and c {c_old}")
		
		rxn.upper_bound = float(max(0, f_goal_high))
		rxn.lower_bound = float(min(f_goal_low, 0))
		if verbose:
			print(f'New bounds are {rxn.lower_bound:.2f} -- {rxn.upper_bound:.2f}')
			print('---\n')
		adjusted_rxns.append(rxn.id)
	if verbose:
		print(f'{len(adjusted_rxns)} reactions adjusted:')
		print(adjusted_rxns)
		print()
		print(f'{len(not_adjusted_rxns)} reactions NOT adjusted:')
		print(not_adjusted_rxns)

	
	return model_new

def create_models_bounds(model, mappings, kms, c_old_dict, c_new_dict, id_type='ec', decision='min', preference=None, obj_frac=.9):
	"""
	Creates a model for the given NAD concentrations. Uses flux variance to update fluxes.
	:param model: Model object to be used
	:param mappings: Mappings of reactions to ECs or Uniprot IDs.
	:param kms: Dataframe containing SabioRK or Brenda data.
	:param c_old_dict: Dictionary of compartment to [NAD].
	:param c_new_dict: Dictionary of compartment to [NAD].
	:param id_type: Identifier to use. 'ec' or 'up'.
	:param decision: How to decide between multiple KMs. Defaults to min (picking minimal Km).
	:param preference: List of scientific species names. How to rank the different species. First one is most preferred.
	:param obj_frac: Fraction of optimal objective flux that needs to be retained. Float, defaults to 0.9 .
	:return: Model with adapted boundaries.
	"""
	model_new = deepcopy(model)
	if not preference:
		preference = ['Homo sapiens', 'Sus scrofa', 'Bos taurus', 'Rattus norvegicus', 'Mus musculus']
	
	for rxn in model_new.reactions:
		km = _get_km(kms, rxn.id, mappings, id_type, decision, preference)
		if km is None:
			continue
		comp = rxn.get_compartments()[0]
		c_old = c_old_dict[comp]
		c_new = c_new_dict[comp]
		ub = rxn.upper_bound
		lb = rxn.lower_bound
		try:
			f_goal_high = ub * (c_new / (km + c_new)) * ((km + c_old) / c_old)
			f_goal_low = lb * (c_new / (km + c_new)) * ((km + c_old) / c_old)
		except ZeroDivisionError:
			raise ZeroDivisionError(f"ZeroDivisionError in {rxn.id} with KM {km}, c' {c_new} and c {c_old}")
		
		rxn.upper_bound = max(0, f_goal_high)
		rxn.lower_bound = min(f_goal_low, 0)
	
	return model_new

def scan_reaction(model, rxn_id, mappings, kms, c_old_dict, fva=False, id_type='ec', decision='min', preference=None):
	"""
	Checks reaction activity and reduced cost over a wide range of NAD concentrations
	:param model: Model object to be used.
	:param rxn_id: ID of reaction to be tracked.
	:param mappings: Mappings of reactions to identifiers.
	:param kms: Dataframe containing Brenda or SabioRK data.
	:param c_old_dict: Orginal concentrations. Dict with compartments as keys and concentrations as values.
	:param fva: Boolean. Use FVA based integration. Defaults to False.
	:param id_type: Type of identifier 'ec' or 'up'. Defaults to 'ec'.
	:param decision: How to choose between multiple KMs. Defaults to 'min' (using the minimal Km found).
	:param preference: Order of preference of species. List of scientific species names.
	:return: List of reaction fluxes.
	"""
	fluxes = []
	reduced_costs = []
	for i in range(1, 99):
		ratio = i / 100
		c_new_dict = {}
		for comp in c_old_dict:
			c_new_dict[comp] = c_old_dict[comp] * ratio
		
		if fva:
			model_new = create_models_fva(model, mappings, kms, c_old_dict, c_new_dict, id_type, decision, preference)
		else:
			model_new = create_models(model, mappings, kms, c_old_dict, c_new_dict, id_type, decision, preference)
		sol_p = cobra.flux_analysis.parsimonious.pfba(model_new)
		sol = model_new.optimize()
		fluxes.append(sol_p.fluxes[rxn_id])
		reduced_costs.append(sol.reduced_costs[rxn_id])
	
	return fluxes, reduced_costs


def scan_all_fluxes(model,  mappings, kms, c_old_dict, fva=False, id_type='ec', decision='min', preference=None):
	"""
	Checks reaction activity and reduced cost over a wide range of NAD concentrations
	:param model: Model object to be used.
	:param rxn_id: ID of reaction to be tracked.
	:param mappings: Mappings of reactions to identifiers.
	:param kms: Dataframe containing Brenda or SabioRK data.
	:param c_old_dict: Orginal concentrations. Dict with compartments as keys and concentrations as values.
	:param fva: Boolean. Use FVA based integration. Defaults to False.
	:param id_type: Type of identifier 'ec' or 'up'. Defaults to 'ec'.
	:param decision: How to choose between multiple KMs. Defaults to 'min' (using the minimal Km found).
	:param preference: Order of preference of species. List of scientific species names.
	:return: List of reaction fluxes.
	"""
	fluxes = []
	reduced_costs = []
	dfs = []
	for i in range(1, 99):
		ratio = i / 100
		c_new_dict = {}
		for comp in c_old_dict:
			c_new_dict[comp] = c_old_dict[comp] * ratio
		
		if fva:
			model_new = create_models_fva(model, mappings, kms, c_old_dict, c_new_dict, id_type, decision, preference)
		else:
			model_new = create_models(model, mappings, kms, c_old_dict, c_new_dict, id_type, decision, preference)
		sol_p = cobra.flux_analysis.parsimonious.pfba(model_new)
		dfs.append(sol_p.fluxes)
	dfs = pd.concat(dfs, axis=1, sort=True)
	dfs.columns = list(range(1, 99))
	
	return dfs


def scan_total_flux(model, mappings, kms, c_old_dict, fva=False, id_type='ec', decision='min', preference=None):
	"""
	Checks reaction activity for a range of different NAD concentrations.
	:param model: Model object to be used.
	:param mappings: Mappings of reactions to identifiers.
	:param kms: Dataframe containing Brenda or SabioRK data.
	:param c_old_dict: Orginal concentrations. Dict with compartments as keys and concentrations as values.
	:param fva: Boolean. Use FVA based integration. Defaults to False.
	:param id_type: Type of identifier 'ec' or 'up'. Defaults to 'ec'.
	:param decision: How to choose between multiple KMs. Defaults to 'min' (using the minimal Km found).
	:param preference: Order of preference of species. List of scientific species names.
	:return: List of total flux values.
	"""
	fluxes = []
	for i in range(1, 99):
		ratio = i / 100
		c_new_dict = {}
		for comp in c_old_dict:
			c_new_dict[comp] = c_old_dict[comp] * ratio
		
		if fva:
			model_new = create_models_fva(model, mappings, kms, c_old_dict, c_new_dict, id_type, decision, preference)
		else:
			model_new = create_models(model, mappings, kms, c_old_dict, c_new_dict, id_type, decision, preference)
		sol = cobra.flux_analysis.parsimonious.pfba(model_new)
		fluxes.append(sol.objective_value)
	
	return fluxes


def scan_all_reactions(model, mappings, kms, c_old_dict, fva=False, id_type='ec', decision='min', preference=None):
	"""
	Checks reaction activity and reduced cost over a wide range of NAD concentrations
	:param model: Model object to be used.
	:param mappings: Mappings of reactions to identifiers.
	:param kms: Dataframe containing Brenda or SabioRK data.
	:param c_old_dict: Orginal concentrations. Dict with compartments as keys and concentrations as values.
	:param fva: Boolean. Use FVA based integration. Defaults to False.
	:param id_type: Type of identifier 'ec' or 'up'. Defaults to 'ec'.
	:param decision: How to choose between multiple KMs. Defaults to 'min' (using the minimal Km found).
	:param preference: Order of preference of species. List of scientific species names.
	:return: List of reaction fluxes.
	"""
	fluxes = {}
	reduced_costs = {}
	rxns = [r.id for r in model.reactions]
	for r in rxns:
		fluxes[r] = []
		reduced_costs[r] = []
	for c in range(1, 99):
		ratio = c / 100
		c_new_dict = {}
		for comp in c_old_dict:
			c_new_dict[comp] = c_old_dict[comp] * ratio
		
		if fva:
			model_new = create_models_fva(model, mappings, kms, c_old_dict, c_new_dict, id_type, decision, preference)
		else:
			model_new = create_models(model, mappings, kms, c_old_dict, c_new_dict, id_type, decision, preference)
		
		sol_p = cobra.flux_analysis.parsimonious.pfba(model_new)
		sol = model_new.optimize()
		for r in rxns:
			fluxes[r].append(sol_p.fluxes[r])
			reduced_costs[r].append(sol.reduced_costs[r])
	
	fluxes = pd.DataFrame(fluxes)
	#fluxes.index = [r.id for r in model.reactions]
	reduced_costs = pd.DataFrame(reduced_costs)
	#reduced_costs.index = [r.id for r in model.reactions]
	
	return fluxes, reduced_costs


def scan_reaction_variance(model, rxn_id, mappings, kms, c_old_dict, tol=.99, fva=False, id_type='ec', decision='min', preference=None):
	"""
	Checks reaction activity and reduced cost over a wide range of NAD concentrations
	:param model: Model object to be used.
	:param rxn_id: ID of reaction to be tracked.
	:param mappings: Mappings of reactions to identifiers.
	:param kms: Dataframe containing Brenda or SabioRK data.
	:param c_old_dict: Original concentrations. Dict with compartments as keys and concentrations as values.
	:param tol: Fractional optimality tolerated  for FVA.
	:param fva: Boolean. Use FVA based integration. Defaults to False.
	:param id_type: Type of identifier 'ec' or 'up'. Defaults to 'ec'.
	:param decision: How to choose between multiple KMs. Defaults to 'min' (using the minimal Km found).
	:param preference: Order of preference of species. List of scientific species names.
	:return: List of reaction fluxes.
	"""
	upper = []
	lower = []
	for i in range(1, 99):
		ratio = i / 100
		c_new_dict = {}
		for comp in c_old_dict:
			c_new_dict[comp] = c_old_dict[comp] * ratio
		
		if fva:
			model_new = create_models_fva(model, mappings, kms, c_old_dict, c_new_dict, id_type, decision, preference)
		else:
			model_new = create_models(model, mappings, kms, c_old_dict, c_new_dict, id_type, decision, preference)
		
		sol_fva = cobra.flux_analysis.flux_variability_analysis(model_new, [rxn_id], fraction_of_optimum=tol)
		upper.append(sol_fva.loc[rxn_id,  'maximum'])
		lower.append(sol_fva.loc[rxn_id, 'minimum'])
	
	return upper, lower


def scan_all_variance(model, mappings, kms, c_old_dict, tol=.99, fva=False, id_type='ec', decision='min', preference=None):
	"""
	Checks reaction activity and reduced cost over a wide range of NAD concentrations
	:param model: Model object to be used.
	:param mappings: Mappings of reactions to identifiers.
	:param kms: Dataframe containing Brenda or SabioRK data.
	:param c_old_dict: Original concentrations. Dict with compartments as keys and concentrations as values.
	:param tol: Fractional optimality tolerated  for FVA.
	:param fva: Boolean. Use FVA based integration. Defaults to False.
	:param id_type: Type of identifier 'ec' or 'up'. Defaults to 'ec'.
	:param decision: How to choose between multiple KMs. Defaults to 'min' (using the minimal Km found).
	:param preference: Order of preference of species. List of scientific species names.
	:return: List of reaction fluxes.
	"""
	values = {}
	reactions = [r.id for r in model.reactions]
	for r in reactions:
		values[r] = []
	for i in range(1, 99):
		ratio = i / 100
		c_new_dict = {}
		for comp in c_old_dict:
			c_new_dict[comp] = c_old_dict[comp] * ratio
		
		if fva:
			model_new = create_models_fva(model, mappings, kms, c_old_dict, c_new_dict, id_type, decision, preference)
		else:
			model_new = create_models(model, mappings, kms, c_old_dict, c_new_dict, id_type, decision, preference)
		
		sol_fva = cobra.flux_analysis.flux_variability_analysis(model_new, fraction_of_optimum=tol)
		for r in reactions:
			values[r].append((sol_fva.loc[r, 'maximum'], sol_fva.loc[r, 'minimum']))
	
	return pd.DataFrame(values)
