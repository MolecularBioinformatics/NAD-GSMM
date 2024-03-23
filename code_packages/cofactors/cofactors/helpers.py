import matplotlib.pyplot as pp


def read_subsystems(file_path):
	"""
	Reads subsystems from csv file into dict.
	:param file_path: path to csv file
	:return: Dictionary of subsystems matching reaction id (key) to subsystem (val).
	"""
	subs = {}
	with open(file_path, 'r') as f:
		for line in f:
			lin = line.split(',')
			lin[0] = lin[0][2:]
			subs[lin[0]] = lin[1]
	return subs


def plot_flux_with_var(vars, fluxes, rxn_id, savefile=None):
	"""
	Plots variability and fluxes of a scan of 1-99% concentration, displaying fluxes as graph, variability as error bars.
	:param vars: flux variability, data frame matching tuples to reactions
	:param fluxes: flux measured at dfferent concentrations for one reactions (list)
	:param rxn_id: ID of reaction
	:param savefile: Optional. File to save graph to.
	:return: No return value.
	"""
	rxn_var = list(vars[rxn_id])
	atp_low = [x[1] for x in rxn_var]
	atp_low = [abs(x1 - x2) for (x1, x2) in zip(atp_low, fluxes[0])]
	atp_high = [x[0] for x in rxn_var]
	atp_high = [abs(x1 - x2) for (x1, x2) in zip(atp_high, fluxes[0])]
	yerr = [atp_low, atp_high]
	
	x = range(1, 99)
	y = fluxes[0]
	fig, ax = pp.subplots()
	ax.errorbar(x, y, xerr=0, yerr=yerr, ecolor='r')
	pp.show()
	
	if savefile:
		fig.savefig(savefile)
	
	return fig
