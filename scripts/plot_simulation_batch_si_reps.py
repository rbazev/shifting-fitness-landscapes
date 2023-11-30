import myutils
import glob
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from itertools import cycle
import argparse

# parse command-line arguments:
#####################################################
parser = argparse.ArgumentParser(description='plot population weighted specialization index trajectories for replicate simulations for a single model')
parser.add_argument('outputdir', type=str, help='main output directory')
parser.add_argument('simsetname', type=str, help='a simulation set name. expects multiple runname subdirectories within it')
parser.add_argument('batches', type=str, nargs='+', help='list of simulation batch/run names, expect within each a ./simdata')
args = parser.parse_args()

# customcolors = ['tab:pink', 'tab:cyan', 'tab:green','tab:purple', 'tab:brown', 'tab:red','tab:gray', 'skyblue','pink',  'tab:blue', 'tab:orange',  ]
# customcolors = ['tab:orange', 'tab:cyan', 'tab:pink',  'tab:green','tab:purple', 'tab:brown', 'tab:red','tab:gray', 'skyblue','pink',  'tab:blue',   ]
# custcomcolorscycle = cycle(customcolors)


for batch in args.batches:

	fig = plt.figure(figsize=(8, 1.5)) # 4.5, 4

	print("beginning processing simulation batch %s" % batch)
	if 'continuous' in batch: batch_label = 'continuous'
	if 'static' in batch: batch_label = batch[batch.find('static') : batch.find('_') ].replace("EvoC","Gen")
	if 'discrete' in batch: batch_label = 'discrete'
	batch_color = 'k'
	batch_dir = args.outputdir + '/simulations/%s/%s/simdata/' % (args.simsetname, batch)
	simfiles = glob.glob(batch_dir + '/sim-*.csv')

	lol_si = [] # list-of-lists of si. each list in this list-of-lists will contain the population average si trajectory of a simulation.

	for simfile in simfiles:
		print("extracting si trajectory for simfile %s" % simfile)
		simdf = pd.read_csv(simfile)

		# extract the list of si for this simulation
		l_si = list(simdf.loc[1:,"landscape"]) # "1:"" specifies all rows from generation 1 onward (0 is starting population), and "landscape" specifies the column label
		l_si = [entry.split(';')[0] for entry in l_si] # only keep the string to the left of the ';' containing 'si=...'
		l_si = [float(entry.split('=')[1]) for entry in l_si] # extract the si from the string 'si=...'
		lol_si.append(l_si)

	print("completed extracting si trajectory arrays for this batch.")

	for simulation_trajectory_si in lol_si:
		plt.plot(simulation_trajectory_si, color=batch_color, alpha = 0.25, linewidth=0.3)

	plt.gca().yaxis.set_tick_params(labelsize=8)
	plt.gca().xaxis.set_tick_params(labelsize=8)
	plt.gca().set_ylim(-1.1,1.1)
	plt.gca().set_title(batch, fontsize=8)
	plt.gca().set_xlabel("Generation", fontsize=8)
	plt.gca().set_ylabel("Pop. SI", fontsize=8)
	plt.gca().spines[['top', 'right']].set_visible(False)

	#plt.legend(fontsize=8)
	plt.tight_layout()
	plt.savefig(args.outputdir + '/simulations/%s/si_simreps-%s.pdf' % (args.simsetname, batch))
	plt.savefig(args.outputdir + '/simulations/%s/si_simreps-%s.png' % (args.simsetname, batch))

