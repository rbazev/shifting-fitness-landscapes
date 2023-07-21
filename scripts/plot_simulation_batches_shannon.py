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
parser = argparse.ArgumentParser(description='plot shannon entropy trajectories for replicate simulations under various models')
parser.add_argument('outputdir', type=str, help='main output directory')
parser.add_argument('simsetname', type=str, help='a simulation set name. expects multiple runname subdirectories within it')
parser.add_argument('pdfprefix', type=str, help='prefix for pdf of figure')
parser.add_argument('batches', type=str, nargs='+', help='list of simulation batch/run names, expect within each a ./simdata')
args = parser.parse_args()

customcolors = ['tab:pink', 'tab:cyan', 'tab:green','tab:purple', 'tab:brown', 'tab:red','tab:gray', 'skyblue','pink',  'tab:blue', 'tab:orange',  ]
customcolors = ['tab:orange', 'tab:cyan', 'tab:pink',  'tab:green','tab:purple', 'tab:brown', 'tab:red','tab:gray', 'skyblue','pink',  'tab:blue',   ]

custcomcolorscycle = cycle(customcolors)

# percentile boundaries for plotting around the medians:
lower_bound_shading = 25
upper_bound_shading = 75

fig = plt.figure(figsize=(4.5, 4))

for i,batch in enumerate(args.batches): 
	print("beginning processing of simulation batch %s" % batch)

	if 'continuous' in batch: batch_label = 'continuous'
	if 'static' in batch: batch_label = batch[batch.find('static') : batch.find('_') ].replace("EvoC","Gen")
	if 'discrete' in batch: batch_label = 'discrete'

	batch_color = next(custcomcolorscycle)
	batch_dir = args.outputdir + '/simulations/%s/%s/simdata/' % (args.simsetname,batch)
	simfiles = glob.glob(batch_dir + '/sim-*.csv')
	loa_diversity = [] # list-of-arrays of diversity. each array in this list is the diversity trajectory of a simulation.
	for simfile in simfiles:
		print("calculating diversity array for simfile %s" % simfile)
		simdf = pd.read_csv(simfile)
		a_diversity = myutils.compute_diversity(simdf, metric='shannon') # returns an array of length(generations) computing the diversity metric at each generation
		loa_diversity.append(a_diversity)

	print("completed calculating diversity arrays for this batch.")

	print("now calculating median and percentile intervals")
	aggregate_medians = []
	aggregate_lower_bound = []
	aggregate_upper_bound = []

	for time in range(len(loa_diversity[0])):
		aggregate_diversity_t = [] # list of diversities at time t (ie array index 'time')
		for diversity_array in loa_diversity: # for each simulation's array
			aggregate_diversity_t.append(diversity_array[time])

		aggregate_medians.append( np.median(aggregate_diversity_t) )
		lower,upper = np.percentile(aggregate_diversity_t, [lower_bound_shading,upper_bound_shading])
		aggregate_lower_bound.append(lower)
		aggregate_upper_bound.append(upper)

	#plot median diversity and shaded upper/lower percentiles
	times = range(len(loa_diversity[0]))
	plt.fill_between(times, aggregate_lower_bound, aggregate_upper_bound,
                 facecolor=batch_color, # The fill color
                 color=batch_color,       # The outline color
                 alpha=0.2)          # Transparency of the fill
	plt.plot(aggregate_medians, color=batch_color, label=batch_label)

plt.gca().yaxis.set_tick_params(labelsize=8)
plt.gca().xaxis.set_tick_params(labelsize=8)
plt.gca().set_xlabel("Generation", fontsize=8)
plt.gca().set_ylabel("Genetic diversity (Shannon entropy)", fontsize=8)
plt.gca().spines[['top', 'right']].set_visible(False)

plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig(args.outputdir + '/simulations/%s/shannon_aggregate-%s.pdf' % (args.simsetname,args.pdfprefix))
plt.savefig(args.outputdir + '/simulations/%s/shannon_aggregate-%s.png' % (args.simsetname,args.pdfprefix))