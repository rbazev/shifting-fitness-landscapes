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
parser = argparse.ArgumentParser(description='plot the virus phenotype outcomes at the end of evolution simulations for various models')
parser.add_argument('outputdir', type=str, help='main output directory')
parser.add_argument('simsetname', type=str, help='a simulation set name, will have multiple runname subdirectories within it')
parser.add_argument('pdfprefix', type=str, help='prefix for pdf of figure')
parser.add_argument('batches', type=str, nargs='+', help='list of simulation batch/run names, each is a directory containing a ./simdata')
args = parser.parse_args()

customcolors = ['tab:pink', 'tab:cyan', 'tab:green','tab:purple', 'tab:brown', 'tab:red','tab:gray', 'skyblue','pink',  'tab:blue', 'tab:orange',  ]

# specify file paths containing genotypes that are L or O specialists:
Lspecs = myutils.read_selectionrate_file(args.outputdir + '/selectionrate/specialists/Lspecialists.csv', as_dict=True)
Ospecs = myutils.read_selectionrate_file(args.outputdir + '/selectionrate/specialists/Ospecialists.csv', as_dict=True)
# define genotypes that are neither L-specialists or O-specialists:
nonspecs = [gt for gt in myutils.construct_genotype_space(length=9) if gt not in Lspecs and gt not in Ospecs]

endpoint_groups_by_batch = {} # keyed by batch_label eg 'continuous' or 'static-vEvoC-mixed' etc; 
	# then keyed by 'single', 'double', 'triple' for number of phenotypes. value is *number of sims ending with that number of phenotypes*.
	# another key 'dualspec' has value of number of sims ending with a combination of L and O specialist.
	# OR key by N, L, O, NL, NO, LO, NLO

batch_labels = []

for i,batch in enumerate(args.batches): 
	print("beginning processing of simulation batch (model) %s" % batch)
	if 'continuous' in batch: batch_label = 'continuous'
	if 'static' in batch: batch_label = batch[batch.find('static') : batch.find('_')].replace("EvoC","Gen")
	if 'discrete' in batch: batch_label = 'discrete'
	batch_labels.append(batch_label)

	#endpoint_groups_by_batch[batch_label] = {'single':0, 'double':0, 'triple':0, 'dualspec':0} # for this batch/model, count these resulting phenotype combinations.
	endpoint_groups_by_batch[batch_label] = {'G':0, 'L':0, 'O':0, 'GL':0, 'GO':0, 'LO':0, 'GLO':0}

	batch_dir = args.outputdir + '/simulations/%s/%s/simdata/' % (args.simsetname,batch)
	simfiles = glob.glob(batch_dir + '/sim-*.csv')
	numsims = len(simfiles)

	cum_endpoints_by_gt = {} #k: genotype v: # of simulations that genotype ends w/ frequency >= 5%
	for gt in myutils.construct_genotype_space(length=9): cum_endpoints_by_gt[gt] = 0

	for simfile in simfiles:
		#print("parsing endpoint genotypes in simulation %s" % simfile)
		simdf = pd.read_csv(simfile)

		# parse genotypes above cutoff at final timepoint:
		tf = simdf.index[-1]
		total_tf = simdf.loc[tf, '[000000000]':'[111111111]'].sum() # total abundance at time t final
		a_densities_tf = simdf.loc[tf, '[000000000]':'[111111111]'] / float(total_tf) # array of genotype densities at time t final
		endpoint_winners = a_densities_tf[a_densities_tf >= 0.025]
		# endpoint_winners is a list of genotypes that are >= 2.5% frequency at end of sim. 
		# now save cumulative endpoints,
		# and figure out if these winners are a combination of phenotypes or not.
		Lspec_present = False
		Ospec_present = False
		Nonspec_present = False

		for winner in endpoint_winners.index:
			cum_endpoints_by_gt[winner] += 1
			if winner in Lspecs:
				Lspec_present = True
			elif winner in Ospecs:
				Ospec_present = True
			else:
				Nonspec_present = True

		num_phenotypes_present = sum([Lspec_present,Ospec_present,Nonspec_present])
		#print("There are %s phenotypes present in this simulation" % num_phenotypes_present)

		if num_phenotypes_present == 1:
			if Lspec_present: key = 'L'
			elif Ospec_present: key = 'O'
			else: key = 'G'
		elif num_phenotypes_present == 2:
			if Lspec_present and Ospec_present: key = 'LO'
			elif Lspec_present: key = 'GL'
			elif Ospec_present: key = 'GO'
		elif num_phenotypes_present == 3: key = 'GLO'
		else: raise phenotypeserror

		endpoint_groups_by_batch[batch_label][key] += 1


	print("completed parsing endpoints for this batch of simulations.")

	# save a file with the cumulative number of simulations each genotype was found at the simulation endpoint
	cum_endpoints_file = open(batch_dir + '/cumulativeEndpoints.csv','w')
	cum_endpoints_file.write('genotype,cumulativeEndpoints\n')
	for gt in myutils.construct_genotype_space(length=9):
		cum_endpoints_file.write('%s,%s\n'%(gt,cum_endpoints_by_gt[gt]))
	cum_endpoints_file.close()

	print(endpoint_groups_by_batch[batch_label])


# adjacent barplots:

if len(args.batches) == 6:
	nrows = 2
	ncols = 3
	figsize = (4.5,4)
else:
	nrows = 1
	ncols = len(args.batches)

fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)

if nrows > 1:
	axes = axes.ravel() # axs.ravel() converts 2-dim axes object to a list in row-major order

plt.rcParams['hatch.linewidth'] = 7
hatch = '/' # '///'

color_assignments = {'G':'darkgray','L':'tab:blue','O':'tab:red'}
label_assignments = {'G': 'generalist', 'L':'L-specialist', 'O':'O-specialist'}

for batch_i, batch_label in enumerate(batch_labels):

	# G/L/O/GL/GO/LO/GLO
	xticks = [1,2,3,4,5,6] # eg., ax.set_xticks([2,4,6,8,10]) # 7
	xticklabels = ['G','L','O','G+L','G+O','L+O'] #eg.,  ax.set_xlabels([‘two’, ‘four’,’six’, ‘eight’, ‘ten’]) # 'G+L+O'

	for xtick, data_label in zip(xticks,['G','L','O','GL','GO','LO']): # GLO

		if len(data_label) == 1:
			color = color_assignments[data_label] 
			label = label_assignments[data_label]
			axes[batch_i].bar(xtick,endpoint_groups_by_batch[batch_label][data_label],width=0.7, color=color, edgecolor=None, lw=0, label=label)

		elif len(data_label) == 2:
			facecolor = color_assignments[data_label[0]]
			hatchcolor = color_assignments[data_label[1]]
			axes[batch_i].bar(xtick,endpoint_groups_by_batch[batch_label][data_label],width=0.7, color=facecolor, hatch=hatch, edgecolor=hatchcolor, lw=0)
		else:
			axes[batch_i].bar(xtick,endpoint_groups_by_batch[batch_label][data_label],width=0.7, color='k', lw=0)

	axes[batch_i].set_xticks(xticks)
	axes[batch_i].set_ylim(0,numsims)
	axes[batch_i].yaxis.set_tick_params(labelsize=8)
	axes[batch_i].set_xticklabels(xticklabels, fontsize=8, rotation=90)
	axes[batch_i].spines[['top', 'right']].set_visible(False)
	axes[batch_i].set_title(batch_label, fontsize=8)

fig.supxlabel('receptor phenotypes present at endpoint', fontsize=8)
fig.supylabel('number of simulations',fontsize=8)
plt.legend(fontsize=6, bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.tight_layout()
plt.savefig(args.outputdir + '/simulations/%s/endpointphenotypes-%s.pdf' % (args.simsetname,args.pdfprefix))
plt.savefig(args.outputdir + '/simulations/%s/endpointphenotypes-%s.png' % (args.simsetname,args.pdfprefix))
plt.close()
