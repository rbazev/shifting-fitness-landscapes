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
parser = argparse.ArgumentParser(description='plot simulated phage populations and parse aggregate metrics over multiple simulations')
parser.add_argument('outputdir', type=str, help='main output directory (needs to contain ./selectionrate/specialists/*specialists.csv)')
parser.add_argument('simulationdir', type=str, help='directory containing the simulation .log (.csv should be in ./simdata/)')
parser.add_argument('nrows', type=int, default=4, help='number of subplot rows')
parser.add_argument('ncols', type=int, default=4, help='number of subplot cols')
parser.add_argument('--minimumabundance', type=float, default=1e8, help='threshold max genotype abundance across entire sim to plot genotype')
parser.add_argument('--legendtype', type=str, default='single', help='single or multiple legends')
args = parser.parse_args()


# parse some info for this batch of simulations
#####################################################

simfiles = glob.glob(args.simulationdir + '/simdata/sim*.csv')
assert args.nrows*args.ncols <= len(simfiles)
myutils.makedir(args.simulationdir + '/pdfs/')

#parse the log file to annotate the plot with various parameters
simlogfile = args.simulationdir + '/log.txt'
with open(simlogfile) as f:
	for (i,line) in enumerate(f):
		if 'landscape_trajectory' in line:
			landscape_trajectory = line.strip('\n')[21:]
		elif 'noisepergen' in line:
			noisepergen = line.strip('\n')
		elif 'noisepersim' in line:
			noisepersim = line.strip('\n')
		elif 'startingpop' in line:
			starting_pop = line.strip('\n')

# define some things important in the analysis/plotting
#####################################################

# map a color cycle to genotypes 
# this has the advantage of providing each genotype a discernable color, as opposed to colormappings where adjacent genotypes have very similar colors
# sometimes, some colors can be repeated; if an issue, may have to expand or contract the cycle.
customcolors = ['tab:gray', 'tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:cyan','skyblue','pink']
custcomcolorscycle = cycle(customcolors)
mapped_custom_colors = {}
for gt in myutils.construct_genotype_space(length=9):
	mapped_custom_colors[gt] = next(custcomcolorscycle)

# import L and O specialist genotype identities from prior analysis
Lspecs = myutils.read_selectionrate_file(args.outputdir + '/selectionrate/specialists/Lspecialists.csv', as_dict=True)
Ospecs = myutils.read_selectionrate_file(args.outputdir + '/selectionrate/specialists/Ospecialists.csv', as_dict=True)
si_dict = myutils.read_selectionrate_file(args.outputdir + '/selectionrate/specialists/si.csv', as_dict=True)

# plot population curves for nrows*ncols simulations 
# as representative examples of the batch of all simulations.
#####################################################
#colortypes = ['colorbygenotype-plasma','colorbygenotype-jet','colorbygenotype-gist_rainbow','colorbygenotype-prism','customcolorcycle','colorbyspecialist']
colortypes = ['customcolorcycle','colorbysi']
for colortype in colortypes:
	fig, axes = plt.subplots(args.nrows, args.ncols, figsize=(5.5*args.ncols,2.5*args.nrows))
	# iterate sims and plot abundances to subplot axes
	for i, axis in enumerate(axes.flat):
		df = pd.read_csv(simfiles[i])
		for gt in myutils.construct_genotype_space(length=9):
			if df[gt].max() > args.minimumabundance:
				if colortype == 'colorbysi':
					si = si_dict[gt]
					cmap = mpl.cm.get_cmap('coolwarm_r')
					norm = mpl.colors.Normalize(vmin=-1, vmax=1)
					this_color = cmap(norm(si))
				elif colortype == 'colorbyspecialist':
					this_color = 'tab:blue' if gt in Lspecs else 'tab:orange' if gt in Ospecs else 'darkgray'
				elif 'colorbygenotype' in colortype:
					cmapname = colortype[16:]
					cmap = mpl.cm.get_cmap(cmapname, 512) # 'gist_rainbow''plasma''tab20''tab20b''jet'...
					this_color = cmap(myutils.gt_to_L9index[gt])[:3]
				elif colortype == 'customcolorcycle':
					this_color = mapped_custom_colors[gt]
				else:
					raise colorerror

				axis.plot(df['Time'],df[gt], color=this_color, label=gt, linewidth=1.5, alpha=0.8) 
				
				if args.legendtype == 'multiple':
					leg = axis.legend(loc="upper right", fontsize=6)
					for legobj in leg.legendHandles:
						legobj.set_linewidth(2.0)
				# if annotate: for t in df['Time']: if t>0 and df.loc[t][gt] > 1e8 and df.loc[t-1][gt] < 1e8:
				axis.spines['top'].set_visible(False)
				axis.spines['right'].set_visible(False)
		axis.axhline(y=0, color='k', linestyle=':')

		# # easy plotting of all genotypes that are observed at least some number of times *cumulatively*:
		# df.drop([col for col, val in df.sum().items() if val < (6e9/1000) and col != 'Time'], axis=1, inplace=True)
		# melted = pd.melt(df, ['Time'], var_name='genotype', value_name='abundance') # var_name = genotype, value_name = abundance 
		# g = sns.lineplot(x='Time', y='abundance', hue='genotype', ax=axis, data=melted)
		# g.legend_.remove()
		# axis.set_xlim(0,250)

	plt.suptitle("landscape trajectory: %s\n%s\n%s\n%s" % (landscape_trajectory,starting_pop,noisepergen,noisepersim), fontsize=10)
	
	if args.legendtype == 'single':
		labels_handles = {
	 	 label: handle for ax in fig.axes for handle, label in zip(*ax.get_legend_handles_labels())
		}
		leg = fig.legend(
		  labels_handles.values(),
		  labels_handles.keys(),
		  bbox_to_anchor=(1,1),
		  loc='upper right',
		  fontsize=6,
		  #bbox_transform = plt.gcf().transFigure,
		)
		for legobj in leg.legendHandles:
			legobj.set_linewidth(2.0)

	plt.tight_layout()
	plt.savefig(args.simulationdir+'/pdfs/%s.pdf' % colortype)
	plt.savefig(args.simulationdir+'/%s.png' % colortype)

	plt.close()




# similar form as the above but plotting diversity metrics instead of population abundances

for diversitymetric in ['shannon']: 
	fig, axes = plt.subplots(args.nrows, args.ncols, figsize=(5.5*args.ncols,2.5*args.nrows))
	# iterate sims and plot abundances to subplot axes
	for i, axis in enumerate(axes.flat):
		df = pd.read_csv(simfiles[i])
		a_diversity = myutils.compute_diversity(df, metric=diversitymetric) # returns an array of length(generations) computing the diversity metric at each generation
		axis.plot(a_diversity, color='darkgray', label=diversitymetric, alpha=1)
		axis.axhline(y=0, color='k', linestyle=':')
		axis.spines['top'].set_visible(False)
		axis.spines['right'].set_visible(False)
		axis.set_ylim(0,1.5)
	plt.suptitle("landscape trajectory: %s\n%s\n%s\n%s" % (landscape_trajectory,starting_pop,noisepergen,noisepersim), fontsize=10)

	if True: # args.legendtype == 'single':
		labels_handles = {
	 	 label: handle for ax in fig.axes for handle, label in zip(*ax.get_legend_handles_labels())
		}
		leg = fig.legend(
		  labels_handles.values(),
		  labels_handles.keys(),
		  bbox_to_anchor=(1,1),
		  loc='upper right',
		  fontsize=6,
		  #bbox_transform = plt.gcf().transFigure,
		)
		for legobj in leg.legendHandles:
			legobj.set_linewidth(2.0)

	plt.tight_layout()
	plt.savefig(args.simulationdir+'/pdfs/%s.pdf' % diversitymetric)
	plt.savefig(args.simulationdir+'/%s.png' % diversitymetric)
	plt.close()














