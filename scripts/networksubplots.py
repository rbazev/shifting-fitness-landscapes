import myutils
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import argparse

parser = argparse.ArgumentParser(description='plot network diagram subplots for all samples')
parser.add_argument('sr_files_dir', type=str, help='path to sr files')
parser.add_argument('plot_dir', type=str, help='path to save network plots')
parser.add_argument('runname', type=str, help='runname prefix')

# if coloring, specify both of these arguments. If multiple colordatafiles given, colortype must be categorical, and this script will not check to ensure that the sets are mutually exclusive:
parser.add_argument('--colortype', type=str, help='categorical or quantitative')
parser.add_argument('--colordatafiles', type=str, nargs='+', help='path to file(s) to use for quantitative or categorical coloring')
parser.add_argument('--colordatalabels', type=str, nargs='+', help='labels associated with the given colordatafiles for categorical coloring')

args = parser.parse_args()
myutils.makedir(args.plot_dir)
subplots_output_dir = args.plot_dir + '/%s/' % args.runname
myutils.makedir(subplots_output_dir) # subdirectory for this runname.

if args.colordatafiles:
	if len(args.colordatafiles) > 1:
		assert args.colortype == 'categorical'
		colortype = 'categorical'
	else:
		colortype = args.colortype
	colordata_dicts = []
	for file in args.colordatafiles:
		colordata_dicts.append(myutils.read_selectionrate_file(file,as_dict=True))
else:
	colordata_dicts = []
	colortype = 'categorical'

# used to either color quantitatively (eg, by values listed in the csv using cmap)
# or categorically (eg. all genotypes in the given file will be given the same color)
# this same colordata will be used for all of the subplots figures created:
#colortype = 'categorical' # v. 'quantitative'. categorig
default_color = 'lightgray' 
quant_cmap = plt.cm.coolwarm.reversed() # used if colortype == 'quantitative'; reversed_map = orig_map.reversed()
(vmin, vmax) = (-1,1)
categorical_colors = ['dodgerblue', 'tab:orange', 'lime']  # used if colortype == 'categorical'
nodesize = 13
nodelinewidths=0.45
jitter = 0.22
edgewidth = 0.005
nodealpha = 1
edgealpha = 0.2
titlefontsize = 8


def plotNetworkSubplots(samples, nrows, ncols, fileprefix, squish=False):
	assert ncols*nrows == len(samples)
	if squish:
		figsize = (9.25,3.25)
	else:
		figsize = (4*ncols,4*nrows)
	fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)

	if len(samples) == 1:
		axes = [axes]
	else:
		axes = axes.ravel() # axs.ravel() converts 2-dim axes object to a list in row-major order

	for axis, sample in zip(axes, samples): 
		
		srdata = myutils.read_selectionrate_file('%s/sr_%s.csv' % (args.sr_files_dir,sample), as_dict=True)

		# is this 8-cube or 9-cube?
		length = 9 if '[000000000]' in srdata else 8 # this assumes that it will only ever see length 8 or 9

		# develop G, networkx graph object containing genotypes observed in this sample's sr file
		G = myutils.BuildHypercubeGraph(length=length)
		for gt in myutils.construct_genotype_space(length=length):
			if np.isnan(srdata[gt]):
				G.remove_node(gt)
		# develop pos, list of (x,y) positions in node order
		pos = myutils.position_bysr(G,srdata,jitter)
		# develop c, list of colors (or values to colormap)
		if colortype == 'quantitative':
			c = [colordata_dicts[0][node] for node in G] # c = [colordata[node] for node in G]
			cmap = quant_cmap
		elif colortype == 'categorical':
			c = []
			for gt in G.nodes:
				found_gt = False
				for colorset_i in range(len(colordata_dicts)):
					if gt in colordata_dicts[colorset_i]: # note that nothing will be colored if an 8-d colorfile given for 9-d graph, and vice-versa
						c.append(categorical_colors[colorset_i])
						found_gt = True
						break
				if not found_gt:
					c.append(default_color)
			cmap = False

		axis.axhline(color = "grey", linestyle = ":", linewidth = 1)
		nx.draw_networkx_edges(G=G, pos = pos, ax=axis, edge_color='k', alpha=edgealpha, width=edgewidth)
		nx.draw_networkx_nodes(G=G, pos = pos, ax=axis, alpha=nodealpha, linewidths=nodelinewidths, edgecolors='k',
                        node_size = nodesize, node_color=c, cmap=cmap, vmin=vmin, vmax=vmax)

		axis.set_title('\n%s' % sample, fontsize=titlefontsize)

		axis.set_xlim(-0.25,length+0.25)
		axis.set_ylim(-10,10)
		plt.axis('on')
		axis.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True, labelsize=8)

		if len(samples) == 1:
			axis.set_xlabel("number of mutations")
			axis.set_ylabel("selection rate")
		axis.set_xticks([i for i in range(length+1)])

		axis.spines['top'].set_visible(False)
		axis.spines['right'].set_visible(False)

	if len(samples) > 1 and not squish:
		fig.supxlabel('Hamming distance from generalist', fontsize=8)
		fig.supylabel('selection rate', fontsize=8)

	if squish:
		for i in range(len(axes)):
			if i == 0:
				axes[i].set_ylabel("selection rate", fontsize=8)
			if i == 2:
				axes[i].set_xlabel("hamming distance", fontsize=8)
			if i>0: 
				axes[i].get_yaxis().set_visible(False)

	plt.tight_layout()
	plt.savefig('%s/%s_subplt-by-%s.pdf' % (subplots_output_dir, fileprefix, args.runname))
	plt.savefig('%s/%s_subplt-by-%s.png' % (subplots_output_dir, fileprefix, args.runname))

	plt.close()




# harcode all the arrangements of which samples are in which orders on which plots

# ALL INDIVIDUAL SAMPLES
samples = [
			'LibA_vGen_1',
			'LibA_vGen_2',
			'LibA_vGen_3',
			'LibA_vGen_LamB_1',
			'LibA_vGen_LamB_2',
			'LibA_vGen_LamB_3',
			'LibA_vGen_OmpF_1',
			'LibA_vGen_OmpF_2',
			'LibA_vGen_OmpF_3',
			'LibA_vLspec_1',
			'LibA_vLspec_2',
			'LibA_vLspec_3',
			'LibA_vOspec_1',
			'LibA_vOspec_2',
			'LibA_vOspec_3',
			'LibB_vGen_1',
			'LibB_vGen_2',
			'LibB_vGen_3',
			'LibB_vGen_LamB_1',
			'LibB_vGen_LamB_2',
			'LibB_vGen_LamB_3',
			'LibB_vGen_OmpF_1',
			'LibB_vGen_OmpF_2',
			'LibB_vGen_OmpF_3',
			'LibB_vLspec_1',
			'LibB_vLspec_2',
			'LibB_vLspec_3',
			'LibB_vOspec_1',
			'LibB_vOspec_2',
			'LibB_vOspec_3',
			 ]
for s in samples:
	(ncols,nrows) = (1,1)
	fileprefix = s

	if True:	plotNetworkSubplots([s], nrows, ncols, fileprefix)

# now a plot with all 30 competition flasks;
# for each 5 conditions plot the 6 replicates (3 w/ a, 3 w/ b)
vgensamps = [s for s in samples if 'vGen' in s and 'LamB' not in s and 'OmpF' not in s]
vgen_lambsamps = [s for s in samples if 'vGen_LamB' in s]
vgen_ompfsamps = [s for s in samples if 'vGen_OmpF' in s]
vlspecsamps = [s for s in samples if 'vLspec' in s]
vospecsamps = [s for s in samples if 'vOspec' in s]

if True:
	plotNetworkSubplots(vgensamps, 1, 6, 'vgen_reps')
	plotNetworkSubplots(vgen_lambsamps, 1, 6, 'vgenLamB_reps')
	plotNetworkSubplots(vgen_ompfsamps, 1, 6, 'vgenOmpF_reps')
	plotNetworkSubplots(vlspecsamps, 1, 6, 'vLspec_reps')
	plotNetworkSubplots(vospecsamps, 1, 6, 'vOspec_reps')


# now the fully averaged for all 5 environments:
averaged_samples = ['LibAvg_vGen_OmpF',
					'LibAvg_vLspec',
					'LibAvg_vGen',
					'LibAvg_vOspec',
					'LibAvg_vGen_LamB']
plotNetworkSubplots(averaged_samples, 1, 5, 'five_environments')
plotNetworkSubplots(averaged_samples, 1, 5, 'five_environments_squish', squish=True)