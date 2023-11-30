import argparse
import myutils
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='plot distribution of fitness effects for each fitness landscape')
parser.add_argument('srfilesdir', type=str, help='directory containing the selection rate files')
parser.add_argument('plotsavedir', type=str, help='directory to save plot')
args = parser.parse_args()

landscape_filenames = ['sr_LibAvg_vGen_OmpF.csv', 'sr_LibAvg_vLspec.csv', 'sr_LibAvg_vGen.csv', 'sr_LibAvg_vOspec.csv', 'sr_LibAvg_vGen_LamB.csv']

# plot histrogram of fitness values for each fitness landscape
fig, axes = plt.subplots(nrows=1, ncols=5, figsize=(9.25,3.25))
axes = axes.ravel() # axs.ravel() converts 2-dim axes object to a list in row-major order


for (axis,filename) in zip (axes, landscape_filenames):
	fitness_dict = myutils.read_selectionrate_file('%s/%s' % (args.srfilesdir,filename), as_dict=True)
	fitness_list = list(fitness_dict.values())
	axis.hist(fitness_list, density=True, bins=40, alpha=0.55, facecolor='k', edgecolor='k', linewidth=0.5, label=filename)
	axis.axvline(x=0, ls=':', color='k')

	# standardize axis limits across subplots
	axis.set_ylim(0,.4)

	axis.set_title(filename, fontsize=8)

	axis.spines['top'].set_visible(False)
	axis.spines['right'].set_visible(False)

# aesthetics:
for i in range(len(axes)):
	if i == 0:
		axes[i].set_ylabel("Density", fontsize=8)
	if i == 2:
		axes[i].set_xlabel("Selection rate", fontsize=8)
	if i>0: 
		axes[i].get_yaxis().set_visible(False)

plt.tight_layout()
plt.savefig('%s/DFE.pdf' % (args.plotsavedir))
plt.savefig('%s/DFE.png' % (args.plotsavedir))

# debug: @ Jlib-pilot-analysis/SFL-revisionsandbox: python ./scripts/plot_DFE.py ./output-notebook/selectionrate/srfiles/ ./output-notebook/selectionrate/

