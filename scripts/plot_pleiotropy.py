import argparse
import myutils
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='analyze pleiotropic effects of mutations on fitness on each receptor')
parser.add_argument('srfilesdir', type=str, help='directory containing the selection rate files')
parser.add_argument('plotsavedir', type=str, help='directory to save plot')

args = parser.parse_args()

sr_on_L_file = '%s/sr_LibAvg_vGen_LamB.csv' % args.srfilesdir
sr_on_O_file = '%s/sr_LibAvg_vGen_OmpF.csv' % args.srfilesdir
fitness_on_LamB = myutils.read_selectionrate_file(sr_on_L_file, as_dict=True)
fitness_on_OmpF = myutils.read_selectionrate_file(sr_on_O_file, as_dict=True)

# for each possible mutation (0 --> 1) there is a change in fitness on LamB,
#  and a change in fitness on OmpF.
# here we will iterate over all mutations (0->1) and keep track of
# the effect that mutation has on change in fitness on LamB (delta on lamB)
# and the effect that mutation has on change in fitness on OmpF (delta on OmpF).

deltas_on_LamB = []
deltas_on_OmpF = []
num_posL_posO = 0 # to count how many mutations fall in each of the four quadrants
num_posL_negO = 0
num_negL_posO = 0
num_negL_negO = 0

# # optionally, for visualization purposes, assign this value to unobserved genotypes as selection rate:
nanval = -10 
# however the choice of value is arbitrary and i think it's more confusing to perform this analysis on 
# genotypes that have a fitness so low they are unobserved, so the more straightforward way to analyze is
# to exclude any genotypes that have an unobserved fitness on either receptor, so that we restrict our
# analysis to strictly observed fitness measurements.
exclude_nan  = True
# we should keep track of how many genotype pairs are rejected from the analysis due to at least one 
# unobserved fitness:
counts = {	"genotype_pairs_total":0, # total examined pairs of genotypes separated by any 0>1 at the 9 sites.
			"rejected_nan_L":0, # pairs of genotypes rejected for unobserved fitness on LamB
			"rejected_nan_O":0, # pairs of genotypes rejected for unobserved fitness on LamB 

			}

for gt in myutils.construct_genotype_space(length = 9):
	for i,pos in enumerate(gt[1:10]): # examine all single mutations from this starting genotype
		if pos == '0':
			newgt = gt[:i+1] + '1'+ gt[i+2:]
			counts["genotype_pairs_total"] += 1

			# replace unobserved selection rates with `nanval` for visualization purposes, or exclude altogether if exclude_nan == True
			if np.isnan(fitness_on_LamB[gt]):
				if exclude_nan:
					counts['rejected_nan_L'] += 1
					continue
				else:
					fitness_on_LamB[gt] = nanval
			if np.isnan(fitness_on_LamB[newgt]):
				if exclude_nan:
					counts['rejected_nan_L'] += 1
					continue
				else:
					fitness_on_LamB[newgt] = nanval
			if np.isnan(fitness_on_OmpF[gt]):
				if exclude_nan:
					counts['rejected_nan_O'] += 1
					continue
				else:
					fitness_on_OmpF[gt] = nanval
			if np.isnan(fitness_on_OmpF[newgt]):
				if exclude_nan:
					counts['rejected_nan_O'] += 1
					continue
				else:
					fitness_on_OmpF[newgt] = nanval

			# compute deltas
			delta_on_LamB = fitness_on_LamB[newgt] - fitness_on_LamB[gt]
			delta_on_OmpF = fitness_on_OmpF[newgt] - fitness_on_OmpF[gt]
			deltas_on_LamB.append(delta_on_LamB)
			deltas_on_OmpF.append(delta_on_OmpF)

			if delta_on_LamB >= 0:
				if delta_on_OmpF >= 0: num_posL_posO += 1 
				else: num_posL_negO += 1
			else:
				if delta_on_OmpF >= 0: num_negL_posO += 1
				else: num_negL_negO += 1

print("Excludenan = %s" % exclude_nan)
print(counts)
num_observed = float(num_posL_posO + num_posL_negO + num_negL_posO + num_negL_negO)
print("num obs = %s" % num_observed)

#===========
# thank you https://matplotlib.org/stable/gallery/lines_bars_and_markers/scatter_hist.html

# Start with a square Figure.
fig = plt.figure(figsize=(3, 3))
# Add a gridspec with two rows and two columns and a ratio of 1 to 4 between
# the size of the marginal axes and the main axes in both directions.
# Also adjust the subplot parameters for a square plot.
gs = fig.add_gridspec(2, 2,  width_ratios=(4, 1), height_ratios=(1, 4),
                      left=0.1, right=0.9, bottom=0.1, top=0.9,
                      wspace=0.05, hspace=0.05)
# Create the Axes.
ax = fig.add_subplot(gs[1, 0])
ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)

# no labels
ax_histx.tick_params(axis="x", labelbottom=False)
ax_histy.tick_params(axis="y", labelleft=False)

# the scatter plot:
ax.scatter(deltas_on_LamB,deltas_on_OmpF,color='dimgray',s=8,linewidth=0.107, edgecolor='k', alpha=0.5)
ax.set_xlim(-10,10)
ax.set_ylim(-20,20)
ax.set_xlabel("Change in selection rate on LamB", fontsize=8)
ax.set_ylabel("Change in selection rate on OmpF", fontsize=8)
ax.spines[['top', 'right']].set_visible(False)
ax.set_xticks([-10,0,10])
ax.set_yticks([-10,0,10])

ax.axhline(y=0, color = "k", linestyle = ":", linewidth = 1.07)
ax.axvline(x=0, color = "k", linestyle = ":", linewidth = 1.07)

ax.text(8,8,"%.0f%%" % (num_posL_posO/num_observed*100), fontsize=8)
ax.text(8,-8,"%.0f%%" % (num_posL_negO/num_observed*100), fontsize=8)
ax.text(-9,-8,"%.0f%%" % (num_negL_negO/num_observed*100), fontsize=8)
ax.text(-9,8,"%.0f%%" % (num_negL_posO/num_observed*100), fontsize=8)

# now determine nice limits by hand:
binwidth = 0.2
xymax = max(np.max(np.abs(deltas_on_LamB)), np.max(np.abs(deltas_on_OmpF)))
lim = (int(xymax/binwidth) + 1) * binwidth
bins = np.arange(-lim, lim + binwidth, binwidth)
ax_histx.hist(deltas_on_LamB, bins=bins, color='dimgray', rwidth=0.8)
ax_histy.hist(deltas_on_OmpF, bins=bins, color='dimgray', orientation='horizontal', rwidth=0.8)
ax_histx.spines[['top', 'right']].set_visible(False)
ax_histy.spines[['top', 'right']].set_visible(False)

ticksize = 8
for tick in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks():
	tick.label.set_fontsize(ticksize)

ax_histx.spines[['top', 'right', 'left']].set_visible(False)
ax_histx.get_yaxis().set_ticks([])
ax_histy.spines[['top', 'right', 'bottom']].set_visible(False)
ax_histy.get_xaxis().set_ticks([])

gs.tight_layout(fig)
plt.savefig('%s/pleiotropy.pdf' % (args.plotsavedir))
plt.savefig('%s/pleiotropy.png' % (args.plotsavedir))
