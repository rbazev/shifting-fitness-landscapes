import argparse
import myutils
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr

parser = argparse.ArgumentParser(description='plot selection rate correlation')
parser.add_argument('srfile1', type=str, help='path to first sr file')
parser.add_argument('srfile2', type=str, help='path to second sr file')
parser.add_argument('label1', type=str, help='label for first sample')
parser.add_argument('label2', type=str, help='label for second sample')
parser.add_argument('figpath', type=str, help='path for output figures')
parser.add_argument('--plotnan', type=float, help='plot NaN (unobserved genotypes) as negative this value (eg 10-->-10)')
args = parser.parse_args()

myutils.makedir(args.figpath)

#print('\ndrawing scatter for %s vs %s;\n%s vs %s.\n' % (args.srfile1, args.srfile2, args.label1, args.label2))
#if args.plotnan: print('plotting NaN values at %s' % args.plotnan)

s1_data_dict = myutils.read_selectionrate_file(args.srfile1, as_dict=True)
s2_data_dict = myutils.read_selectionrate_file(args.srfile2, as_dict=True)

xs = [] # for plotting (can include genotypes only observed in one experiment to be plotted but not included in correlations)
ys = []
colors = []
marked_gt = []

xs_for_correlation = [] # only include those genotypes observed in both experiments.
ys_for_correlation = []

gts = myutils.construct_genotype_space()

for gt in gts:
	if np.isfinite(s1_data_dict[gt]) and np.isfinite(s2_data_dict[gt]):
		# this genotype is observed in both conditions
		xs.append(s1_data_dict[gt])
		ys.append(s2_data_dict[gt])
		xs_for_correlation.append(s1_data_dict[gt])
		ys_for_correlation.append(s2_data_dict[gt])
		colors.append(myutils.encode_genotype_color(gt, schema="black"))

	else:
		# this genotype is unobserved in at least one of the two conditions
		if args.plotnan:
			xs.append(s1_data_dict[gt] if np.isfinite(s1_data_dict[gt]) else -args.plotnan)
			ys.append(s2_data_dict[gt] if np.isfinite(s2_data_dict[gt]) else -args.plotnan)
			# now find outliers that have sr > X in observed condition, 
			# and NaN (unobserved presumably due to intense negative selection) in the other condition
			if np.isfinite(s1_data_dict[gt]):
				assert np.isnan(s2_data_dict[gt])
				colors.append('gray')
			else:
				assert np.isnan(s1_data_dict[gt])
				colors.append('gray')
	
fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(3,3))
axs.scatter(xs,ys,color=colors,s=14,linewidth=0.107*3.5, edgecolor='k', alpha=0.7)
axs.set_xlabel(args.label1, fontsize=8)
axs.set_ylabel(args.label2, fontsize=8)

if args.plotnan:
	axs.set_xlim(-11,10)
	axs.set_ylim(-11,10)
else:
	axs.set_xlim(-11,10)
	axs.set_ylim(-11,10)

plt.xticks(fontsize=8)
plt.yticks(fontsize=8)

axs.axline((0, 0), slope=1, color = "grey", linestyle = ":", linewidth = 1)

result = pearsonr(xs_for_correlation,ys_for_correlation)
axs.text(-10,8,'r=%.2f'%result.statistic,fontsize=10)
axs.spines[['top', 'right']].set_visible(False)
axs.set_title("selection rate", fontsize=8)
plt.tight_layout()
plt.savefig(args.figpath+'/%s-vs-%s.pdf' % (args.label1, args.label2))
plt.savefig(args.figpath+'/%s-vs-%s.png' % (args.label1, args.label2))
plt.close()