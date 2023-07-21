import argparse
import myutils
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr

parser = argparse.ArgumentParser(description='plot selection rate correlations across replicate experiments')
parser.add_argument('srfilesdir', type=str, help='path to directory containing sr files')
parser.add_argument('figpath', type=str, help='path for output figures')
parser.add_argument('--plotnan', type=float, default=10, help='plot NaN (unobserved genotypes) as the negative of this value')
args = parser.parse_args()
myutils.makedir(args.figpath)

plotnan = -args.plotnan
nrows = 5 # 5 rows for the five experimental conditions
ncols = 7 # 7 columns showing 7 different correlation plots for each condition
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4*ncols,4*nrows))

row_designations = {0:'vGen_OmpF',
					1:'vLspec',
					2:'vGen',
					3:'vOspec',
					4:'vGen_LamB'}

for irow in range(0,nrows):

	environment = row_designations[irow]

	for icol in range(0,ncols):
		
		# define the two sr files used for this axis' plot
		if icol == 0:
			# LibA 1 v 2:
			f1 = "sr_LibA_%s_1.csv" % environment 
			f2 = "sr_LibA_%s_2.csv" % environment 
		elif icol == 1:
			# LibA 1 v 3:
			f1 = "sr_LibA_%s_1.csv" % environment 
			f2 = "sr_LibA_%s_3.csv" % environment 
		elif icol == 2:
			# LibA 2 v 3:
			f1 = "sr_LibA_%s_2.csv" % environment 
			f2 = "sr_LibA_%s_3.csv" % environment 
		elif icol == 3:
			# LibB 1 v 2:
			f1 = "sr_LibB_%s_1.csv" % environment 
			f2 = "sr_LibB_%s_2.csv" % environment 
		elif icol == 4:
			# LibB 1 v 3:
			f1 = "sr_LibB_%s_1.csv" % environment 
			f2 = "sr_LibB_%s_3.csv" % environment 
		elif icol == 5:
			# LibB 2 v 3:
			f1 = "sr_LibB_%s_2.csv" % environment 
			f2 = "sr_LibB_%s_3.csv" % environment 
		elif icol == 6:
			# LibA avg v LibB avg:
			f1 = "sr_LibA_%s_avg.csv" % environment 
			f2 = "sr_LibB_%s_avg.csv" % environment 

		s1_data_dict = myutils.read_selectionrate_file(args.srfilesdir + '/%s' % f1, as_dict=True)
		s2_data_dict = myutils.read_selectionrate_file(args.srfilesdir + '/%s' % f2, as_dict=True)

		xs_plot = [] # for plotting coordinates; unobserved genotypes are plotted at -plotnan
		ys_plot = []
		xs_corr = [] # for calculating correlation coefficient; only included if observed on both files.
		ys_corr = []

		for gt in myutils.construct_genotype_space():

			val1 = s1_data_dict[gt] if np.isfinite(s1_data_dict[gt]) else plotnan
			val2 = s2_data_dict[gt] if np.isfinite(s2_data_dict[gt]) else plotnan
			xs_plot.append(val1)
			ys_plot.append(val2)

			if np.isfinite(s1_data_dict[gt]) and np.isfinite(s2_data_dict[gt]):
				xs_corr.append(val1)
				ys_corr.append(val2)

		axes[irow,icol].scatter(xs_plot,ys_plot,color='k',s=14,linewidth=0.107*3.5, edgecolor='k', alpha=0.7)
		axes[irow,icol].set_xlabel(f1, fontsize=16)
		axes[irow,icol].set_ylabel(f2, fontsize=16)
		axes[irow,icol].set_xlim(-11,10)
		axes[irow,icol].set_ylim(-11,10)

		axes[irow,icol].xaxis.set_tick_params(labelsize=14)
		axes[irow,icol].yaxis.set_tick_params(labelsize=14)

		axes[irow,icol].axline((0, 0), slope=1, color = "grey", linestyle = ":", linewidth = 1)

		if len(xs_corr) > 2:
			result = pearsonr(xs_corr,ys_corr)
			axes[irow,icol].text(-10,8,'r=%.2f'%result.statistic,fontsize=28)

		axes[irow,icol].spines[['top', 'right']].set_visible(False)

plt.tight_layout()
plt.savefig(args.figpath+'/sr-correlations.pdf')
plt.savefig(args.figpath+'/sr-correlations.png')


# remove ticks/axis labels, maybe add row labels in inkscape:
for irow in range(0,nrows):
	environment = row_designations[irow]
	for icol in range(0,ncols):
		axes[irow,icol].set_xticks([])
		axes[irow,icol].set_yticks([])
		axes[irow,icol].xaxis.label.set_visible(False)
		axes[irow,icol].yaxis.label.set_visible(False)
plt.tight_layout()
plt.savefig(args.figpath+'/sr-correlations-clean.pdf')
plt.savefig(args.figpath+'/sr-correlations-clean.png')

plt.close()