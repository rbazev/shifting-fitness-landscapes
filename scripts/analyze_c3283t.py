import myutils
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='plot correlations in genotypes contining c3283t mutation')
parser.add_argument('srfilesdir', type=str, help='path to sr files')
parser.add_argument('plotdir', type=str, help='path to save output plots')
args = parser.parse_args()

# scatter selection rate of genotypes without and with the c3283t mutation

samples = [
			'LibA_vGen_1_unprogrammed',
			'LibA_vGen_2_unprogrammed',
			'LibA_vGen_3_unprogrammed',
			'LibA_vGen_LamB_1_unprogrammed',
			'LibA_vGen_LamB_2_unprogrammed',
			'LibA_vGen_LamB_3_unprogrammed',
			'LibA_vGen_OmpF_1_unprogrammed',
			'LibA_vGen_OmpF_2_unprogrammed',
			'LibA_vGen_OmpF_3_unprogrammed',
			'LibA_vLspec_1_unprogrammed',
			'LibA_vLspec_2_unprogrammed',
			'LibA_vLspec_3_unprogrammed',
			'LibA_vOspec_1_unprogrammed',
			'LibA_vOspec_2_unprogrammed',
			'LibA_vOspec_3_unprogrammed',
			'LibB_vGen_1_unprogrammed',
			'LibB_vGen_2_unprogrammed',
			'LibB_vGen_3_unprogrammed',
			'LibB_vGen_LamB_1_unprogrammed',
			'LibB_vGen_LamB_2_unprogrammed',
			'LibB_vGen_LamB_3_unprogrammed',
			'LibB_vGen_OmpF_1_unprogrammed',
			'LibB_vGen_OmpF_2_unprogrammed',
			'LibB_vGen_OmpF_3_unprogrammed',
			'LibB_vLspec_1_unprogrammed',
			'LibB_vLspec_2_unprogrammed',
			'LibB_vLspec_3_unprogrammed',
			'LibB_vOspec_1_unprogrammed',
			'LibB_vOspec_2_unprogrammed',
			'LibB_vOspec_3_unprogrammed',
			 ]

srfiles = ["%s/sr_%s.csv" % (args.srfilesdir, s) for s in samples]

plotnan = -10
nrows = 6
ncols = 5
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4*ncols,4*nrows))
if len(samples) == 1:
		axes = [axes]
else:
	axes = axes.ravel() # axs.ravel() converts 2-dim axes object to a list in row-major order

for i,f in enumerate(srfiles):
	sr_dict = myutils.read_selectionrate_file(f, as_dict=True)
	axis = axes[i]
	title = f[f.find("Lib"):-4]
	xs = [] # sr of gt without LF
	ys = [] # sr of gt with LF

	for gt in myutils.construct_genotype_space():
		if gt in sr_dict and np.isfinite(sr_dict[gt]):
			xs.append(sr_dict[gt])
		else:
			xs.append(plotnan)

		gt_LF = gt + '+183A' # this is the c3283t mutation in "longset" amplicon setting in parse_FASTQ.py. I call it "LF" in this code, it's a Leu>Phe
		if gt_LF in sr_dict and np.isfinite(sr_dict[gt_LF]):
			ys.append(sr_dict[gt_LF])
		else:
			ys.append(plotnan)

	axis.scatter(xs,ys,color='k',s=14,linewidth=0.107*3.5, edgecolor='k', alpha=0.6)
	axis.set_xlabel("without c3283t", fontsize=8)
	axis.set_ylabel("with c3283t", fontsize=8)
	axis.set_title(title, fontsize=8)
	axis.spines[['top', 'right']].set_visible(False)
	axis.set_xlim(-11,10)
	axis.set_ylim(-11,10)
	axis.axline((0, 0), slope=1, color = "grey", linestyle = ":", linewidth = 1)

plt.tight_layout()
plt.savefig('%s/c3283t.pdf' % args.plotdir)
plt.close()




