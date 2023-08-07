import myutils
import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='plot correlations in genotypes contining G2, G3, G2,3')
parser.add_argument('srfilesdir', type=str, help='path to sr files')
parser.add_argument('plotdir', type=str, help='path to save output plots')
args = parser.parse_args()

def convert_g2_to_g3(gt):
	assert gt[2] == '1' and gt[3] == '0'
	return gt[:2] + '01' + gt [4:]

def convert_g2_to_g23(gt):
	assert gt[2] == '1' and gt[3] == '0'
	return gt[:2] + '11' + gt [4:]

def convert_g3_to_g23(gt):
	assert gt[3] == '1' and gt[2] == '0'
	return gt[:2] + '11' + gt [4:]
 
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

srfiles = ["%s/sr_%s.csv" % (args.srfilesdir, s) for s in samples]


# plot G2 vs. G3:
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
	# iterate over genotypes G2
	# scatter G2 sr vs. G3 sr
	xs = [] # fitness w G2
	ys = [] # fitness w G3
	for gt in myutils.construct_genotype_space():
		if gt[2] == '1' and gt[3] == '0':
			if gt in sr_dict and np.isfinite(sr_dict[gt]):
				xs.append(sr_dict[gt])
			else:
				xs.append(plotnan)
			newgt = convert_g2_to_g3(gt)
			if newgt in sr_dict and np.isfinite(sr_dict[newgt]):
				ys.append(sr_dict[newgt])
			else:
				ys.append(plotnan)
	axis.scatter(xs,ys,color='k',s=14,linewidth=0.107*3.5, edgecolor='k', alpha=0.6)
	axis.set_xlabel("G2=1, G3=0", fontsize=8)
	axis.set_ylabel("G2=0, G3=1", fontsize=8)
	axis.set_title(title, fontsize=8)
	axis.spines[['top', 'right']].set_visible(False)
	axis.set_xlim(-11,10)
	axis.set_ylim(-11,10)
	axis.axline((0, 0), slope=1, color = "grey", linestyle = ":", linewidth = 1)
plt.tight_layout()
plt.savefig('%s/g2g3.pdf' % args.plotdir)
plt.savefig('%s/g2g3.png' % args.plotdir)
plt.close()


# plot G2 vs. G23:
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
	# iterate over genotypes G2
	# scatter G2 sr vs. G3 sr
	xs = [] # fitness w G2
	ys = [] # fitness w G3
	for gt in myutils.construct_genotype_space():
		if gt[2] == '1' and gt[3] == '0':
			if gt in sr_dict and np.isfinite(sr_dict[gt]):
				xs.append(sr_dict[gt])
			else:
				xs.append(plotnan)
			newgt = convert_g2_to_g23(gt)
			if newgt in sr_dict and np.isfinite(sr_dict[newgt]):
				ys.append(sr_dict[newgt])
			else:
				ys.append(plotnan)
	axis.scatter(xs,ys,color='k',s=14,linewidth=0.107*3.5, edgecolor='k', alpha=0.6)
	axis.set_xlabel("G2=1, G3=0", fontsize=8)
	axis.set_ylabel("G2=1, G3=1", fontsize=8)
	axis.set_title(title, fontsize=8)
	axis.spines[['top', 'right']].set_visible(False)
	axis.set_xlim(-11,10)
	axis.set_ylim(-11,10)
	axis.axline((0, 0), slope=1, color = "grey", linestyle = ":", linewidth = 1)
plt.tight_layout()
plt.savefig('%s/g2g23.pdf' % args.plotdir)
plt.savefig('%s/g2g23.png' % args.plotdir)
plt.close()



