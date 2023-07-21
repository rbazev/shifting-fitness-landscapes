import argparse
import myutils
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

'''
This script examines the selection rate data to identify genotypes
with large differences in fitness between conditions on L-only host (OmpF knockout)
and O-only host (LamB knockoyt). I specify that genotypes are a specialist on a
receptor if they have much better fitness on that receptor than the other receptor. 

A specialization index is computed for every genotype; this continuous variable is useful
in subsequent analyses (eg simulated evolution). The most extreme positive and negative 
specialization index genotypes are (with arbitrary cutoff) categorized as L- or O- specialists.
'''

parser = argparse.ArgumentParser(description='parse specialist genotypes from single-host competition data')
parser.add_argument('selectionrateonL', type=str, help='path to the selection rate file for LamB')
parser.add_argument('selectionrateonO', type=str, help='path to the selection rate file for OmpF')
parser.add_argument('specialistdir', type=str, help='output directory to save specialist genotype files/plots')
args = parser.parse_args()
myutils.makedir(args.specialistdir)

print('\nIdentifying specialist genotypes using arguments:\n%s\n' % (args))
print('Using selection rate data from files:')
print('Sr on L:%s' % args.selectionrateonL)
print('Sr on O:%s' % args.selectionrateonO)
print('Saving results to %s' % args.specialistdir)

si_fileout = open(args.specialistdir + '/si.csv', 'w')
si_fileout.write('genotype,si\n')
Lspec_fileout = open(args.specialistdir + '/Lspecialists.csv', 'w')
Lspec_fileout.write('genotype,Lspecialist_si\n')
Ospec_fileout = open(args.specialistdir + '/Ospecialists.csv', 'w')
Ospec_fileout.write('genotype,Ospecialist_si\n')
sr_on_L = myutils.read_selectionrate_file(args.selectionrateonL, as_dict=True)
sr_on_O = myutils.read_selectionrate_file(args.selectionrateonO, as_dict=True)

genotypes = myutils.construct_genotype_space()

########## Compute a specialization index for every genotype:
# the specialization index Si is OmpF-centric; +1 means fully O specialist and -1 means fully L specialist
si_list = [] # in gt order.
si_dict = {}
onanval = -99 # for purposes of calculating Si (although for plotting, to see them on the plot, a different value used below)
lnanval = -99
for gt in genotypes:
	if np.isfinite(sr_on_L[gt]):
		srl = sr_on_L[gt]
	else:
		srl = lnanval
	if np.isfinite(sr_on_O[gt]):
		sro = sr_on_O[gt]
	else:
		sro = onanval
	si = (np.exp(sro)-np.exp(srl))/(np.exp(sro)+np.exp(srl))
	si_list.append(si)
	si_dict[gt] = si

	#output si file:
	si_fileout.write('%s,%s\n' % (gt,si_dict[gt]))

# make histogram of si
fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(3,3))
axs.hist(si_list,30)
axs.set_xlim(-1,1)
axs.spines[['top', 'right']].set_visible(False)
axs.set_title("specialization index distribution", fontsize=8)
plt.tight_layout()
plt.savefig(args.specialistdir + '/specialization_index_dist.pdf')
plt.close()

# write specialist files with si > threshold = O-specialists and < -threshold = L-specialists
threshold = 0.33
for gt in myutils.construct_genotype_space():
	if si_dict[gt] > threshold:
		# O specialist:
		Ospec_fileout.write(gt+',%.7f\n' % si_dict[gt])
	elif si_dict[gt] < -threshold:
		# L specialist:
		Lspec_fileout.write(gt+',%.7f\n' % si_dict[gt])

Ospec_fileout.close()
Lspec_fileout.close()

# Will build up lists for plotting as we identify the specialists.
xs = []
ys = []

# aesthetic parameters:
plotnanval = -9

for gt in genotypes:
	if np.isfinite(sr_on_L[gt]) and np.isfinite(sr_on_O[gt]):
		xs.append(sr_on_L[gt])
		ys.append(sr_on_O[gt])
	elif np.isfinite(sr_on_L[gt]):
		# this genotype only observed on L receptor
		xs.append(sr_on_L[gt])
		ys.append(plotnanval)
	elif np.isfinite(sr_on_O[gt]):
		# this genotype only observed on O receptor
		xs.append(plotnanval)
		ys.append(sr_on_O[gt])
	else:
		xs.append(plotnanval)
		ys.append(plotnanval)

fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(3.25,3))
axs.axline((0, 0), slope=1, color = "grey", linestyle = ":", linewidth = 1.07/1.5)
sc = axs.scatter(xs,ys,c=si_list, cmap='coolwarm',vmin=-1,vmax=1,s=8,linewidth=0.107, edgecolor='k', alpha=0.9) # 'PuOr_r'

# overlay points specifically for the three competitors used in the experiment:
Lspec = ('[100000010]', 'b') 
Ospec = ('[101001011]', 'r')
Genst = ('[000000000]', '#888888')

for  (gt,c) in [Lspec,Ospec, Genst]:
	si_list = []
	if np.isfinite(sr_on_L[gt]) and np.isfinite(sr_on_O[gt]):
		x = sr_on_L[gt]
		y = sr_on_O[gt]
	elif np.isfinite(sr_on_L[gt]):
		# this genotype not observed on Omp-F-only hosts
		x = sr_on_L[gt]
		y = plotnanval
	elif np.isfinite(sr_on_O[gt]):
		# this genotype not observed on LamB-only hosts
		x = plotnanval
		y = sr_on_O[gt]
	si_list.append(si_dict[gt])
	axs.scatter(x,y,c=si_list,cmap='coolwarm',vmin=-1,vmax=1,s=16,linewidth=0.107*11, edgecolor='k', alpha=1.0 ,marker='o')

axs.set_xlabel("Selection rate on LamB", fontsize=8)
axs.set_ylabel("Selection rate on OmpF", fontsize=8)
axs.set_xlim(-10,4)
axs.set_ylim(-10,10)
plt.yticks(fontsize=8)
plt.xticks(fontsize=8)
ticks = [-6,-3,0,3,6]
axs.set_xticks(ticks)
axs.set_yticks(ticks)

cbar = plt.colorbar(sc)
sc.set_clim(-1, 1)
cbar.set_ticks([-1,0,1])
cbar.set_ticklabels(['-1',"0","+1"])
cbar.ax.tick_params(labelsize=8) 

axs.spines[['top', 'right']].set_visible(False)
axs.set_title("specialization index", fontsize=8)
plt.tight_layout()
plt.savefig(args.specialistdir + '/specialists.pdf')
plt.savefig(args.specialistdir + '/specialists.png')
plt.close()