import argparse
import myutils
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

def hybrid(gt1, gt2):
	# "combine" mutations from two genotypes
	# eg, any mutations from the ancester generalist sequence in either genotype are present in the hybrid
	hybridgt = ''
	for pos in [1,2,3,4,5,6,7,8,9]:
		if int(gt1[pos]) + int(gt2[pos])  > 0:
			hybridgt = hybridgt + '1'
		else:
			hybridgt = hybridgt + '0'
	#print("the hybrid between:\n%s\n%s\n%s" % (gt1[1:-1], gt2[1:-1], hybridgt))
	return hybridgt

parser = argparse.ArgumentParser(description='analyze hybrids from pairs of specialists')
parser.add_argument('srfilesdir', type=str, help='directory containing the selection rate files')
parser.add_argument('specialistdir', type=str, help='directory containing specialist genotype classifications')
args = parser.parse_args()

# L and O specialists:
L_specialists_file = args.specialistdir + '/Lspecialists.csv'
O_specialists_file = args.specialistdir + '/Ospecialists.csv'
L_specialists_d = myutils.read_selectionrate_file(L_specialists_file, as_dict=True)
L_specialists_l = myutils.read_selectionrate_file(L_specialists_file, as_dict=False)
O_specialists_d = myutils.read_selectionrate_file(O_specialists_file, as_dict=True)
O_specialists_l = myutils.read_selectionrate_file(O_specialists_file, as_dict=False)
print("There are %s L specialists and %s O specialists\n" % (len(L_specialists_l), len(O_specialists_l)))

# landscapes on L and O:
sr_on_L_file = '%s/sr_LibAvg_vGen_LamB.csv' % args.srfilesdir
sr_on_O_file = '%s/sr_LibAvg_vGen_OmpF.csv' % args.srfilesdir
sr_on_L = myutils.read_selectionrate_file(sr_on_L_file, as_dict=True)
sr_on_O = myutils.read_selectionrate_file(sr_on_O_file, as_dict=True)


# subset on specialists that are:
# 1. measured on both L and O, and
# 2. equal number between the two; this will be defined by the number of L specialists since there are fewer.

# 1. keep specialists that have a measured selection rate on each receptor:
l_specs_dual = []
o_specs_dual = []
for ls in L_specialists_d:
	if np.isfinite(sr_on_O[ls]):
		l_specs_dual.append(ls)
for os in O_specialists_d:
	if np.isfinite(sr_on_L[os]):
		o_specs_dual.append(os)
print("there are %s L-specialists that have a measured selection rate on OmpF" % len(l_specs_dual))
print("there are %s O-specialists that have a measured selection rate on LamB" % len(o_specs_dual))

# 2. narrow the list of O-specialists to the top X
num_keep = len(l_specs_dual)
sorted_O_specialists_list = sorted(O_specialists_l, key=lambda x: x[1], reverse=True)[:num_keep]
o_specs_dual = [tup[0] for tup in sorted_O_specialists_list]

counts = {'total':0, # total number of O/L specialist pairs analyzed
		  'redundant_hybrid_is_Lsp':0, # the hybrid of these two specialists is one of these specialists (mutation set of one is subset of other)
		  'redundant_hybrid_is_Osp':0, # the hybrid of these two specialists is one of these specialists (mutation set of one is subset of other)
		  'spec_nan':0, # at least one of the specialists has an unobserved sr on at least one receptor; cannot compute
		  'hybrid_nan_on_either':0, # the hybrid is unobserved on at least one receptor, so cannot compute vector
		  'hybrid_nan_on_L':0, # subset of hybrid_nan, if nan on L (nonexclusive) -
		  'hybrid_nan_on_O':0, # subset of hybrid_nan, if nan on O (nonexclusive)
		  'hybrid_nan_on_both':0, # subset of hybrid_nan, if nan on L and O
		  'hybrid_nan_is_Lspec':0, # some hybrids cannot have epsilon computed because they are NaN on OmpF but they are Lspecialists
		  'hybrid_nan_is_Ospec':0, # some hybrids cannot have epsilon computed because they are NaN on LamB but they are Ospecialists
		  'triad_analyzed':0,
		  'hybrid_is_Lsp':0,
		  'hybrid_is_Osp':0}

slope_m_l = [] # slope m of y=mx+b connecting two specialists
xs_ospecs = []
ys_ospecs = []
xs_lspecs = []
ys_lspecs = []

# bookkeeping calculations of the vector epsilon (geometric interpretation of change in fitness in hybrid from two parental specialists)
epsilon_i_l = []
epsilon_j_l = []
epsilon_abs_l = []
epsilon_signed_l = []
epsilon_heron_abs_l = []
epsilon_heron_signed_l = []

distribution_avgH_d_avgS = []

hybrid_delta_on_L_l = []
hybrid_delta_on_O_l = []
hybrid_delta_colors = []

print("Analyzing %s Ls and %s O sp" % (len(l_specs_dual), len(o_specs_dual)))
for Lspec in l_specs_dual:
	for Ospec in o_specs_dual:
		counts['total'] += 1
		hybridgt = '[' + hybrid(Lspec, Ospec) + ']'

		if hybridgt == Lspec: 
			counts['redundant_hybrid_is_Lsp'] += 1
			continue
		if hybridgt == Ospec:
			counts['redundant_hybrid_is_Osp'] += 1
			continue

		if np.isnan(sr_on_L[Lspec]) or \
			np.isnan(sr_on_O[Lspec]) or \
			np.isnan(sr_on_L[Ospec]) or \
			np.isnan(sr_on_O[Ospec]):
			counts['spec_nan'] += 1 # should be zero, already checked for this when building the lists l_specs_dual, o_specs_dual.
			continue

		if np.isnan(sr_on_L[hybridgt]) or np.isnan(sr_on_O[hybridgt]): 
			counts['hybrid_nan_on_either'] += 1 
			if np.isnan(sr_on_L[hybridgt]):
				counts['hybrid_nan_on_L'] += 1
			if np.isnan(sr_on_O[hybridgt]):
				counts['hybrid_nan_on_O'] += 1
			if np.isnan(sr_on_L[hybridgt]) and np.isnan(sr_on_O[hybridgt]): 
				counts['hybrid_nan_on_either'] += 1
			# some of these hybrids will actually be L or O specialists that are nan on opposite receptor; not counted now
			if hybridgt in L_specialists_d:
				counts['hybrid_nan_is_Lspec'] += 1
			if hybridgt in O_specialists_d:	
				counts['hybrid_nan_is_Ospec'] += 1
			continue 

		# ok now we can fully analyze this triad by computing epsilon:
		counts['triad_analyzed'] += 1

		# calculate epsilon:
		
		# 1. orthogonal vector method
		# notation: A is Ospec, B is Lspec, H is hybrid; X-axis is sr_on_L; Y-axis is sr_on_O
		xa = sr_on_L[Ospec]
		ya = sr_on_O[Ospec]
		xb = sr_on_L[Lspec]
		yb = sr_on_O[Lspec]
		xh = sr_on_L[hybridgt]
		yh = sr_on_O[hybridgt]

		# separately: the effect on hybridization on fitness on each receptor:
		# eg, how much more or less fit is the hybrid than the average of the specialists on that receptor:
		hybrid_delta_on_L = xh - (xa+xb)/2
		hybrid_delta_on_O = yh - (ya+yb)/2
		hybrid_delta_on_L_l.append(hybrid_delta_on_L)
		hybrid_delta_on_O_l.append(hybrid_delta_on_O)

		# is the hybrid a specialist?
		if hybridgt in L_specialists_d:
			counts['hybrid_is_Lsp'] += 1
			c = 'tab:blue'
		elif hybridgt in O_specialists_d:
			counts['hybrid_is_Osp'] += 1
			c = 'tab:red'
		else:
			c = 'lightgrey'
		hybrid_delta_colors.append(c)

		m = (yb-ya)/(xb-xa)

		spec_average_fitness = (xa+xb+ya+yb)/4.0
		hybrid_average_fitness = (xh+yh)/2.0
		distribution_avgH_d_avgS.append(hybrid_average_fitness - spec_average_fitness)

		xs_ospecs.append(xa)
		ys_ospecs.append(ya)
		xs_lspecs.append(xb)
		ys_lspecs.append(yb)
		slope_m_l.append(m)

		b = ya-m*xa
		# x2,y2 is the point on the line connecting the specialists that the orthogonal vector to H intersects
		x2 =(xh+m*(yh-b))/(1+m*m)
		y2 =(m*xh+m*m*yh+b)/(1+m*m)
		epsilon_i = xh-x2
		epsilon_j = yh-y2
		epsilon_abs = np.sqrt(epsilon_i*epsilon_i + epsilon_j*epsilon_j) 
		d = (xh-xa)*(yb-ya) - (yh-ya)*(xb-xa)
		epsilon_signed = epsilon_abs*np.sign(-d)
		epsilon_i_l.append(epsilon_i)
		epsilon_j_l.append(epsilon_j)
		epsilon_abs_l.append(epsilon_abs)
		epsilon_signed_l.append(epsilon_signed)

		# 2. Heron's formula method; this arrives at the same results as the orthogonal vector method
		heron_a = np.sqrt(np.square(ya-yh)+np.square(xa-xh))
		heron_b = np.sqrt(np.square(yb-ya)+np.square(xb-xa))
		heron_c = np.sqrt(np.square(yb-yh)+np.square(xb-xh))
		heron_x = (heron_a+heron_b+heron_c)/2.0
		heron_epsilon = 2.0*np.sqrt( heron_x*(heron_x - heron_a)*(heron_x - heron_b)*(heron_x - heron_c) )/heron_b
		epsilon_heron_abs_l.append(heron_epsilon)
		# determine sign of epsilon using y=mx+b form of line connecting two specialists; plug in xh and compare to yh
		y_test = m*xh+b
		if y_test < yh:
			sign = 1
		elif y_test > yh:
			sign = -1
		else:
			sign = 0
			assert epsilon_heron == 0

		epsilon_heron_signed_l.append(sign*heron_epsilon)

print(counts)

# plot histrogram of epsilon for main text figure
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(3,3))

# histogram signed epsilons
num_neg_epsilon = sum(i < 0 for i in epsilon_heron_signed_l)
num_pos_epsilon = sum(i > 0 for i in epsilon_heron_signed_l)
total_epsilon = num_neg_epsilon + num_pos_epsilon
frac_neg_epsilon = num_neg_epsilon/total_epsilon
frac_pos_epsilon = num_pos_epsilon/total_epsilon
axes.hist(epsilon_heron_signed_l, density=True, bins=30, alpha=0.55, facecolor='k', edgecolor='k', linewidth=0.5, label='signed epsilon, heron')
axes.axvline(x=0, ls=':', color='k')
axes.text(-1.75,0.3,"%.0f%%" % (frac_neg_epsilon*100), fontsize=8)
axes.text(1,0.3,"%.0f%%" % (frac_pos_epsilon*100), fontsize=8)
axes.set_xlabel("hybrid change in fitness\nfrom specialist average", fontsize=8)
axes.set_ylabel("density", fontsize=8)
axes.spines[['top', 'right']].set_visible(False)
#axes.legend(fontsize=6)
axes.tick_params(labelsize=8) 
plt.tight_layout()
plt.savefig("%s/hybrid_change_in_fitness.pdf" % args.specialistdir)
plt.savefig("%s/hybrid_change_in_fitness.png" % args.specialistdir)



# plot scatter of specialists included in hybrid analyses and hybrid change scatter for supplement
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(6,3))

# The L and O specs used for triad analysis:
axes[0].axline((0, 0), slope=1, color = "grey", linestyle = ":", linewidth = 1.07/1.5)
axes[0].scatter(xs_lspecs,ys_lspecs,s=12,linewidth=0.5, alpha=0.25, edgecolor='k', c='tab:blue', label='Lspecialists')
axes[0].scatter(xs_ospecs,ys_ospecs,s=12,linewidth=0.5, alpha=0.25, edgecolor='k', c='tab:red', label='Ospecialists')
axes[0].set_xlabel("sr on LamB")
axes[0].set_ylabel("sr on OmpF")
axes[0].set_xlim(-10,10)
axes[0].set_ylim(-10,10)
axes[0].spines[['top', 'right']].set_visible(False)
axes[0].legend(fontsize=6)
axes[0].tick_params(labelsize=8) 

# The effect of hybridization on fitness for each receptor compared to the average fitness on that receptor for the two specialists
axes[1].scatter(hybrid_delta_on_L_l, hybrid_delta_on_O_l, s=10, linewidth=0.5, alpha=0.35, edgecolor='k',c=hybrid_delta_colors)
axes[1].axvline(x=0, ls=':', color='k')
axes[1].axhline(y=0, ls=':', color='k')
axes[1].set_xlabel('hybrid change in Sr_L from avg', fontsize=8)
axes[1].set_ylabel('hybrid change in Sr_O from avg', fontsize=8)
axes[1].tick_params(labelsize=8) 
axes[1].set_xlim(-10,10)
axes[1].set_ylim(-10,10)

plt.tight_layout()
plt.savefig("%s/hybridsummary.pdf" % args.specialistdir)
plt.savefig("%s/hybridsummary.png" % args.specialistdir)