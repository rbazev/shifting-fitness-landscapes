import myutils
import numpy as np
import argparse
parser = argparse.ArgumentParser(description='merge data from multiple selection rate files')

parser.add_argument('--mergedfile', type=str, help='path to merged sr file')
parser.add_argument('--srfiles', type=str, help='paths to srfiles to merge', nargs = '*')
args = parser.parse_args()
outfile = open(args.mergedfile, 'w')
outfile.write('genotype,s_r_merged\n')

merged_data = {} # the merged selection rate, keyed by genotype. 

numfiles = len(args.srfiles)
srdicts = [myutils.read_selectionrate_file(srfile, as_dict=True) for srfile in args.srfiles]

# First, build a list of all genotypes; this is the programmed genotypes plus any additional genotypes with selection rate measured
# in at least one of the files being merged
gts = myutils.construct_genotype_space()
gts_set = set(gts)
union = gts_set.union(gts_set) # this is silly making a union with itself but it works. could just be union = gts_set
for srd in srdicts:
	union = union.union(set(srd.keys())) # collect all observed genotypes into the set 'union'
prog_and_addition_genotypes = [] 
for gt in myutils.construct_genotype_space():
	prog_and_addition_genotypes.append(gt) # first populate with the programmed genotypes,
for gt in union: # then add non-programmed genotypes:
	if gt not in prog_and_addition_genotypes: prog_and_addition_genotypes.append(gt)

for gt in prog_and_addition_genotypes: 
	# a list of all of the observations for selection rate for this genotype:
	observations = []
	for srdict in srdicts:
		if gt in srdict and not np.isnan(srdict[gt]):
			observations.append(srdict[gt])

	if len(observations) < 1: 
		merged_data[gt] = np.nan
	else:
		merged_data[gt] = np.mean(observations)
	
	outfile.write('%s,%s\n' % (gt,merged_data[gt]))

outfile.close()