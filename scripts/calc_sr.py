# calculate selection rates for all genotypes, from pre- and post-competition sequencing counts
# the pre- and post- counts files are hardcoded at end

import itertools
import numpy as np
import myutils
import argparse

parser = argparse.ArgumentParser(description='calculate selection rates')
parser.add_argument('outputdir', type=str, help='top output directory containing genotypecounts directory')
parser.add_argument('--includeunprogrammed', type=int, default=0, help='set as 1 to include UNPROGRAMMED genotypes')
args = parser.parse_args()
print('about to calculate selection rates. using output directory %s' % args.outputdir)


def calc_sR(outputdir,sample_post,sample_pre):

	print('computing sr for post-selection %s vs. pre-selection %s' % (sample_post,sample_pre))
	# selection rate, the difference of malthusian growth rates
	# where f_i,t = genotype frequency of variant i at time t
	# sR of variant i = ln(f_i,4/f_i,0) - ln(f_wt,4/f_wt,0)

	# sample_pre is the sample name for pre-selection 
	# sample_post is the sample name for post-selection 

	# output a csv file of genotype,sR
	myutils.makedir('%s/selectionrate/' % outputdir)
	myutils.makedir('%s/selectionrate/srfiles' % outputdir)

	if args.includeunprogrammed:
		outfile_name = '%s/selectionrate/srfiles/sr_%s_unprogrammed.csv' % (outputdir, sample_post)
	else:
		outfile_name = '%s/selectionrate/srfiles/sr_%s.csv' % (outputdir, sample_post)
	outfile = open(outfile_name, 'w')
	outfile.write('genotype,s_r vs. %s\n' % sample_pre)

	sample_pre_file = '%s/genotypecounts/countsfiles/%s_formattedcounts.csv' % (outputdir, sample_pre)
	sample_post_file = '%s/genotypecounts/countsfiles/%s_formattedcounts.csv' % (outputdir, sample_post)

	genotype_percent_pre = myutils.read_genotype_counts_file(sample_pre_file, as_dict=True)
	genotype_percent_post = myutils.read_genotype_counts_file(sample_post_file, as_dict=True)

	genotype_counts_pre = myutils.read_genotype_counts_file(sample_pre_file, as_counts_dict=True)
	genotype_counts_post = myutils.read_genotype_counts_file(sample_post_file, as_counts_dict=True)

	assert '[000000000]' in genotype_percent_pre
	if '[000000000]' not in genotype_percent_post:
		print("CANNOT MEASURE SR BECAUSE [000000000] IS NOT OBSERVED IN POST") 
		# cannot measure selection rate (as defined) if the wt is not observed

	gts = myutils.construct_genotype_space()
	
	if args.includeunprogrammed == 1:
		additional_gts_post = [gt for gt in genotype_percent_post if gt not in gts ] # extra gts obsvered in post
		additional_gts_both = [gt for gt in additional_gts_post if gt in genotype_percent_pre] # extra gts observed in pre and post
		gts = gts + additional_gts_both

	for gt in gts : 

		try:
			f_i4 = genotype_percent_post[gt]
			f_i0 = genotype_percent_pre[gt]
			f_wt4 = genotype_percent_post['[000000000]']
			f_wt0 = genotype_percent_pre['[000000000]']
			#print("f_i4, f_i0, f_wt4, f_wt0: %s, %s, %s, %s" % (f_i4, f_i0, f_wt4, f_wt0))
			#input()
			sr = np.log(f_i4/f_i0) - np.log(f_wt4/f_wt0)

		except KeyError:
			# this genotype wasn't observed in either the pre- or post-
			# selection sequencing experiments. 
			# in most cases, the given genotype was wiped out below LOD by competition; 
			# in some cases it may be possible the genotype is not observed in the pre-selection library.
			# will give sr = NaN
			sr = np.nan

		outfile.write('%s,%s\n' % (gt,sr))

	pass
	outfile.close()

samples = [
			'LibA_pre',
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
			'LibB_pre',
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

for sample in samples:
	if 'pre' in sample:
		continue
	comparison_sample = sample[:4] + '_pre'
	calc_sR(args.outputdir, sample, comparison_sample)