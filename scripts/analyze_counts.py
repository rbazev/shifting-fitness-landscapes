import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
import myutils

parser = argparse.ArgumentParser(description='summarize information on read/barcode/genotype counts across a list of samples ')
parser.add_argument('outputdir', type=str, help='path to output directory')
parser.add_argument('sampleids', type=str, help='list of sample IDs', nargs='+')
args = parser.parse_args()

outputsubdir = args.outputdir + '/genotypecounts/' # should already be made by parse_FASTQ run with the given outputdir
countsfilesdir = args.outputdir + '/genotypecounts/countsfiles/' # should already be made by parse_FASTQ
depthstatsdir = args.outputdir + '/genotypecounts/depthstats/'
myutils.makedir(depthstatsdir)

samples = args.sampleids

# bookkeeping for total/retained reads/barcodes counting:
depthdata = {} # keyed by sample, value is another dict keyed by 'TR', 'RR', 'TB', 'RB' for total/retained reads/barcodes
depthdata_outfile = depthstatsdir+'depth_summary.csv'
depthdata_out = open(depthdata_outfile, 'w')
depthdata_out.write("sample,total_reads,retained_reads,total_barcodes,retained_barcodes\n")

# bookkeeping for x/256 or x/512 variants observed in each sample:
librarysampling_outfile = depthstatsdir + '/librarysampling.txt'
librarysampling_out = open(librarysampling_outfile, 'w')

for i_sample, sample in enumerate(samples): 
	# extract total/retained read and barcode counts from the sample's logfile
	# summarise this in the depth_summary.csv and also will use the data to plot later.
	filename = countsfilesdir + '%s_logfile.txt' % sample
	file = open(filename,'r')
	lines = file.readlines()
	total_reads = 0
	retained_reads = 0
	total_barcodes = 0
	retained_barcodes = 0
	for line in lines:
		if "Out of " in line[:10]:
			total_reads = [int(s) for s in line.split() if s.isdigit()][0]
		elif "Overall retained " in line[:20]:
			retained_reads = [int(s) for s in line.split() if s.isdigit()][0]
			total_barcodes = [int(s) for s in line.split() if s.isdigit()][1]
		elif "number of retained barcodes" in line:
			retained_barcodes = [int(s) for s in line.split() if s.isdigit()][0]
	depthdata[sample] = { 	'TR':total_reads, 'RR':retained_reads,
							'TB':total_barcodes, 'RB':retained_barcodes }
	depthdata_out.write("%s,%s,%s,%s,%s\n" % (sample, total_reads, retained_reads, total_barcodes, retained_barcodes))

	# summarize mutation sampling (X/256, X/512) from the genotype counts files
	genotype_counts_datadict = myutils.read_genotype_counts_file("%s/%s_formattedcounts.csv" % (countsfilesdir, sample), as_dict=True)

	if '[00000000]' in genotype_counts_datadict: length=8
	elif '[000000000]' in genotype_counts_datadict: length=9
	elif '[010000001]' in genotype_counts_datadict: length=9
	else: 
		print("could not detect length of binary genotype, skipping analysis of mutation sampling and some plots for %s" % sample)
		continue

	total_possible = 512 if length==9 else 256
	muts_sampled = 0
	muts_missing = 0
	for gt in myutils.construct_genotype_space(length=length):
		if gt not in genotype_counts_datadict:
			muts_missing += 1
		else:
			muts_sampled += 1
	assert (muts_missing + muts_sampled) == total_possible 
	librarysampling_out.write("%s: %i/%i genotype variants observed.\n" % (sample,muts_sampled,total_possible))


