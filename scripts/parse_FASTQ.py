# parse FASTQ files, build consensus sequence for each unique molecular barcode observed at least 3 times,
# and save a file of all genotypes in library with their corresponding number of observations (number of unique molecular barcodes observed with that genotype)

import pysam
import glob
import argparse
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import myutils

# parse command-line arguments
# ============================
parser = argparse.ArgumentParser(description='parse barcodes and genotypes from FASTQ files')
parser.add_argument('sampleid', type=str, help='Sample ID which is prefix for associated FASTQ files')
parser.add_argument('fastqdir', type=str, help='Directory containing the FASTQ files')
parser.add_argument('outputdir', type=str, help='Output directory')
parser.add_argument('--purgereads', type=float, help='purge read frequency')
parser.add_argument('--purgebarcodes', type=float, help='purgebarcode frequency')
parser.add_argument('--seed',help='seed for random.seed()')
parser.add_argument('--maxreads', type=int, help='max reads to parse') 
parser.add_argument('--maxbarcodes', type=int, help='max barcoded consensus sequences to build') 
parser.add_argument('--makebarcodeplots', type=int, default=0, help='set to 1 to output barcode stats pdfs')
parser.add_argument('--ampliconversion', type=str, help='There are three amplicons to align to, depending on the PCR primers used and the trimming parameters', default='longest') 
args = parser.parse_args()
if args.seed:
	random.seed(args.seed)

# Identify the FASTQ files for this sample ID
# ===========================================
r1files = sorted(glob.glob("%s/SRR*-%s_1.fastq.gz" % (args.fastqdir, args.sampleid)))
r2files = sorted(glob.glob("%s/SRR*-%s_2.fastq.gz" % (args.fastqdir, args.sampleid)))
assert len(r1files) == len(r2files) == 1

# Specify alignment/search and trim parameters
# ============================================
# strings used to find perfect match of the initial J sequence for each read; 
# the barcode preceeds this sequence and is of variable length (5-12 by design)
# using a longer or shorter search string may change how many reads are included in the analysis
# note there is a separate barcode for each read direction; these will be concatenated into a single barcode.

# pertinent encoding of all ampliconversion details is here in parse_FASTQ and in myutils.format_genotype()
if args.ampliconversion == 'short':
	refseq = "tcgtcggggaaattgtaaaggCggcgagcgcggcttttccgcgccagcgtgaaAgCGgtgTggacctggcatgtcaacaatacgggagaacacctgtGCcGcctcgttcg"
	r1startseq = "tcgtcggggaaattg".upper()
	r2startseq = "ctggcatgtcaacaa".upper()
	# length of sequence to keep from each read; the kept sequence from reads 1 and 2 wuill be concatenated into a single 'seq'.
	r1keeplength = 65
	r2keeplength = 45
elif args.ampliconversion == 'long':
	refseq = "tcgtcggggaAattgtaaaggCggcgagcgcggcttttccgcgccagcgtgaaAgCggtgTggaccggcggaatTtTtgccgaatGccgtgtggacgtaagcgtgaacgtcaggatcacgtttccccgacccgctggcatgtcaacaatacgggagaacacctgtgCcGcct"
	r1startseq = "tcgtcggggaaattg".upper()
	r2startseq = "cggcggaatTtTtgccg".upper()
	r1keeplength = 65
	r2keeplength = 107
elif args.ampliconversion == 'longest':
	refseq = "tcgtcggggaAattgtaaaggCggcgagcgcggcttttccgcgccagcgtgaaAgCggtgTggactggccgtcaggtacccgtactgtcaccgtgaccgatgaccatccttttgatcgccagatagtggtgcttccggcggaatTtTtgccgaatGccgtgtggacgtaagcgtgaacgtcaggatcacgtttccccgacccgctggcatgtcaacaatacgggagaacacctgtgCcGcctcgttcgccgcgccatcataaatcaccgc"
	r1startseq = "tcgtcggggaaattg".upper()
	r2startseq = "cggcggaatTtTtgccg".upper()
	r1keeplength = 135 
	r2keeplength = 135 
else:
	print("not a valid amplicon version, must be 'short' or 'long' or 'longest'")
	exit()

# open files to write output
# ==========================
myutils.makedir(args.outputdir)
myutils.makedir(args.outputdir+'/genotypecounts/')
myutils.makedir(args.outputdir+'/genotypecounts/countsfiles/')

logfile = open('%s/genotypecounts/countsfiles/%s_logfile.txt' % (args.outputdir, args.sampleid), 'w') 
genotypesfile = open('%s/genotypecounts/countsfiles/%s_genotypecounts.txt' % (args.outputdir, args.sampleid), 'w')


# set up the paired read iterables 
# ================================
'''
re: pysam iterator for paired files. Will need update this code if >1 fastq file per sample per read direction.
for now can just use the globbed r1files[0], BUT
in the future if >1 file per sample... will need to make an iterater function that zips the r1 and r2 files and yields the reads that way,
or just manually merge the fastq files into a single file in preprocessing would probably be easier; 
eg see https://wiki.bits.vib.be/index.php/Combine_the_content_of_several_fastq_files_into_one_big_fastq_file
''' 
print("Opening r1 and r2 files:\n%s\n%s\n" % (r1files[0],r2files[0]))
r1_feed = pysam.FastxFile(r1files[0])
r2_feed = pysam.FastxFile(r2files[0])
paired_read_feed = zip(r1_feed, r2_feed)


# establish bookkeeping data structures
# =====================================
# barcodes is dict of dicts keyed by barcode;
# each barcode key returns a value dictionary containing two keys:
# 	'seqs' returns list of sequences observed with the barcode;
# 	'consensus' returns the consensus sequence for the barcode
barcodes = {}
# i am also curious about the distribution of barcode lengths.
bc1len_dist = {}
bc2len_dist = {}
# counters for number of read pairs either kept or discarded for various reasons:
n_reads = { 
	'total':0,
	'contain N':0,
	'failed substring search':0, # search for perfect J-sequence read start failed
	'improper barcode length':0, # <5 or >12
	'retained':0, # should equal total minus the rest but will tally anyway
	'purged':0, # for optional purging of reads at a specified frequency to subsample
	'fail filter':0
	}

# iterate through the paired reads:
# =================================
for (read1, read2) in paired_read_feed:

	n_reads['total'] += 1
	if n_reads['total'] % 5000000 == 0: 
		print("Parsing read pair #%i" % n_reads['total'])

	if args.maxreads:
		if n_reads['total'] == args.maxreads + 1:
			print("Reached maxreads of %s and will stop parsing reads." % args.maxreads)
			break

	if args.purgereads:
		if random.random() < args.purgereads:
			n_reads['purged'] += 1
			continue
	
	assert read1.name == read2.name 
	
	# discard any read failing the illumina filter:
	# read.comment[2] designates illumina filter (Y if fails filter)
	if read1.comment[2] == 'Y' or read2.comment[2] == 'Y':
		n_reads['fail filter'] += 1
		continue

	# discard any read containing an N:
	if "N" in read1.sequence or "N" in read2.sequence:
		n_reads['contain N'] += 1
		continue
	
	# find the beginning of J sequence to identify barcode length,
	# extract the barcode and the amplicon sequence for each read
	try: 
		bc1len = read1.sequence.index(r1startseq)
	except ValueError:
		n_reads['failed substring search'] += 1
		continue
	bc1 = read1.sequence[:bc1len]
	r1seq = read1.sequence[bc1len:] # ie barcode trimmed
	try:
		bc2len = read2.sequence.index(r2startseq)
	except ValueError:
		n_reads['failed substring search'] += 1
		continue
	bc2 = read2.sequence[:bc2len]
	r2seq = read2.sequence[bc2len:] # ie barcode trimmed

	# tally distributions of barcode lengths observed (but will remove some as below)
	if bc1len in bc1len_dist:
		bc1len_dist[bc1len] += 1
	else:
		bc1len_dist[bc1len] = 1
	if bc2len in bc2len_dist:
		bc2len_dist[bc2len] += 1
	else:
		bc2len_dist[bc2len] = 1

	# rarely, length of either side barcode will be less than 5 or greater than 12;
	# I kept track of these in the barcode length distribution above to see how often, but
	# I will not include them downstream from here in the analysis:
	if bc1len > 12 or bc1len < 5:
		n_reads['improper barcode length'] += 1
		continue
	if bc2len > 12 or bc2len < 5:
		n_reads['improper barcode length'] += 1
		continue

	# Concatenate the read 1 and read 2 barcodes into a single barcode, 
	# concatenate the R1 and R2 J sequence into a single seq,
	# and save barcode : seq in barcodes dict.
	# note we will have multiple sequences for most barcodes (with sufficient sequencing depth),
	# which are all derived from different reads, and are mostly going to be the same, except some
	# will have errors during round 2 PCR and/or during the sequencer run which we will eventually 
	# correct for by creating a consensus sequence.
	concat_barcode = bc1 + bc2
	concat_seq = r1seq[:r1keeplength] + r2seq[:r2keeplength] # only retain this amount of each seq, and concatenate
	
	n_reads['retained'] += 1
	if concat_barcode in barcodes:
		barcodes[concat_barcode]['seqs'].append(concat_seq)
	else:
		barcodes[concat_barcode] = {'seqs':[concat_seq],'consensus':""}

# Done iterating through readpairs; we now have a dict of barcodes with associated sets of sequences
# ==================================================================================================
sorted_bc1len_dist = sorted(bc1len_dist.items(), key=lambda x: x[1], reverse=True)
sorted_bc2len_dist = sorted(bc2len_dist.items(), key=lambda x: x[1], reverse=True)
logfile.write("after parsing FASTQ_files, initially saw these distributions for barcode lengths -- but note that if a read had barcode lengths outside of 5-12 on either side, that read pair is discarded:\n")
logfile.write("barcode 1 length distribution:\n%s\n" % sorted_bc1len_dist)
logfile.write("barcode 2 length distribution:\n%s\n" % sorted_bc2len_dist)
logfile.write("Discarded %s reads that failed illumina filter\n" % n_reads['fail filter'])
logfile.write("\nOut of %i read pairs:\n\tpurged %s,\n\tdiscarded %i for containing N,\n\tdiscarded %i for improper barcode length,\n\tdiscarded %i for failed substring search.\n" % 
		(n_reads['total'], n_reads['purged'], n_reads['contain N'], n_reads['improper barcode length'], n_reads['failed substring search']))
logfile.write("\nOverall retained %i read pairs containing %i unique barcodes, prior to consensus building.\n" % (n_reads['retained'], len(barcodes)))

	
n_reads_per_bc_dist = {} # key = # of reads; value = # of barcodes with that many reads
genotype_counts = {} # key = consensus sequence (genotype); value = # of barcodes counted for that genotype (counts)
num_barcodes_analyzed = 0
num_barcodes_retained = 0

# Now iterate through barcodes and build consensus sequences
for barcode in barcodes:

	if args.purgebarcodes:
		if random.random() < args.purgebarcodes:
			continue

	num_barcodes_analyzed += 1
	
	if num_barcodes_analyzed % 5000000 == 0: 
		print("Analyzing barcode #%i" % num_barcodes_analyzed)

	if args.maxbarcodes:
		if num_barcodes_retained == args.maxbarcodes + 1:
			print("Reached maxbarcodes of %s and will stop building consensus sequences" % args.maxbarcodes)
			break
	
	n_reads_this_bc = len(barcodes[barcode]['seqs'])
	if n_reads_this_bc in n_reads_per_bc_dist:
		n_reads_per_bc_dist[n_reads_this_bc] += 1
	else:
		n_reads_per_bc_dist[n_reads_this_bc] = 1

	if n_reads_this_bc >= 3:
		num_barcodes_retained += 1
		# build consensus sequence
		for pos in range(len(barcodes[barcode]['seqs'][0])):
			char_counts = {}
			for seq in barcodes[barcode]['seqs']:
				if seq[pos] in char_counts:
					char_counts[seq[pos]] += 1
				else:
					char_counts[seq[pos]] = 1
			barcodes[barcode]['consensus'] += max(char_counts, key=char_counts.get) # tie goes to the runner

		if barcodes[barcode]['consensus'] in genotype_counts:
			genotype_counts[barcodes[barcode]['consensus']] += 1
		else:
			genotype_counts[barcodes[barcode]['consensus']] = 1

sorted_n_reads_per_bc_dist = sorted(n_reads_per_bc_dist.items(), key=lambda x: x[1], reverse=True)
sorted_genotypes = sorted(genotype_counts.items(), key=lambda x: x[1], reverse=True)

logfile.write("\ndistribution of number of reads per observed barcode (note consensus was only built for >= 3 reads per barcode):\n%s\n" % sorted_n_reads_per_bc_dist)
logfile.write("\nnumber of retained barcodes, ie, number of consensus sequences built: %i\n" % num_barcodes_retained)

for (genotype,count) in sorted_genotypes:
	genotypesfile.write("%s %i\n" % (genotype,count))

logfile.close()
genotypesfile.close()


if args.makebarcodeplots:
	# Plot summary statistics for barcode lengths and reads per barcode:
	#===================================================================
	fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(9, 3))
	fig.suptitle('%s barcode stats' % args.sampleid)
	xs = [x for (x,y) in sorted_bc1len_dist]
	ys = [y for (x,y) in sorted_bc1len_dist]
	axs[0].bar(xs,ys)
	axs[0].set_title("R1 bc lengths")
	axs[0].set_xlim(3.5,13.5)
	xs = [x for (x,y) in sorted_bc2len_dist]
	ys = [y for (x,y) in sorted_bc2len_dist]
	axs[1].bar(xs,ys)
	axs[1].set_title("R2 bc lengths")
	axs[1].set_xlim(3.5,13.5)
	xs = [x for (x,y) in sorted_n_reads_per_bc_dist]
	ys = [y for (x,y) in sorted_n_reads_per_bc_dist]
	axs[2].bar(xs,ys)
	axs[2].set_title("# reads per bc")
	axs[2].set_xlim(0,20)
	for ax in axs:
		ax.spines[['top', 'right']].set_visible(False)
	plt.tight_layout()
	myutils.makedir(args.outputdir+'/genotypecounts/barcodestats/')
	plt.savefig('%s/genotypecounts/barcodestats/%s_barcodestats.pdf' % (args.outputdir, args.sampleid))

# convert the raw sequence/counts file into 
# a formatted genotype string counts & frequency csv
#===================================================

infilepath = '%s/genotypecounts/countsfiles/%s_genotypecounts.txt' % (args.outputdir, args.sampleid)
outfilepath = '%s/genotypecounts/countsfiles/%s_formattedcounts.csv' % (args.outputdir, args.sampleid)
print("now formatting genotypes from input counts file %s.\nwill write results to %s." % (infilepath, outfilepath))
outfile_feed = open(outfilepath, "w")
outfile_feed.write("genotype,counts,frequency(%)\n")
infile_feed = open(infilepath, "r")

# first iterate all genotypes to get total counts
totalcounts = 0
for line in infile_feed:
	(seq,counts) = line.split()
	counts = int(counts)
	totalcounts += counts

infile_feed.close() # i think this needs to be closed and re-opened to iterate through again
infile_feed = open(infilepath, "r")

for line in infile_feed:
	(seq,counts) = line.split()  
	counts = float(counts)
	percentage = (counts/totalcounts)*100
	counts = int(counts)

	(this_formatted_genotype, this_formatted_additional_muts) = myutils.format_genotype(seq, refseq, args.ampliconversion)
	if len(this_formatted_additional_muts) > 0:
		outfile_feed.write("[%s]+%s, %s, %f\n" % 
			(this_formatted_genotype, this_formatted_additional_muts, counts,percentage))
	else:
		outfile_feed.write("[%s], %s, %f\n" % 
			(this_formatted_genotype, counts, percentage))

outfile_feed.close()
infile_feed.close()