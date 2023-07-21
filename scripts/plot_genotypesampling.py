import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
import myutils

parser = argparse.ArgumentParser(description='plot library genotype sampling for a set of samples.')
parser.add_argument('outputdir', type=str, help='path to output directory')
parser.add_argument('sampleids', type=str, help='list of sample IDs', nargs='+')
parser.add_argument('nrows', type=int, help='number of rows for subplots')
parser.add_argument('ncols', type=int, help='number of cols for subplots')
parser.add_argument('--nplots', type=int, default=4, help='number of subplots to make for rankedorder mutation plots.')
parser.add_argument('--plotfileprefix', type=str, default='', help='prefix for the saved plots files')
parser.add_argument('--xlim', type=int, default=500, help='x axis limit (number of genotypes)')
parser.add_argument('--ylim', type=int, default=5, help='y axis limit (genotype frequency)')


args = parser.parse_args()

outputsubdir = args.outputdir + '/genotypecounts/' # should already be made by parse_FASTQ run with the given outputdir
countsfilesdir = args.outputdir + '/genotypecounts/countsfiles/' # should already be made by parse_FASTQ
rankordereddir = outputsubdir + '/rankorderplots/'
myutils.makedir(rankordereddir)

samples = args.sampleids


# set up subplots for the rank-ordered genotype frequency plots:
assert args.nrows * args.ncols == args.nplots
fig, axes = plt.subplots(nrows=args.nrows, ncols=args.ncols, figsize=(4*args.ncols,4*args.nrows)) # access with axes[irow,icol] unless ravel.
if len(samples) == 1:
		axes = [axes]
else:
	axes = axes.ravel() # axs.ravel() converts 2-dim axes object to a list in row-major order
schema = 'prog-stop-LF' # coloring scheme based on genotype, defined in myutils.encode_genotype_color(). by default, library=purple, stopcodon=red, L>F mutation teal
schemadir = rankordereddir + '/' + schema 
myutils.makedir(schemadir)

for i_sample, sample in enumerate(samples): 
	if i_sample < args.nplots:
		# plot the frequencies of each of the observed genotype in a rank-ordered plot
		genotype_counts_datalist = myutils.read_genotype_counts_file("%s/%s_formattedcounts.csv" % (countsfilesdir, sample))
		xs = [i for i in range(1,len(genotype_counts_datalist)+1)]
		ys = [f for (g,c,f) in genotype_counts_datalist]

		colors = [myutils.encode_genotype_color(g,schema) for (g,c,f) in genotype_counts_datalist]
		axes[i_sample].bar(xs,ys,color=colors,width=1)
		axes[i_sample].spines[['top', 'right']].set_visible(False)
		axes[i_sample].set_title(sample, fontsize=14)

		axes[i_sample].xaxis.set_tick_params(labelsize=12)
		axes[i_sample].yaxis.set_tick_params(labelsize=12)
	
fig.supxlabel('Ranked genotypes')
fig.supylabel('Frequency (%)')

for scale in ['linear']: 
	outputdir = schemadir + '/' + scale 
	myutils.makedir(outputdir)
	#[axis.set_yscale(scale) for axis in axes]
	[axis.set_xlim(0,args.xlim) for axis in axes]
	[axis.set_ylim(0,args.ylim) for axis in axes]
	plt.tight_layout()
	plt.savefig('%s/genotypesampling-%s-%s.pdf' % (outputdir, args.plotfileprefix, scale))
	plt.savefig('%s/genotypesampling-%s-%s.png' % (outputdir, args.plotfileprefix, scale))
