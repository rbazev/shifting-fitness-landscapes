# misc functions often used in multiple scripts

import itertools
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['font.family'] = 'Arial'
import numpy as np
import networkx as nx

# RELATING TO THE EXPERIMENTAL CONDITIONS / LIBRARY DESIGN / GENOTYPE SPACE:
# L-spec competitior = [10000001] (100000010) G1, G8
# O-spec competitor = G1,3,6,8,9; [10100101]1 == (101001011)

def competition_factors(sample):

	factors = {	'host' : '',
				'competitor' : '',
				'library' : '' }

	# sample is string in the form 'LibXvCompX_(singlehost)_rep#'

	if 'OmpF' in sample:
		factors['host'] = 'L-'
	elif 'LamB' in sample:
		factors['host'] = 'O-'
	else:
		factors['host'] = 'L-+O-'

	if 'vGen' in sample:
		factors['competitor'] = 'EvoC'
	elif 'vLspec' in sample:
		factors['competitor'] = 'Lspec'
	elif 'vOspec' in sample:
		factors['competitor']  = 'Ospec'
	else:
		raise hostascertainmenterror

	if 'LibA' in sample:
		factors['library'] = 'LibA'
	elif 'LibB' in sample:
		factors['library'] = 'LibB'
	else:
		raise libraryascertainmenterror

	assert factors['host']
	assert factors['competitor']
	assert factors['library']

	return factors

def format_genotype(genotype_seq, refseq, ampliconversion):
	'''
	Takes a genotype sequence in the format parsed by parse_FASTQ and
	formats a genotype string as a binary yes (1) / no (0) for carrying mutation
	at each of the mutation sites in mut_sites and retuns a tuple containing:
	1. a string designating mutation status across the designated mutation sites
		eg. '000000000' = WT; '001100010' = mutation at sites 3, 4, and 8, etc
	2. a list of additional mutations in the format (position,mutation)
		eg [(1,'A'),(3,'C')]

	Thus return structure is:
		('110001010',[(4,'A'),(5,'T')])

	If a mutation site is something besides the WT/engineered mutation, 
	the string of 9 will list the nucleotide instead of a 0 (ancestral) or 1 (engineered mutation).
	'''
	# note the read 2 portion of the read is always in reverse complement with respect to the J gene

	#print("The reference sequence, and observed genotype from sequencing, are as below:")
	#print("ref>%s" % refseq)
	#print("gen>%s" % genotype_seq)
	refseq = refseq.upper()

	# mutsites is a list of tuples defining:
	# 1. the index of the mutation site (we call them site 1, site 2, ... site 8, (site 9))
	# 2. the position within the concatenated read; dependent on the parameters r1keeplength = 65 and r2keeplength = 45 or 107 (and the amplicon version).
	# 3. the "wild-type" nucleotide (EvoC sequence); this should be the nucleotide listed in the referense sequence.
	# 4. the mutation nucleotide in the library
	# note that for sites 6, 7, and 8 (and 9); they are on read 2 which is reverse-complement.

	# pertinent encoding of all ampliconversion details is in parse_FASTQ.py and here in myutils.format_genotype()
	if ampliconversion == 'short': #96A is the stops-down mutation (makes cag>Tag in forward notation)
		mutsites = [
			(1,22,'C','T'), # J numbering 2999, codon gcg(Ala) > gTg(Val)
			(2,54,'A','C'), # J numbering 3031, codon agc(Ser) > Cgc(Arg)
			(3,56,'C','A'), # J numbering 3033, codon agc(Ser) > agA(Arg) 
							# *** note G2 and G3 are in same codon; G2G3 == codon agc(Ser) > CgA(Arg)
			(4,57,'G','A'), # J numbering 3034
			(5,61,'T','C'), # J numbering 3038
			(6,101,'G','T'), # J numbering 3227; Forward notation C>A
			(7,99,'C','T'), # J numbering 3229; Forward notation G>A
			(8,98,'G','A') # J numbering 3230; Forward notation T>C
			]
	elif ampliconversion == 'long': # 164A is stops-down; 113A is a L>F
		mutsites = [
			(1,22,'C','T'), # J numbering 2999, codon gcg(Ala) > gTg(Val)
			(2,54,'A','C'), # J numbering 3031, codon agc(Ser) > Cgc(Arg)
			(3,56,'C','A'), # J numbering 3033, codon agc(Ser) > agA(Arg) 
							# *** note G2 and G3 are in same codon; G2G3 == codon agc(Ser) > CgA(Arg)
			(4,57,'G','A'), # J numbering 3034
			(5,61,'T','C'), # J numbering 3038

			(6,169,'G','T'), # J numbering 3227; Forward notation C>A
			(7,167,'C','T'), # J numbering 3229; Forward notation G>A
			(8,166,'G','A'), # J numbering 3230; Forward notation (actually C>T; this is an EvoC reversion)
			(9,86,'G','A') # J numbering 3310; Forward notation C>T
			]
	elif ampliconversion == 'longest':
		mutsites = [
			(1,22,'C','T'), # J numbering 2999, codon gcg(Ala) > gTg(Val)
			(2,54,'A','C'), # J numbering 3031, codon agc(Ser) > Cgc(Arg)
			(3,56,'C','A'), # J numbering 3033, codon agc(Ser) > agA(Arg) 
							# *** note G2 and G3 are in same codon; G2G3 == codon agc(Ser) > CgA(Arg)
			(4,57,'G','A'), # J numbering 3034
			(5,61,'T','C'), # J numbering 3038

			(6,239,'G','T'), # J numbering 3227; Forward notation C>A
			(7,237,'C','T'), # J numbering 3229; Forward notation G>A
			(8,236,'G','A'), # J numbering 3230; Forward notation (actually C>T; this is an EvoC reversion)
			(9,156,'G','A') # J numbering 3310; Forward notation C>T
			]

	# Note:
	# L-spec competitior = [10000001] (100000010) G1, G8
	# O-spec competitor = G1,3,6,8,9; [10100101]1 == (101001011)

	# it was easier to type out the mutsides above by hand, but a dict will be more useful in the code
	mutsites_dict = {} #keyed by concat read position, value is tup (ordinal site, wt, mut)
	for tup in mutsites:
		mutsites_dict[tup[1]] = (tup[0],tup[2],tup[3])

	# will assign 0, 1, or other to each of the 8 mutation sites to build the "binary" genotype
	binary_genotype = { 1:'',2:'',3:'',4:'',5:'',6:'',7:'',8:''}
	if ampliconversion == 'long' or ampliconversion == 'longest':
		binary_genotype[9] = ''
	# will tally additional mutations as a list of tuples of (position,mutant nucleotide)
	additional_muts = []

	assert len(genotype_seq) == len(refseq)

	# iterate the sequence by position and tally that position's info, 
	# either for the binary genotype over 8/9 mutation sites or as an additional site mutation.
	for pos in range(len(genotype_seq)): # note pos is 0-indexed

		if pos+1 in mutsites_dict: 
			#print("this is a mutsite position (1-indexing position %i), %s > %s" % (pos+1, refseq[pos],genotype_seq[pos]))
			#print("this is mutation site# %i" % mutsites_dict[pos+1][0])
			if genotype_seq[pos] == mutsites_dict[pos+1][1]:
				#print("this is WT")
				binary_genotype[mutsites_dict[pos+1][0]] = 0
			elif genotype_seq[pos] == mutsites_dict[pos+1][2]:
				#print("this is mut")
				binary_genotype[mutsites_dict[pos+1][0]] = 1
			else:
				#print("NOT WT OR MUT") # the binary genotype is not actually binary; it can also list a nucleotide for non-engineered mutations.
				binary_genotype[mutsites_dict[pos+1][0]] = genotype_seq[pos]
		else:
			if genotype_seq[pos] != refseq[pos]:
				additional_muts.append((pos+1,genotype_seq[pos]))

	if ampliconversion == 'short':
		binary_genotype_string = ''.join(str(binary_genotype[p]) for p in [1,2,3,4,5,6,7,8])
	elif ampliconversion == 'long' or ampliconversion == 'longest':
		binary_genotype_string = ''.join(str(binary_genotype[p]) for p in [1,2,3,4,5,6,7,8,9])
	additional_muts_string = '+'.join(str(m[0])+str(m[1]) for m in additional_muts)
	
	return(binary_genotype_string, additional_muts_string)

def construct_genotype_space(include_96A = False, length = 9):
	# Returns a list of all 2^length genotypes as strings of '[XXXXXXXX]';
	# optionally adds all 2^5=32 genotypes of the form [XXXXX000]+96A since these are all present in the lysogen libraries as incomplete recombinations
	genotypes = [i for i in itertools.product([0,1],repeat=length)]
	genotypes_strings = []
	for g in genotypes: # there has got to be a better way to do this
		g=[str(i) for i in g] 
		g=''.join(g)  
		genotypes_strings.append('[%s]' % g) 
		
	if include_96A:
		for g in genotypes:
			g=[str(i) for i in g] 
			g=''.join(g)  
			if g[5:8] == '000': 
				genotypes_strings.append('[%s]+96A' % g)
	return genotypes_strings

gt_to_L9index= {}
L9index_to_gt = {}
for (i,gt) in enumerate(construct_genotype_space(length=9)):
	gt_to_L9index[gt] = i
	L9index_to_gt[i] = gt

def hamming(s1,s2):
	# returns hamming distance between two equal-length strings.
	assert len(s1)==len(s2)
	h = 0
	for pos in range(len(s1)):
		if s1[pos] != s2[pos]:
			h += 1
	return h

def compute_diversity(df,metric='shannon'): # return array of length(generations) of diversity, for one simulation's data imported to df
	a_diversities = [] # array of diversities of length(generations)
	#print(df.loc[0, '[000000000]':'[111111111]'])  # syntax df.loc[row1:row2, col1:col2]
	if metric=='shannon':
		for t in df.index:
			total_t = df.loc[t, '[000000000]':'[111111111]'].sum() # total abundance at time t
			a_densities_t = df.loc[t, '[000000000]':'[111111111]'] / float(total_t) # array of genotype densities at time t
			diversity_t = sum([-i*np.log(i) if i>0 else 0 for i in a_densities_t]) # diversity at time t
			a_diversities.append(diversity_t)
	else:
		raise notavalidmetricerror
	return a_diversities

def grouper(n, iterable, fillvalue=None):
    # grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return itertools.zip_longest(fillvalue=fillvalue, *args)


# I/O:

def read_genotype_counts_file(file, as_dict=False, as_counts_dict=False):
	# Returns genotype counts/frequency data from a given file.
	# By default, as a list of tuples of (genotype, counts, frequency(%)); 
	# optionally, as a dict with key=genotype, value=frequency(%)
	instream = open(file,'r')
	data = {} if as_dict else []

	if as_dict:
		data = {} # values are frequencies for as_dict
	elif as_counts_dict:
		data = {} # values are counts for as_counts_dict
	else:
		data = [] # list of tups

	for line in instream:
		if 'genotype,counts,frequency(%)' in line:
			continue # discard the header
		(gt, counts, freq) = line.strip().split(',')
		if as_dict:
			data[str(gt)] = float(freq) 
		elif as_counts_dict:
			data[str(gt)] = int(counts) 
		else:
			data.append((str(gt), int(counts), float(freq)))
	return data

def read_selectionrate_file(file, as_dict=False):
	# Returns genotype selection rate from a given file.
	# the file is a csv with header containing "genotype,s_r"
	# By default, returns a list of tuples of (genotype, selection rate); 
	# optionally, as a dict with key=genotype, value=selection rate
	# This is used for multiple purposes besides reading selection rate files in the code.
	# obviously this can be used in the abstract to read any csv into (k,v) dict entries
	# or list of tupples of (a,b) but it is hardcoded to see a header containing 'genotype,*'

	instream = open(file,'r')
	data = {} if as_dict else []
	for line in instream:
		if 'genotype,' in line:
			continue # discard the header
		(gt, sr) = line.strip().split(',')
		if as_dict:
			data[str(gt)] = float(sr) 
		else:
			data.append((str(gt), float(sr)))
	return data

def makedir(path):
	if not os.path.isdir(path):
		os.mkdir(path)


# plotting/aesthetics:

def encode_genotype_color(gt, schema='default'):

	length = gt.find(']') - 1 # this assumes all gt are strings beginning with '['; length should result as 8 or 9`

	if length == 9:
		if schema == 'default':
			lib1_c = '#e9a3c9'  #lib1 programmed genotypes
			lib1_stop = '#c51b7d' #lib1 containing stop (but not LF)
			lib2_c = '#5ab4ac' 
			lib2_stop = '#01665e' 
			lib1_LF = '#a1d76a' 
			lib1_LF_stop = '#4d9221' 
			lib2_LF = '#fdae61' 
			lib2_LF_stop = '#f45043'
			other = '#bababa'
			if (len(gt) == 11) and not any([mut in gt for mut in ['A','T','G','C'] ] ):
				if gt[9] == '0':
					return lib1_c
				elif gt[9] == '1':
					return lib2_c
				else:
					return other
			elif '183A' in gt: # contains L>F mutation, which is 183A in read position numbering
				if gt[9] == '0':
					if '234A' in gt: # 234A is the stop codon in pre-mutagenesis template sequence
						return lib1_LF_stop
					else:
						return lib1_LF
				elif gt[9] == '1':
					if '234A' in gt: 
						return lib2_LF_stop
					else:
						return lib2_LF
				else:
					return other
			else: # does not contain L>F but has some nonprogrammed mutation
				if gt[9] == '0':
					if '234A' in gt: # stops
						return lib1_stop
					else:
						return other
				elif gt[9] == '1':
					if '234A' in gt: 
						return lib2_stop
					else:
						return other
				else:
					# appears to have a nonprogrammed mutation at position 9
					return other

		elif schema == 'stops_only':
			stops_c = '#e41a1c' # 
			other_c = '#f4cae4' #
			if '234A' in gt:
				return stops_c
			else:
				return other_c

		elif schema == 'lib1vlib2':
			lib1_c = '#FFCB05' # all lib1 programmed b=0/1
			lib2_c = '#00274C' # all lib2 programmed b=0/1
			other_c = '#bababa' # Ethan's pick. for any gt besides [bbbbbbbbX]
			if (len(gt) == 11) and not any([mut in gt for mut in ['A','T','G','C'] ] ):
				if gt[9] == '0': 
					return lib1_c
				elif gt[9] == '1': 
					return lib2_c
				else: return other_c
			else:
				return other_c

		elif schema == 'prog-stop-LF':
			stops_c = '#d95f02' # stopcodon color overrides anything below
			prog_c = '#7570b3' # only the 512 programmed genotypes
			LF_c = '#1b9e77' # after stopcodon and programmed, any L>F containing genotype
			other_c = '#dadada' # Ethan suggested this one
			if '234A' in gt:
				return stops_c
			elif (len(gt) == 11) and not any([mut in gt for mut in ['A','T','G','C'] ] ):
				return prog_c
			elif '183A' in gt:
				return LF_c
			else:
				return other_c

		elif schema == 'black':
			return 'k'

	elif length == 8:

		if schema == 'prog-stop-LF':
			stop_color = '#e41a1c' # stops color overrides anything below
			binarygt_color = '#7570b3' # only the 256 programmed
			other_color = '#dadada' # Ethan picked this one as well.
			stop_color = '#d95f02'
			binarygt_color = '#7570b3'
		else:
			stop_color = 'tab:red' # for all [XXXXX000]+96A
			binarygt_color = 'tab:blue'
			other_color = 'tab:orange'

		if (len(gt) == 10) and (not any([mut in gt for mut in ['A','T','G','C'] ] ) ):
			return binarygt_color
		elif (gt[6:9] == '000') and (']+96A' in gt[-5:]):
			return stop_color
		else:
			return other_color


# NETWORK FUNCTIONS

def BuildHypercubeGraph(length=8):
	G = nx.Graph()
	genotype_space = construct_genotype_space(length=length)
	# add nodes
	for gt in genotype_space:
		G.add_node(gt)
	# add edges
	for gt in genotype_space:
		for gt2 in genotype_space:
			if hamming(gt,gt2) == 1:
				G.add_edge(gt,gt2)
	return G

def BuildSubGraph(genotypes):
	G = nx.Graph()
	# add nodes
	for gt in genotypes:
		G.add_node(gt)
	# add edges
	for gt in genotypes:
		for gt2 in genotypes:
			if hamming(gt,gt2) == 1:
				G.add_edge(gt,gt2)
	return G

def position_bysr(G,srdata,jitter):
	# build positions (x,y) which are x=num muts, y= selection rate
	pos = {}
	for node in G.nodes:
		# how many mutations from [00000000] is this node?
		num_muts = sum( int(i) for i in node[1:-1] )
		if not jitter:
			pos[node] = (num_muts, srdata[node])
		else:
			if num_muts == 0:
				pos[node] = (num_muts, srdata[node])
			else:
				xpos = num_muts + np.random.uniform(low=-jitter,high=jitter)
				pos[node] = (xpos, srdata[node])
	return pos

def position_stacked(G):
	# build positions (x,y) which are x = num muts, y= stacked (not tied to quantity)
	pos = {}
	counter = {} # keyed by nummuts below: value is used to build y-position
	for node in G.nodes:
		# how many mutations from [00000000] is this node?
		num_muts = sum( int(i) for i in node[1:-1] )
		if num_muts not in counter:
			counter[num_muts] = 0
		else:
			counter[num_muts] += 1
		pos[node] = (num_muts, 4*counter[num_muts])
	return pos

def position_stacked_tight(G):
	# build positions (x,y) per node:
		# x = hamming, with some displacement scaled by adjacent_offset, 
		# y = stacking nodes scaled by node_spacing
	pos = {}
	counter = {} # keyed by nummuts below: value is used to build y-position
	node_spacing = 15
	adjacent_offset = 0.15

	for node in G.nodes:
		hamming = sum( int(i) for i in node[1:-1] )
		if hamming not in counter: counter[hamming] = 0 
		else: counter[hamming] += 1

		if hamming == 0 or hamming == 9:
			pos[node] = (hamming,0)
		elif hamming == 1 or hamming == 8:
			pos[node] = (hamming, node_spacing*counter[hamming]-node_spacing*9/2.)
		elif hamming == 2 or hamming == 7:
			if counter[hamming] <= 17:
				pos[node] = (hamming - adjacent_offset/2., node_spacing*counter[hamming]-node_spacing*18/2.)
			else:
				pos[node] = (hamming + adjacent_offset/2., node_spacing*(counter[hamming]-18)-node_spacing*18/2.)
		elif hamming == 3 or hamming == 6:
			if counter[hamming] <= 27:
				pos[node] = (hamming - adjacent_offset, node_spacing*counter[hamming]-node_spacing*28/2.)
			elif counter[hamming] <= 27+28:
				pos[node] = (hamming, node_spacing*(counter[hamming]-28)-node_spacing*28/2.)
			else:
				pos[node] = (hamming + adjacent_offset, node_spacing*(counter[hamming]-28-28)-node_spacing*28/2.)
		elif hamming == 4 or hamming == 5:
			if counter[hamming] <= 41:
				pos[node] = (hamming - adjacent_offset, node_spacing*counter[hamming]-node_spacing*42/2.)
			elif counter[hamming] <= 41+42:
				pos[node] = (hamming, node_spacing*(counter[hamming]-42)-node_spacing*42/2.)
			else:
				pos[node] = (hamming + adjacent_offset, node_spacing*(counter[hamming]-42-42)-node_spacing*42/2.)
	return pos
