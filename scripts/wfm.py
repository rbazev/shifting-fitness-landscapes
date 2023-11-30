# perform simulations of lambda evolution using a modified wright-fisher model incorporating experimentally-measured fitness landscapes 

import myutils
import numpy as np
from numpy.random import binomial
from numpy.random import multinomial
from numpy.random import randint
from numpy.random import normal
import pandas as pd
import argparse
import time

np.random.seed(1)

def load_sr_landscape_to_w(file, nan_conversion=-10, noise=0.0):
	'''
	file provides path to the selection rate file. it should be in 9-site
	genotype format as produced by merge_srdata.py in the correct order
	of genotypes given by myutils.construct_genotype_space(length=9).

	nan_conversion specifies the selection rate assigned to NaN values.

	noise sets the scaling parameter (st dev of normal dist) of normally-distributed
	noise added to the selection rates before transforming into w.

	returns a numpy array of exponentiated selection rates in genotype order.
	'''
	srdf = pd.read_csv(file, index_col=0, nrows=512) # only reads the 512 programmed genotypes in case a file has additional, unprogrammed genotypes

	# add normally-distributed noise to the experimentally measured selection rates;
	# the scale of this noise is tuned such that it produces a correlation with the unnoised landscape
	# similar to the correlation between replicate experiments. 
	noise_array = normal(size=512,loc=0, scale=noise)
	srdf['s_r_merged'] = srdf['s_r_merged'] + noise_array

	# the experimental selection rate is per 4 hours; 
	# but here we want to model per 1.5 hours (single lambda infection cycle) 
	# using a scaled selection rate as below
	srdf /= 2.6 
	wdf = srdf.fillna(value=nan_conversion).apply(np.exp, axis=1).set_axis(['w'], axis=1)
	assert list(wdf.index) == myutils.construct_genotype_space(length=9)
	return wdf.to_numpy().flatten()

def draw_noised_ws_dict(sr_files_dir, noisepersim):
	print('drawing a new set of landscapes with noise-per-simulation %s' % noisepersim)
	ws_dict = {}
	# have data from 5 environments
	# these each use merged data over 6 flasks (lib A x 3 flasks, lib B x 3 flasks)
	ws_dict['vEvoC-mixed'] = load_sr_landscape_to_w(sr_files_dir + 'sr_LibAvg_vGen.csv', noise=noisepersim)
	ws_dict['vLspec-mixed'] = load_sr_landscape_to_w(sr_files_dir + 'sr_LibAvg_vLspec.csv', noise=noisepersim)
	ws_dict['vOspec-mixed'] = load_sr_landscape_to_w(sr_files_dir + 'sr_LibAvg_vOspec.csv', noise=noisepersim)
	ws_dict['vEvoC-onlyOmpF'] = load_sr_landscape_to_w(sr_files_dir + 'sr_LibAvg_vGen_OmpF.csv', noise=noisepersim)
	ws_dict['vEvoC-onlyLamB'] = load_sr_landscape_to_w(sr_files_dir + 'sr_LibAvg_vGen_LamB.csv', noise=noisepersim)
	return ws_dict

def initialize_population(startingpopulation,popsize):
	'''
	startingpopulation is a path to csv file containing the starting population. 
	the csv should be in the format genotype,density.
	density is computed based on the total across all lines, so the density column in the csv does not need to sum to 1.
	returns a flattned numpy array of abundances in genotype order based on computed density from the file and the specified popsize.
	'''	

	if startingpopulation:
		density_dict = myutils.read_selectionrate_file(startingpopulation, as_dict=True)
		total = sum([density_dict[gt] for gt in density_dict])
	else:
		#default starting population is all wildtype/EvoC
		density_dict = {'[000000000]':1.0}
		total = sum([density_dict[gt] for gt in density_dict])

	# Initialize the starting population 
	abundance = pd.DataFrame(index=myutils.construct_genotype_space(length=9))
	abundance['0'] = 0.0 # initialize column '0' with all 0.0's
	# assign population to genotypes:
	for gt in density_dict:
		abundance.at[gt,'0'] = int(float(popsize)*float(density_dict[gt])/float(total))

	# now that the starting population has been initialized as a genotype-ordered dataframe,
	# the rest of the simulation will be indexless in numpy arrays for optimization
	abundance = abundance.to_numpy().flatten()
	return abundance

def compute_average_specialization_index(abundance,si_d):
	'''
	abundance is a numpy array in genotype order, giving the number of particles for each genotype
	si_d is a dictionary keyed by genotype of specialization index

	returns a weighted average of specialization index across the population; 
		eg returns sum_gt(abundance_gt*si_gt)/512
	'''
	weighted_sum = 0.0
	abundance_sum = 0.0
	for i,gt in enumerate(myutils.construct_genotype_space()):
		#print("on genotype #%s:%s; abundance is %s and SI is %s" % (i,gt, abundance[i],si_d[gt]))
		weighted_sum += abundance[i]*si_d[gt]
		abundance_sum += abundance[i]
	return weighted_sum/abundance_sum


def compute_pvals(abundance, w, noise=False):
	'''
	abundance and w are both numpy arrays (in genotype order) 
	noise scales the added noise to w representing uncertainties in experimentally measured fitness

	the return value is a numpy array (in genotype order) of success probability to be observed in the progeny generation; 
	the probability for genotype i = (abundance of i)*(w of i)/sum(abundance*w for all genotypes)
	'''
	if noise:
		# add normally-distributed noise to the selection rate; note this is noise-per-generation (not used in the analysis in the paper.)
		noise_array = normal(size=512,loc=0, scale=noise) # the noise scale had been tuned to experiments measuring sr over 4 hours.
		noise_array /= 2.6 # this scales the noise to the per-generation sr used in the simulations.
		w = np.exp(np.log(w)+noise_array)

	nw = np.multiply(abundance,w)
	nwsum = float(nw.sum())
	pvals = nw/nwsum
	return pvals.flatten()


def reproduce(pvals, population_size):
	# pvals is numpy array in genotype order.
	# population_size is the number of draws from the multinomial distribution to populate the progency population
	# returns numpy array of abundance. it is implicitly in genotype order since pvals is assumed to be in genotype order.
	progeny = multinomial(population_size,pvals)
	return progeny


def mutate_genotype(genotype):
	# for a given genotype, introduces a random mutation (1>0 or 0>1 at a random site)
	# genotype is string in form '[101010101]'; returns the mutated genotype.
	mut_site = randint(1,10) # returns random integer between 1 and 9 inclusive.
	mutation = '1' if genotype[mut_site] == '0' else '0'
	return genotype[:mut_site] + mutation + genotype[mut_site + 1:]

def mutate_population(abundance,mutation_rate=7.7e-8):
	# abundance is a pandas series indexed by genotype providing the abundances of genotypes after a replication
	# cycle [derived by drawing on a multinomial distribution by reproduce()]; 
	# mutation rate is the per-bp per-replication cycle estimate of mutation rate in lambda (7.7e-8 per bp per replication)
	# 
	# mutrate is divided by three since the mutations at each site are one of the three possible nucleotide mutations.
	# this is multiplied by 9 to grossly estimate the mutation rate over 9 sites
	# then, the probability of one mutation (ignoring multiple mutations) by poisson distribution is
	# P(1 mutation) = mutrate * exp(-mutrate) which is essentially the mutation rate
	#
	# For each genotype, samples are drawn from a binomial distribution (n trials = progeny population for that genotype),
	# probability of success p = mutrate
	#
	# the return value is an array of abundances, based on the starting abundance array, 
	# where individuals are moved at some rate based on the mutation rate from their parental genotype to a neighboring genotype. 

	mu_c = mutation_rate*9.0/3.0 # effective mu based on 9 sites possible and 1 of 3 nucleotide mutations possible
	mutated_from = binomial(abundance, mu_c)
	abundance -= mutated_from
	
	for (index,num_mut_from) in enumerate(mutated_from):
		gt = myutils.L9index_to_gt[index]
		for i in range(num_mut_from):
			mutated_to_genotype = mutate_genotype(gt)
			mutated_to_index = myutils.gt_to_L9index[mutated_to_genotype]
			abundance[mutated_to_index] += 1

	return abundance






# parse command-line arguments:
parser = argparse.ArgumentParser(description='simulate evolution on fitness landscapes')
parser.add_argument('outputdir', type=str, help='main output directory, should contain ./genotypecounts and ./selectionrate directories.')
parser.add_argument('simsetname', type=str, help='a simulation set name, will have multiple runname subdirectories within it')
parser.add_argument('runname', type=str, help='short label for this simulation run used for subdirectory organization etc')
parser.add_argument('--numgenerations', type=int, default=1000, help='number of generations to simulate')
parser.add_argument('--numsims', type=int, default=5, help='number of independent simulations')
parser.add_argument('--landscapetrajectory', type=str, default='continuous', help='method of determining which fitness landscape(s) are used')
parser.add_argument('--noisepergen', type=float, default=0.0, help='guassian noise added to selection rates PER GENERATION')
parser.add_argument('--noisepersim', type=float, default=0.8, help='guassian noise added to selection rates PER SIMULATION')
parser.add_argument('--popsize', type=int, default=5e9, help='fixed population size')
parser.add_argument('--mutrate', type=float, default=7.7e-8, help='mutation rate per bp per replication cycle')
parser.add_argument('--startingpopulation', type=str, default=None, help='path to csv listing abundances for starting population genotypes')
args = parser.parse_args()

myutils.makedir(args.outputdir + '/simulations/')

sim_set_subdir = args.outputdir + '/simulations/%s/' % args.simsetname
myutils.makedir(sim_set_subdir)
sim_subdir = sim_set_subdir + '/%s/' % args.runname
myutils.makedir(sim_subdir)
simdata_subdir = sim_subdir + '/simdata/'
myutils.makedir(simdata_subdir)

# Landscape(s): load sr landscapes into w arrays. 
sr_files_dir = args.outputdir + '/selectionrate/srfiles/'

# import L and O specialist genotype identities from prior analysis
specialist_files_dir = args.outputdir + '/selectionrate/specialists/'
Lspecs = myutils.read_selectionrate_file(specialist_files_dir  + 'Lspecialists.csv', as_dict=True)
Ospecs = myutils.read_selectionrate_file(specialist_files_dir + 'Ospecialists.csv', as_dict=True)

# import specialization index
si_file = specialist_files_dir + 'si.csv'
si_dict = myutils.read_selectionrate_file(si_file, as_dict=True)

# Simulation parameters now taken from argparse:
population_size = args.popsize
mutation_rate = args.mutrate
num_generations = args.numgenerations
num_sims = args.numsims
landscape_trajectory = args.landscapetrajectory
	# string describing landscape trajectory in the form:
	# 'static:x' where x is a string key for ws_dict as above; or 
	# 'continuous' for continuous weighted-average transitions across 5 landscapes 
	# 'discrete' for discrete jumps between experimental landscapes, using the nearest landscape based on population SI
noisepergen = args.noisepergen # scaling factor for introducing normally-distributed noise on the experimental selection rates (ie, stdev of normally-distributed, 0-centered noise)
noisepersim = args.noisepersim

logout = open(sim_subdir + '/log.txt', 'w')
logout.write(
	'population_size=%s\n' % population_size +
	'startingpop=%s\n' % args.startingpopulation +
	'mutation_rate=%s\n' % mutation_rate +
	'num_generations=%s\n' % num_generations +
	'num_sims=%s\n' % num_sims +
	'landscape_trajectory=%s\n' % landscape_trajectory +
	'noisepergen=%s\n' % noisepergen +
	'noisepersim=%s\n' % noisepersim
	)

start_time = time.time()

for sim_i in range(1,num_sims+1):

	print("===================================== beginning simulation %s" % sim_i)

	# re-draw selection rate data each simulation with noise-per-simulation 
	ws_dict = draw_noised_ws_dict(sr_files_dir, args.noisepersim)

	# Initialize the starting population 
	abundance = initialize_population(args.startingpopulation,args.popsize)
	
	# for bookkeeping: 
	abundance_record = {} # dict of abundance nparrays keyed by generation #
	abundance_record['0'] = abundance # the initialized abundances at time 0.
	landscape_record = {'0':['NA']} # dict of landscapes keyed by generation #

	for gen_i in range(1,num_generations+1):

		if gen_i % 20 ==0: print("Beginning generation %s" % gen_i)

		# 1. determine the landscape to be used to generate probabilities of reproductive success
		if landscape_trajectory == 'continuous':
			
			# compute the average specindex (weighted by genotype abundance):
			weighted_si = compute_average_specialization_index(abundance,si_dict)

			Lspec_competitor_si = si_dict['[100000010]']
			Ospec_competitor_si = si_dict['[101001011]']
			# use the si of the competitors, which were at 90% at the start of the competition, 
			# to place the landscapes on the "net si" axis from -1 to 1.
			Lspec_competitor_coord = Lspec_competitor_si*0.9
			Ospec_competitor_coord = Ospec_competitor_si*0.9
			assert weighted_si >= -1.0 and weighted_si <= 1.0
			
			if weighted_si == -1.0:
				w = ws_dict['vEvoC-onlyLamB']
				landscape_record['%s' % gen_i] = ['si= %.2f; 1.0 vEvoC-onlyLamB' % weighted_si]

			elif weighted_si > -1.0 and weighted_si < Ospec_competitor_coord:
				wkey1 = 'vEvoC-onlyLamB'
				wkey2 = 'vOspec-mixed'
				m = np.abs(weighted_si - (Ospec_competitor_coord))
				l = np.abs(-1.0 - weighted_si)
				weight1 = m/(l+m)
				weight2 = 1.0 - weight1
				w = weight1*ws_dict[wkey1] + weight2*ws_dict[wkey2]
				landscape_record['%s' % gen_i] = ['si= %.2f; %.2f %s + %.2f %s' % (weighted_si, weight1, wkey1, weight2, wkey2)]

			elif weighted_si == Ospec_competitor_coord:
				w = ws_dict['vOspec-mixed']
				landscape_record['%s' % gen_i] = ['si= %.2f; 1.0 vOspec-mixed' % weighted_si]

			elif weighted_si > Ospec_competitor_coord and weighted_si < 0.0:
				wkey1 = 'vOspec-mixed'
				wkey2 = 'vEvoC-mixed'
				m = np.abs(weighted_si - (0.0))
				l = np.abs(Ospec_competitor_coord - weighted_si)
				weight1 = m/(l+m)
				weight2 = 1.0 - weight1
				w = weight1*ws_dict[wkey1] + weight2*ws_dict[wkey2]
				landscape_record['%s' % gen_i] = ['si= %.2f; %.2f %s + %.2f %s' % (weighted_si, weight1, wkey1, weight2, wkey2)]

			elif weighted_si == 0.0:
				w = ws_dict['vEvoC-mixed']
				landscape_record['%s' % gen_i] = ['si= %.2f; 1.0 vEvoC-mixed' % weighted_si]

			elif weighted_si > 0.0 and weighted_si < Lspec_competitor_coord:
				wkey1 = 'vEvoC-mixed'
				wkey2 = 'vLspec-mixed'
				m = np.abs(Lspec_competitor_coord - weighted_si)
				l = np.abs(weighted_si - 0.0)
				weight1 = m/(l+m)
				weight2 = 1.0 - weight1
				w = weight1*ws_dict[wkey1] + weight2*ws_dict[wkey2]
				landscape_record['%s' % gen_i] = ['si= %.2f; %.2f %s + %.2f %s' % (weighted_si, weight1, wkey1, weight2, wkey2)]

			elif weighted_si == Lspec_competitor_coord:
				w = ws_dict['vLspec-mixed']
				landscape_record['%s' % gen_i] = ['si= %.2f; 1.0 vLspec-mixed' % weighted_si]

			elif weighted_si > Lspec_competitor_coord and weighted_si < 1.0:
				wkey1 = 'vLspec-mixed'
				wkey2 = 'vEvoC-onlyOmpF'
				m = np.abs(1.0 - weighted_si)
				l = np.abs(weighted_si - Lspec_competitor_coord)
				weight1 = m/(l+m)
				weight2 = 1.0 - weight1
				w = weight1*ws_dict[wkey1] + weight2*ws_dict[wkey2]
				landscape_record['%s' % gen_i] = ['si= %.2f; %.2f %s + %.2f %s' % (weighted_si, weight1, wkey1, weight2, wkey2)]

			elif weighted_si == 1.0:
				w = ws_dict['vEvoC-onlyOmpF']
				landscape_record['%s' % gen_i] = ['si= %.2f; 1.0 vEvoC-onlyOmpF' % weighted_si]

			else:
				raise landscapeaxiserror

		elif landscape_trajectory.count('static') == 1:
			wkey = landscape_trajectory[7:] # 'static:xxxx' > 'xxxx'
			w = ws_dict[wkey]

			weighted_si = compute_average_specialization_index(abundance,si_dict) # only calculating for bookkeeping; not used to select landscape
			landscape_record['%s' % gen_i] = ['si= %.2f; static:%s' % (weighted_si, wkey)]


		elif landscape_trajectory == 'discrete':
			# use the experimental fitness landscape closest to the simulated population's specialization index, `weighted_si`.
			weighted_si = compute_average_specialization_index(abundance,si_dict)

			Lspec_competitor_si = si_dict['[100000010]']
			Ospec_competitor_si = si_dict['[101001011]']
			# use the si of the competitors, which were at 90% at the start of the competition, 
			# to place the landscapes on the "net si" axis from -1 to 1. (same math as in the continuous shifting model above)
			Lspec_competitor_coord = Lspec_competitor_si*0.9
			Ospec_competitor_coord = Ospec_competitor_si*0.9
			assert weighted_si >= -1.0 and weighted_si <= 1.0

			# define the 5 experimental landscapes by their respective SI coordinates.
			wkeys_dict = { -1.0 : 'vEvoC-onlyLamB',
						   Ospec_competitor_coord : 'vOspec-mixed',
						   0.0 : 'vEvoC-mixed',
						   Lspec_competitor_coord : 'vLspec-mixed',
						   1.0 : 'vEvoC-onlyOmpF'}

			# find the discrte SI coordinate corresponding to an experimental landscape above, 
			# which is closest to the simulation's SI (`weighted_si`):
			coords = [-1.0, Ospec_competitor_coord, 0.0, Lspec_competitor_coord, 1.0]
			discrete_coord = min(coords, key=lambda x:abs(x - weighted_si)) 

			# use the experimental fitness landscape nearest to the weighted_si, without interpolation between two landscapes
			wkey = wkeys_dict[discrete_coord]
			w = ws_dict[wkey] 

			# keep record of the weighted_si this generation, and the discrete fitness landscape used:
			landscape_record['%s' % gen_i] = ['si= %.2f; landscape=%s' % (weighted_si, wkey)]

		else:
			raise landscapetrajectoryerror

		# 2. compute probabilities of reproductive success
		pvals = compute_pvals(abundance, w, noisepergen)

		# 3. draw progeny using given probabilities
		# draw progeny from multinomial distribution with the given probability values array
		# for simplicity I have separated replication and mutation functions; this is just replication:
		time_i_progeny = reproduce(pvals,population_size) 

		# 4. introduce mutations into the replicated genomes 
		# in the course of replication mutations will accumulate so we introduce those mutations now:
		# mutate time_i_progeny and insert as ultimate column in abundances with column name step_i
		new_abundance = mutate_population(time_i_progeny, mutation_rate) 

		# 5. bookkeeping end of generation and reset abundance array for next generation
		abundance_record['%s' % gen_i] = new_abundance
		abundance = new_abundance

	
	abundances = pd.DataFrame.from_dict(abundance_record)
	# now transpose into plot-friendly format and save to CSV
	abundances = abundances.T
	abundances.columns = myutils.construct_genotype_space(length = 9)
	abundances.insert(loc = 0,
          column = 'Time',
          value = abundances.index)

	landscapes = pd.DataFrame.from_dict(landscape_record)
	landscapes = landscapes.T
	landscapes.columns = ['landscape']
	landscapes.insert(loc=0, column='Time', value=landscapes.index)
	# could consider saving this as a separate csv if desired, but i'll just append it to the abundances to save csv with all data
	abundances = abundances.join(landscapes['landscape'])

	abundances.to_csv('%s/sim-%s.csv' % (simdata_subdir,sim_i), index=False)

logout.write('execution time = %s' % (time.time() - start_time))
logout.close()
