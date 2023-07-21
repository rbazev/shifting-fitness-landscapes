import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LassoLarsIC
from sklearn.linear_model import LassoCV
from sklearn.pipeline import make_pipeline
import statsmodels.api as sm
import myutils
import itertools
import argparse

parser = argparse.ArgumentParser(description='build linear regression models and perform LASSO regression for feature selection')
parser.add_argument('outputdir', type=str, help='path to top-level output directory') # from here, find sr files in ./selectionrate/srfiles/*
parser.add_argument('runname', type=str, help='nickname for this run of mlr.py, results will go into (outputdir)/selectionrate/regression/(runname)/')
parser.add_argument('sampleids',type=str, help='list of sample IDs to include', nargs='+')
args = parser.parse_args()
srfiles_dir = args.outputdir + '/selectionrate/srfiles/'

myutils.makedir(args.outputdir + '/selectionrate/regression/')
this_run_outputdir = args.outputdir + '/selectionrate/regression/%s/' % args.runname
myutils.makedir(this_run_outputdir)

logfile = open('%s/%s-logfile.txt' % (this_run_outputdir,args.runname), 'w')
logfile.write("Doing linear regression and model selection within the main data output directory:\n%s\n" % args.outputdir)
logfile.write("Will include the following samples:\n%s\n" % args.sampleids)

# will use these to build and access columns of the feature matrix for regression:
# I have reduced G2 and G3 features to a single feature G23R representing either G2, or G3, or the combination;
# as all of these result in the same codon mutation to arginine and appear to have consistent selection rates
G_terms = ['G1', 'G23R', 'G4', 'G5', 'G6', 'G7', 'G8', 'G9']
GxG_terms_unjoined = list(itertools.combinations(G_terms, 2))
GxG_terms = []
for (t1,t2) in GxG_terms_unjoined:
	GxG_terms.append('%sx%s' % (t1,t2))
C_terms = ["CE", "CL", "CO"] # competitor is EvoC, L-specialist, or O-specialist
CxG_terms = ["%sx%s" % (c,g) for c in C_terms for g in G_terms]
CxGxG_terms = ["%sx%s" % (c,g) for c in C_terms for g in GxG_terms]

# build the dictionary with keys for each column that will eventually go into a pandas dataframe
feature_matrix = {'genotype vs Cx' : []} 
for gx in G_terms: feature_matrix[gx] = []
for gxgx in GxG_terms: feature_matrix[gxgx] = []
for cx in C_terms: feature_matrix[cx] = []
for cxgx in CxG_terms: feature_matrix[cxgx] = []
for cxgxgx in CxGxG_terms: feature_matrix[cxgxgx] = []
feature_matrix['fitness'] = []

for sample in args.sampleids:

	logfile.write( "\n\nParsing sr data from %s to include in feature matrix" % sample)
	experiment_factors = myutils.competition_factors(sample)
	lib = experiment_factors['library']
	comp = experiment_factors['competitor']
	host = experiment_factors['host']
	logfile.write("For %s, Library %s vs. Competitor %s on Host %s" % (sample,lib,comp,host))

	# selection rate data to incorporate into feature matrix
	sr_data = myutils.read_selectionrate_file("%s/sr_%s.csv" % (srfiles_dir, sample), as_dict=True)
	logfile.write( "read selection rates from %s/sr_%s.csv" % (srfiles_dir, sample))
	
	# iterate through the selection rate data for this experiment;
	# add the appropriate data to every column of the feature_matrix 
	# (ie, append to the list contained at every key, as it is currently a dict and later will be a df)
	for gt in sr_data:
		if gt not in myutils.construct_genotype_space():
			continue
		feature_matrix['genotype vs Cx'].append(gt+'vs.'+comp[0])
		# now need to parse the genotype and add 0 or 1 to all of the indicator variables Gx, GxGx, Cx, CxGx, CxGxGx..
		stripped = gt.strip('[]')
		mut_positions = [int(pos)+1 for pos, char in enumerate(stripped) if char == '1']

		for gx in G_terms: # G1, G23R, G3, G4, G5, G6...
			if gx == 'G23R':
				feature_matrix[gx].append(1) if (2 in mut_positions) or (3 in mut_positions) else feature_matrix[gx].append(0)
			else:
				feature_matrix[gx].append(1) if int(gx[1]) in mut_positions else feature_matrix[gx].append(0)
		
		for gxgx in GxG_terms: # G1xG23R, G1xG4, G1xG5, ... 
			if 'G23R' in gxgx:
				remaining_term = gxgx.strip('G23R').strip('x').strip('G') # what remains is just an integer specifying the other G-position
				#input('detected g23r in this gxgx and stripped remainder: %s, %s' % (gxgx,remaining_term))
				feature_matrix[gxgx].append(1) if int(remaining_term) in mut_positions and (2 in mut_positions or 3 in mut_positions) else feature_matrix[gxgx].append(0)
			else:
				#input('did not detect g23r in this gxgx: %s' % gxgx)
				feature_matrix[gxgx].append(1) if int(gxgx[1]) in mut_positions and int(gxgx[4]) in mut_positions else feature_matrix[gxgx].append(0)
		
		for cx in C_terms:
			# comp[0] gives either E, L, or O to designate competitor
			# cx is either "CE", "CL", "CO"
			feature_matrix[cx].append(1) if comp[0] in cx[1] else feature_matrix[cx].append(0)

		for cxgx in CxG_terms:
			if 'G23R' in cxgx:
				if comp[0] in cxgx[1] and (2 in mut_positions or 3 in mut_positions):
					feature_matrix[cxgx].append(1)
				else:
					feature_matrix[cxgx].append(0)
			else:
				if comp[0] in cxgx[1] and int(cxgx[4]) in mut_positions:
					feature_matrix[cxgx].append(1)
				else:
					feature_matrix[cxgx].append(0)

		for cxgxgx in CxGxG_terms:
			if 'G23R' in cxgxgx:
				if cxgxgx[4] == '2': # G23R is the first G, so the second G is the end of the string
					remaining_mutsite = int(cxgxgx[-1])
					#print("For cxgxgx term %s; G23R is the first G; found remaining_mutsite as: %s" % (cxgxgx, remaining_mutsite))
				else: # G23R is the second G, so the first G site is at index 4
					remaining_mutsite = int(cxgxgx[4])
					#print("For cxgxgx term %s; G23R is the second; found remaining_mutsite as: %s" % (cxgxgx, remaining_mutsite))
				if comp[0] in cxgxgx[1] and (2 in mut_positions or 3 in mut_positions) and remaining_mutsite in mut_positions:
					feature_matrix[cxgxgx].append(1)
				else:
					feature_matrix[cxgxgx].append(0)
			else:
				if comp[0] in cxgxgx[1] and int(cxgxgx[4]) in mut_positions and int(cxgxgx[7]) in mut_positions:
					feature_matrix[cxgxgx].append(1)
				else:
					feature_matrix[cxgxgx].append(0)

		feature_matrix['fitness'].append(sr_data[gt])


# wrangling the full feature matrix X and the fitnesses y, together they are within the dataframe df created from feature_matrix
df = pd.DataFrame(feature_matrix)
df.to_csv('%s/%s-fullmatrix.csv' % (this_run_outputdir,args.runname))
df = df.dropna() # drop rows with NaN, unobserved fitness for that genotype against that competitor in that replicate of the experiment
logfile.write("\n\nThere are %s rows (observed fitness measurements) and %s columns (terms of model, including a fitness column)\n\n" % df.shape)
# create feature matrix X with just a column for constant (intercept) and interaction terms
X = df.drop('genotype vs Cx', axis=1).drop('fitness',axis=1)
X = sm.add_constant(X)
X.to_csv('%s/%s-Xmatrix-with-const.csv' % (this_run_outputdir,args.runname))
y = df['fitness']
y.to_csv('%s/%s-Ymatrix.csv' % (this_run_outputdir,args.runname))

# bookkeeping dataframe resultsdf used later:
allcoefs = list(X) # list of coefficient terms for the regression
resultsdf = pd.DataFrame(allcoefs, index=allcoefs) # dataframe with rows for each coefficient in the model, will store results here
resultsdf['human-readable'] = np.arange(1,resultsdf.shape[0]+1) # number the terms in order for sorting later if needed
resultsdf.drop(resultsdf.columns[[0]], axis=1, inplace=True)

# see "Selecting Lasso via an information criterion" on https://scikit-learn.org/stable/auto_examples/linear_model/plot_lasso_model_selection.html#sphx-glr-auto-examples-linear-model-plot-lasso-model-selection-py
# data accession see https://stackoverflow.com/questions/58615904/how-to-extract-coefficients-from-fitted-pipeline-for-penalized-logistic-regressi
logfile.write("\n\ntuning alpha for AIC...\n\n")
lasso_lars_ic = make_pipeline(StandardScaler(), LassoLarsIC(criterion="aic")).fit(X, y) 
results = pd.DataFrame(
    {
        "alphas": lasso_lars_ic[-1].alphas_,
        "AIC criterion": lasso_lars_ic[-1].criterion_,
    }
).set_index("alphas")
alpha_aic = lasso_lars_ic[-1].alpha_
logfile.write("\n\nalpha chosen by tuning to AIC: %s\n\n" % alpha_aic)
aic_tuned_coefs = lasso_lars_ic['lassolarsic'].coef_

# repeat using BIC
logfile.write("tuning alpha for BIC...")
lasso_lars_ic.set_params(lassolarsic__criterion="bic").fit(X, y)
results["BIC criterion"] = lasso_lars_ic[-1].criterion_
alpha_bic = lasso_lars_ic[-1].alpha_
logfile.write("\n\nalpha chosen by tuning to BIC: %s\n\n" % alpha_bic)
bic_tuned_coefs = lasso_lars_ic['lassolarsic'].coef_

# plot tuning alpha to min AIC and BIC:
ax = results.plot()
ax.vlines(
    alpha_aic,
    results["AIC criterion"].min(),
    results["AIC criterion"].max(),
    label="alpha: AIC estimate",
    linestyles="--",
    color="tab:blue",
)
ax.vlines(
    alpha_bic,
    results["BIC criterion"].min(),
    results["BIC criterion"].max(),
    label="alpha: BIC estimate",
    linestyle="--",
    color="tab:orange",
)
ax.set_xlabel(r"$\alpha$")
ax.set_ylabel("criterion")
ax.set_xscale("log")
ax.legend()
ax.set_title("Information-criterion for model selection: %s" % args.runname)
plt.savefig("%s/%s-modelselection-AICBICtuning.pdf" % (this_run_outputdir,args.runname))
plt.close()

# plot coefficients in lasso models tuned with AIC and BIC
# X.columns contains the feature names
width = 0.35       # the width of the bars
rects1 = plt.bar(X.columns, aic_tuned_coefs, -width, align='edge',
                color='tab:blue', label='AIC-tuned')
rects2 = plt.bar(X.columns, bic_tuned_coefs, +width, align='edge',
                color='tab:orange', label='BIC-tuned')
plt.xlabel('Feature')
plt.ylabel('Coefficient')
plt.xticks(rotation=90, fontsize=3)
plt.title('Tuned model by AIC or BIC; %s' % args.runname)
plt.legend()
plt.savefig("%s/%s-coefs-lasso_AIC_BIC.pdf" % (this_run_outputdir,args.runname))
plt.close()

# Now, Perform Lasso model selection by cross-validation
logfile.write("\n\nPerforming Lasso model selection by cross-validation\n\n")
lasso_cv_model = make_pipeline(StandardScaler(), LassoCV(cv=5, max_iter=10000)).fit(X, y) 
alpha_cv = lasso_cv_model[-1].alpha_
logfile.write("\n\nalpha chosen by CV:%s\n\n"%alpha_cv)
cv_tuned_coefs = lasso_cv_model['lassocv'].coef_

lasso = lasso_cv_model[-1]
plt.semilogx(lasso.alphas_, lasso.mse_path_, linestyle=":")
plt.plot(
    lasso.alphas_,
    lasso.mse_path_.mean(axis=-1),
    color="black",
    label="Average across the folds",
    linewidth=2,
)
plt.axvline(lasso.alpha_, linestyle="--", color="black", label="alpha: CV estimate")
plt.xlabel(r"$\alpha$")
plt.ylabel("Mean square error")
plt.legend()
plt.title("Mean square error on each fold: coordinate descent for %s" % args.runname)
plt.savefig("%s/%s-modelselection-CVtuning.pdf" % (this_run_outputdir,args.runname))
plt.close()

# plot coefficients in lasso models tuned with AIC and CV
# X.columns contains the feature names
width = 0.35       # the width of the bars
rects1 = plt.bar(X.columns, aic_tuned_coefs, -width, align='edge',
                color='tab:blue', label='AIC-tuned')
rects2 = plt.bar(X.columns, cv_tuned_coefs, +width, align='edge',
                color='tab:red', label='cv-tuned')
plt.xlabel('Feature')
plt.ylabel('Coefficient')
plt.xticks(rotation=90,fontsize=3)
plt.title('Tuned model by AIC or cross-validation for %s' % args.runname)
plt.legend()
plt.savefig("%s/%s-coefs-lasso_AIC_CV.pdf" % (this_run_outputdir,args.runname))
plt.close()

resultsdf['lasso-AIC-coefs'] = aic_tuned_coefs.tolist()
resultsdf['lasso-BIC-coefs'] = bic_tuned_coefs.tolist()
resultsdf['lasso-CV-coefs'] = cv_tuned_coefs.tolist()
resultsdf.to_csv('%s/%s-allmodels-unsorted.csv' % (this_run_outputdir,args.runname))

resultsdf.sort_values(by='lasso-AIC-coefs',inplace=True)
resultsdf.to_csv('%s/%s-allmodels-sortedbyAICtunedCoefs.csv' % (this_run_outputdir,args.runname))
resultsdf.plot(kind="bar", y=['lasso-AIC-coefs','lasso-BIC-coefs','lasso-CV-coefs'])
plt.xticks(rotation=90, fontsize=3)
plt.tick_params(axis="x",direction="in", pad=-150)
plt.legend(fontsize=6)
plt.title('Coefficients for %s' % args.runname)
plt.savefig("%s/%s-coefs-all-models-sortedbyAICtunedCoefs.pdf" % (this_run_outputdir,args.runname))
plt.close()

resultsdf.sort_values(by='human-readable',inplace=True)
resultsdf.to_csv('%s/%s-allmodels-by_terms.csv' % (this_run_outputdir,args.runname))
resultsdf.plot(kind="bar", y=['lasso-AIC-coefs','lasso-BIC-coefs','lasso-CV-coefs'])
plt.xticks(rotation=90, fontsize=3)
plt.tick_params(axis="x",direction="in", pad=-107)
plt.legend(fontsize=6)
plt.title('Coefficients for %s' % args.runname)
plt.savefig("%s/%s-coefs-all-models-by_terms.pdf" % (this_run_outputdir,args.runname))
plt.close()