#! /usr/bin/env python

#Adapted from:
#https://github.com/agoldberglab/CV_DuffySelection
#https://doi.org/10.7554/eLife.63177

import msprime, pyslim
import numpy as np
import pandas as pd
import re
import sys


infile = sys.argv[1]
outfile = re.search("(.*).trees", infile).group(1)

ts = pyslim.load(infile).simplify()

#create zero arrays for genomic positions
#and corresponding ancestry proportions
breaks = np.zeros(ts.num_trees + 1)
ancestry1 = np.zeros(ts.num_trees + 1)
ancestry2 = np.zeros(ts.num_trees + 1)
ancestry3 = np.zeros(ts.num_trees + 1)
ancestry4 = np.zeros(ts.num_trees + 1)
ancestry5 = np.zeros(ts.num_trees + 1)
ancestry6 = np.zeros(ts.num_trees + 1)
ancestry7 = np.zeros(ts.num_trees + 1)
ancestry8 = np.zeros(ts.num_trees + 1)

#iterate over trees (1 unique geneology/genomic interval)
for tree in ts.trees():
	subpop_sum, subpop_sum2, subpop_sum3, subpop_sum4, subpop_sum5, subpop_sum6, subpop_sum7, subpop_sum8, subpop_weights = 0, 0, 0, 0, 0, 0, 0, 0, 0

	#iterate over ancestors [roots]
	for root in tree.roots:
		#for each ancestor, count descendents
		leaves_count = tree.num_samples(root) - 1

		#get the pop of ancestor (p1=0, p2=1)
		#subpop_sum = num chromosomes for each genomic position w/ p1 ancestry
		if tree.population(root) == 0:
			subpop_sum += (tree.population(root) + 1) * leaves_count

		elif tree.population(root) == 1:
			subpop_sum2 += tree.population(root) * leaves_count

		elif tree.population(root) == 2:
			subpop_sum3 += tree.population(root) * leaves_count / 2

		elif tree.population(root) == 3:
			subpop_sum4 += tree.population(root) * leaves_count / 3

		elif tree.population(root) == 4:
			subpop_sum5 += tree.population(root) * leaves_count / 4

		elif tree.population(root) == 5:
			subpop_sum6 += tree.population(root) * leaves_count / 5

		elif tree.population(root) == 6:
			subpop_sum7 += tree.population(root) * leaves_count / 6

		elif tree.population(root) == 7:
			subpop_sum8 += tree.population(root) * leaves_count / 7

		#subpop_weights = num of chromosomes in sample
		subpop_weights += leaves_count
	#genomic position corresponds to left end of tree interval
	breaks[tree.index] = tree.interval[0]

	#ancestry proportion = num chroms with p1 ancestry / num chroms
	ancestry1[tree.index] = subpop_sum / subpop_weights
	ancestry2[tree.index] = subpop_sum2 / subpop_weights
	ancestry3[tree.index] = subpop_sum3 / subpop_weights
	ancestry4[tree.index] = subpop_sum4 / subpop_weights
	ancestry5[tree.index] = subpop_sum5 / subpop_weights
	ancestry6[tree.index] = subpop_sum6 / subpop_weights
	ancestry7[tree.index] = subpop_sum7 / subpop_weights
	ancestry8[tree.index] = subpop_sum8 / subpop_weights

#fill chromosome end positions
breaks[-1] = ts.sequence_length
ancestry1[-1] = ancestry2[-2]
ancestry2[-1] = ancestry2[-2]
ancestry3[-1] = ancestry3[-2]
ancestry4[-1] = ancestry2[-2]
ancestry5[-1] = ancestry2[-2]
ancestry6[-1] = ancestry3[-2]
ancestry7[-1] = ancestry2[-2]
ancestry8[-1] = ancestry3[-2]

#transform arrays to pandas dataframe
ancestry_arr = np.stack((breaks, ancestry1, ancestry2, ancestry3, ancestry4, ancestry5, ancestry6, ancestry7, ancestry8), axis = -1)
ancestry_df = pd.DataFrame(ancestry_arr, columns = ["GenomicPosition", "P1AncestryProportion", "P2AncestryProportion", "P3AncestryProportion", "P4AncestryProportion", "P5AncestryProportion", "P6AncestryProportion", "P7AncestryProportion", "P8AncestryProportion"])

outname = f"{outfile}_ancestryproportions.csv"

#save to csv file
ancestry_df.to_csv(outname)
