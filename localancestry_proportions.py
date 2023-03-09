#! /usr/bin/env python

# Adapted from:
# https://github.com/agoldberglab/CV_DuffySelection
# https://doi.org/10.7554/eLife.63177

import msprime, pyslim
import numpy as np
import pandas as pd
import re
import sys
import tskit

infile = sys.argv[1]
outfile = re.search("(.*).trees", infile).group(1)

ts = tskit.load(infile)

popA_nodes = ts.samples(population=1)
popB_nodes = ts.samples(population=2)
popC_nodes = ts.samples(population=3)
popD_nodes = ts.samples(population=4)
admix_inds = pyslim.individuals_alive_at(ts, 0, population=9)
max_time = ts.max_root_time - 1

def get_admixture_proportions(ts, admix_inds, popA_nodes, popB_nodes, popC_nodes, popD_nodes, max_time):
	proportions_A = []
	proportions_B = []
	proportions_C = []
	proportions_D = []
	proportions_Total = []
	admix_nodes = np.array([ts.individual(ind).nodes for ind in admix_inds]).flatten()
	edges = ts.tables.link_ancestors(admix_nodes, np.concatenate((popA_nodes, popB_nodes, popC_nodes, popD_nodes)))
	for ind in admix_inds:
		child_nodes = ts.individual(ind).nodes
		A_edges = edges[np.isin(edges.child, child_nodes) & np.isin(edges.parent, popA_nodes)]
		B_edges = edges[np.isin(edges.child, child_nodes) & np.isin(edges.parent, popB_nodes)]
		C_edges = edges[np.isin(edges.child, child_nodes) & np.isin(edges.parent, popC_nodes)]
		D_edges = edges[np.isin(edges.child, child_nodes) & np.isin(edges.parent, popD_nodes)]
		span_A = np.sum(A_edges.right - A_edges.left)
		span_B = np.sum(B_edges.right - B_edges.left)
		span_C = np.sum(C_edges.right - C_edges.left)
		span_D = np.sum(D_edges.right - D_edges.left)
		proportions_A.append(span_A/(span_A + span_B + span_C + span_D))
		proportions_B.append(span_B / (span_A + span_B + span_C + span_D))
		proportions_C.append(span_C / (span_A + span_B + span_C + span_D))
		proportions_D.append(span_D / (span_A + span_B + span_C + span_D))
		proportions_Total.append((span_A + span_B + span_C + span_D)/(span_A + span_B + span_C + span_D))
	ancestry_arr = np.stack((proportions_A, proportions_B, proportions_C, proportions_D), axis=-1)
	ancestry_df = pd.DataFrame(ancestry_arr, columns=["P1AncestryProportion", "P2AncestryProportion", "P3AncestryProportion", "P4AncestryProportion"])
	outname = f"{outfile}.csv"
	# save to csv file
	ancestry_df.to_csv(outname)

get_admixture_proportions(ts, admix_inds, popA_nodes, popB_nodes, popC_nodes, popD_nodes, max_time)

# python "localancestry_proportions - Copy.py" "E:/git/SLiM_stickles/output/sticklebacks.trees"
