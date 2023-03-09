#! /usr/bin/env python

# Adapted from:
# https://github.com/agoldberglab/CV_DuffySelection
# https://doi.org/10.7554/eLife.63177

import msprime, pyslim
import numpy as np
import pandas as pd
import re
import sys

infile = sys.argv[1]
outfile = re.search("(.*).trees", infile).group(1)

ts = pyslim.load(infile)

popA_nodes = ts.samples(population=1)
popB_nodes = ts.samples(population=2)
popC_nodes = ts.samples(population=3)
popD_nodes = ts.samples(population=4)
popE_nodes = ts.samples(population=5)
popF_nodes = ts.samples(population=6)
popG_nodes = ts.samples(population=7)
popH_nodes = ts.samples(population=8)
admix_inds = ts.individuals_alive_at(0, population=9)
max_time = ts.max_root_time - 1

def get_admixture_proportions(ts, admix_inds, popA_nodes, popB_nodes, popC_nodes, popD_nodes, popE_nodes, popF_nodes, popG_nodes, popH_nodes, max_time):
	proportions_A = []
	proportions_B = []
	proportions_C = []
	proportions_D = []
	proportions_E = []
	proportions_F = []
	proportions_G = []
	proportions_H = []
	proportions_Total = []
	admix_nodes = np.array([ts.individual(ind).nodes for ind in admix_inds]).flatten()
	edges = ts.tables.link_ancestors(admix_nodes, np.concatenate((popA_nodes, popB_nodes, popC_nodes, popD_nodes, popE_nodes, popF_nodes, popG_nodes, popH_nodes)))
	for ind in admix_inds:
		child_nodes = ts.individual(ind).nodes
		A_edges = edges[np.isin(edges.child, child_nodes) & np.isin(edges.parent, popA_nodes)]
		B_edges = edges[np.isin(edges.child, child_nodes) & np.isin(edges.parent, popB_nodes)]
		C_edges = edges[np.isin(edges.child, child_nodes) & np.isin(edges.parent, popC_nodes)]
		D_edges = edges[np.isin(edges.child, child_nodes) & np.isin(edges.parent, popD_nodes)]
		E_edges = edges[np.isin(edges.child, child_nodes) & np.isin(edges.parent, popE_nodes)]
		F_edges = edges[np.isin(edges.child, child_nodes) & np.isin(edges.parent, popF_nodes)]
		G_edges = edges[np.isin(edges.child, child_nodes) & np.isin(edges.parent, popG_nodes)]
		H_edges = edges[np.isin(edges.child, child_nodes) & np.isin(edges.parent, popH_nodes)]
		span_A = np.sum(A_edges.right - A_edges.left)
		span_B = np.sum(B_edges.right - B_edges.left)
		span_C = np.sum(C_edges.right - C_edges.left)
		span_D = np.sum(D_edges.right - D_edges.left)
		span_E = np.sum(E_edges.right - E_edges.left)
		span_F = np.sum(F_edges.right - F_edges.left)
		span_G = np.sum(G_edges.right - G_edges.left)
		span_H = np.sum(H_edges.right - H_edges.left)
		proportions_A.append(span_A / (span_A + span_B + span_C + span_D + span_E + span_F + span_G + span_H))
		proportions_B.append(span_B / (span_A + span_B + span_C + span_D + span_E + span_F + span_G + span_H))
		proportions_C.append(span_C / (span_A + span_B + span_C + span_D + span_E + span_F + span_G + span_H))
		proportions_D.append(span_D / (span_A + span_B + span_C + span_D + span_E + span_F + span_G + span_H))
		proportions_E.append(span_E / (span_A + span_B + span_C + span_D + span_E + span_F + span_G + span_H))
		proportions_F.append(span_F / (span_A + span_B + span_C + span_D + span_E + span_F + span_G + span_H))
		proportions_G.append(span_G / (span_A + span_B + span_C + span_D + span_E + span_F + span_G + span_H))
		proportions_H.append(span_H / (span_A + span_B + span_C + span_D + span_E + span_F + span_G + span_H))
		proportions_Total.append((span_A + span_B + span_C + span_D + span_E + span_F + span_G + span_H)/(span_A + span_B + span_C + span_D + span_E + span_F + span_G + span_H))
	ancestry_arr = np.stack((proportions_A, proportions_B, proportions_C, proportions_D, proportions_E, proportions_F, proportions_G, proportions_H), axis=-1)
	ancestry_df = pd.DataFrame(ancestry_arr, columns=["P1AncestryProportion", "P2AncestryProportion", "P3AncestryProportion", "P4AncestryProportion", "P5AncestryProportion", "P6AncestryProportion", "P7AncestryProportion", "P8AncestryProportion"])
	outname = f"{outfile}_ancestryproportions.csv"
	# save to csv file
	ancestry_df.to_csv(outname)

get_admixture_proportions(ts, admix_inds, popA_nodes, popB_nodes, popC_nodes, popD_nodes, popE_nodes, popF_nodes, popG_nodes, popH_nodes, max_time)

# python "localancestry_proportions - F.py" "E:/git/SLiM_stickles/output/sticklebacks.trees"
