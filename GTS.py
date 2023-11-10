# input: treeXY table with pairwise dXY values
# process:
# 1) pairwise dXY vals to matrix
# 2) distance matrices to UPGMA trees
# 3) grouping tree scan (take random tree, compare to all trees, group)
# output: GTS forests

from argparse import ArgumentParser
import GTS_funcs as gts

# **** TO REMOVE
import itertools
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as shc
# from scipy.spatial import distance_matrix
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram
import matplotlib.pyplot as plt
import networkx as nx

##########################
# command line arguments #
##########################
parser = ArgumentParser(prog="GTS",
                        description="Grouping tree scan")

parser.add_argument("-f", "--file",
                    required=True,
                    help="Input treeXY file to be processed.",
                    metavar="treeXY_file")

parser.add_argument("-m", "--min_sites",
                    type=int,
                    default=15,
                    help="Minimum number of sites for a window to be included")

args = parser.parse_args()

#######
# GTS #
#######

with open(args.file) as file:
    # derive col indices from header line
    h_line = file.readline()
    h_line = h_line.strip("\n")
    h_line = h_line.split(",")
    # piw_inds = [i for i, e in enumerate(h_line) if 'piw' in e]
    dxy_inds = [i for i, e in enumerate(h_line) if 'dXY' in e]
    dxy_h = [h_line[i] for i in dxy_inds]
    # da_inds = [i for i, e in enumerate(h_line) if 'D' in e]
    dxy_mat_list = []
    srb_list = []
    for line in file:
        # strip trailing newline
        line = line.strip("\n")
        line = line.split(",")
        scaff = line[0]
        w_start = line[1]
        w_end = line[2]
        # **** arg to filter on n_sites
        n_sites = line[3]
        dxy = [float(line[i]) for i in dxy_inds]
        dxy_mat = gts.tab_to_matrix(dxy_h, dxy)
        dxy_mat_list.append(dxy_mat)
        # get maximum height at which UPGMA merges populations i.e. the tree height
        # **** implement piw and Nei's D
        tree_height = max(dxy_mat[0:, 2])
        second_max = dxy_mat[0:, 2][-2]
        # shortest root branch = tree height - penultimate cluster height
        srb = tree_height - second_max
        srb_list.append(srb)

# 67270 (30, 1946) nan
# 67271 (30, 1947) nan
# 67272 (30, 1948) nan
# 67273 (30, 1949) nan
# 67274 (30, 1950) nan

# shc.dendrogram(Z=dxy_mat_list[30])
# plt.show()
# exit()
#
# quit()

print(srb_list[0])

all_comps = itertools.combinations([i for i in range(0, len(dxy_mat_list))], 2)
tot_comps = len([i for i in all_comps])
all_comps = itertools.combinations([i for i in range(0, len(dxy_mat_list))], 2)

coph_list = []
for i, comp in enumerate(all_comps):
    tree_1 = comp[0]
    tree_2 = comp[1]
    coph_cor = gts.coph_cor(dxy_mat_list[tree_1], dxy_mat_list[tree_2])
    coph_list.append([tree_1, tree_2, coph_cor])

df = pd.DataFrame(coph_list, columns=["tree_1", "tree_2", "coph_cor"])
# df = df.pivot_table(index='tree_1', columns='tree_2', values='coph_cor')
# df = df.combine_first(df.T)
# set diagonal to 0
# np.fill_diagonal(df.to_numpy(), 0)
# convert distance matrix to minimal representation for scipy linkage compatibility
# dm = squareform(df)

G = nx.from_pandas_edgelist(df, source="tree_1", target="tree_2", edge_attr="coph_cor")
nx.draw(G)
plt.show()

# coph_list = [str(i) for i in coph_list]

# print(",".join(coph_list))
