import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as shc
# from scipy.spatial import distance_matrix
from scipy.spatial.distance import squareform


#############
# functions #
#############
def tab_to_matrix(dxy_header, pairwise_dxy):
    pop_pairs = []
    pop_names = []
    for i, col in enumerate(dxy_header):
        col = col.split("_")
        pop1 = int(col[1])
        pop2 = int(col[2])
        pop_names.append(pop1)
        pop_names.append(pop2)
        dxy = pairwise_dxy[i]
        pop_pairs.append([pop1, pop2, dxy])

    pop_names = list(set(pop_names))

    # read pairwise dXY into pandas data frame, and convert to distance matrix
    df = pd.DataFrame(pop_pairs, columns=["pop1", "pop2", "dXY"])
    df = df.pivot_table(index='pop1', columns='pop2', values='dXY')
    df = df.combine_first(df.T)
    # set diagonal to 0
    np.fill_diagonal(df.to_numpy(), 0)
    # convert distance matrix to minimal representation for scipy linkage compatibility
    dm = squareform(df)

    # UPGMA hierarchical clustering
    # can be visualised as a tree using shc.dendrogram
    clusters = shc.linkage(dm, method="average", metric="euclidean")

    return clusters


def coph_cor(cluster_1, cluster_2):
    coph1 = shc.cophenet(cluster_1)
    coph2 = shc.cophenet(cluster_2)
    coph_mat = np.corrcoef(coph1, coph2)

    return coph_mat[0][1]


def grouping_tree_scan(seed, tree_dict, coph_threshold):
    seed_tree = tree_dict[seed]
    forest = []
    for tree_key in tree_dict.keys():
        cc = coph_cor(seed_tree, tree_dict[tree_key])
        if cc >= coph_threshold:
            forest.append(tree_key)

    return forest
