import numpy as np
import pandas as pd
import random
import scipy.cluster.hierarchy as shc
# from scipy.spatial import distance_matrix
from scipy.spatial.distance import squareform

# **** DELETE
import scipy.cluster.hierarchy as shc
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram
import matplotlib.pyplot as plt


#############
# functions #
#############
def get_pop_names(dxy_header):
    pop_names = []
    for i, col in enumerate(dxy_header):
        col = col.split("_")
        pop1 = int(col[1])
        pop2 = int(col[2])
        pop_names.append(pop1)
        pop_names.append(pop2)

    pop_names = list(set(pop_names))

    return pop_names


def tab_to_df(dxy_header, pairwise_dxy):
    pop_pairs = []
    for i, col in enumerate(dxy_header):
        col = col.split("_")
        pop1 = int(col[1])
        pop2 = int(col[2])
        dxy = pairwise_dxy[i]
        pop_pairs.append([pop1, pop2, dxy])

    # read pairwise dXY into pandas data frame, and convert to distance matrix
    df = pd.DataFrame(pop_pairs, columns=["pop1", "pop2", "dXY"])
    df = df.pivot_table(index='pop1', columns='pop2', values='dXY')
    df = df.combine_first(df.T)
    # set diagonal to 0
    np.fill_diagonal(df.to_numpy(), 0)

    return df


def cluster_df(df):
    # convert distance matrix to minimal representation for scipy linkage compatibility
    dm = squareform(df)
    cluster = shc.linkage(dm, method="average", metric="euclidean")

    return cluster


def coph_cor(cluster_1, cluster_2):
    coph1 = shc.cophenet(cluster_1)
    coph2 = shc.cophenet(cluster_2)
    coph_mat = np.corrcoef(coph1, coph2)

    return coph_mat[0][1]


def test_coph_cor(n_comps, tree_dict):
    # run a number of test cophenetic correlation calculations between random trees
    # can be used to estimate a suitable cutoff for the GTS
    seed_list = random.sample(tree_dict.keys(), n_comps)
    coph_list = []
    for seed in seed_list:
        seed_tree = tree_dict[seed]
        for comp_tree in tree_dict.values():
            coph_list.append(coph_cor(seed_tree, comp_tree))

    return coph_list


def grouping_tree_scan(seed, tree_dict, coph_threshold):
    seed_tree = tree_dict[seed]
    forest = []
    for tree_key in tree_dict.keys():
        cc = coph_cor(seed_tree, tree_dict[tree_key])
        if cc >= coph_threshold:
            forest.append(tree_key)

    return forest


def run_gts(seed_list, tree_dict, coph_threshold):
    r_tree_dict = tree_dict
    forest_list = []
    while len(seed_list) > 0:
        curr_seed = seed_list[0]
        forest = grouping_tree_scan(curr_seed, r_tree_dict, coph_threshold)
        forest_list.append(forest)
        for f_tree in forest:
            seed_list.remove(f_tree)
            del r_tree_dict[f_tree]

    return forest_list


def get_newick(node, parent_dist, leaf_names, newick='') -> str:
    """
    From SO user jfn (https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format)
    Convert scipy.cluster.hierarchy.to_tree()-output to Newick format.

    :param node: output of sciply.cluster.hierarchy.to_tree()
    :param parent_dist: output of sciply.cluster.hierarchy.to_tree().dist
    :param leaf_names: list of leaf names
    :param newick: leave empty, this variable is used in recursion.
    :returns: tree in Newick format
    """
    if node.is_leaf():
        return "%s:%.6f%s" % (leaf_names[node.id], parent_dist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.6f%s" % (parent_dist - node.dist, newick)
        else:
            newick = ");"
        newick = get_newick(node.get_left(), node.dist, leaf_names, newick=newick)
        newick = get_newick(node.get_right(), node.dist, leaf_names, newick=",%s" % (newick))
        newick = "(%s" % (newick)
        return newick
