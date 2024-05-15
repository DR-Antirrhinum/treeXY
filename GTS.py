# input: treeXY table with pairwise dXY values
# process:
# 1) pairwise dXY vals to matrix
# 2) distance matrices to UPGMA trees
# 3) grouping tree scan (take random tree, compare to all trees, group)
# output: GTS forests

import datetime
import random
from argparse import ArgumentParser
import GTS_funcs as gts

# **** DELETE
import scipy.cluster.hierarchy as shc
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram
import matplotlib.pyplot as plt

##########################
# command line arguments #
##########################
parser = ArgumentParser(prog="GTS",
                        description="Grouping tree scan")

parser.add_argument("-f", "--files",
                    required=True,
                    nargs="+",
                    help="Input treeXY file(s) to be processed.",
                    metavar="treeXY_file")

parser.add_argument("-m", "--min_coverage",
                    type=float,
                    default=0.1,
                    help="Minimum read coverage for a window to be included."
                         "Read coverage is calculated as number of sites passing depth / window size.")

args = parser.parse_args()

#######
# GTS #
#######

f_date = datetime.date.today().strftime('%y_%m_%d')
f_id = str(random.getrandbits(32))
f_name_stats = "GTS_forest_stats" + "_" + str(args.min_coverage) + "_" + f_date + "_" + f_id + ".csv"
f_name_trees = "GTS_forests" + "_" + str(args.min_coverage) + "_" + f_date + "_" + f_id + ".csv"
f_name_mean_trees = "GTS_forest_mean_trees" + "_" + str(args.min_coverage) + "_" + f_date + "_" + f_id + ".tsv"

# open forest stats output file
with open(f_name_stats, "w") as out_file:
    # write header
    out_file.write("forest_ID" + "," + "forest_size" + "," + "mean_forest_SRB" + "\n")

# open forest members output file
open(f_name_trees, "w").close()

# open forest mean trees output file
# open(f_name_mean_trees, "w").close()

dxy_mat_dict = {}
dxy_clus_dict = {}
srb_dict = {}
for in_tab in args.files:
    with open(in_tab) as file:
        # derive col indices from header line
        h_line = file.readline()
        h_line = h_line.strip("\n")
        h_line = h_line.split(",")
        # piw_inds = [i for i, e in enumerate(h_line) if 'piw' in e]
        dxy_inds = [i for i, e in enumerate(h_line) if 'dXY' in e]
        dxy_h = [h_line[i] for i in dxy_inds]
        # da_inds = [i for i, e in enumerate(h_line) if 'D' in e]
        for line in file:
            # strip trailing newline
            line = line.strip("\n")
            line = line.split(",")
            scaff = line[0]
            w_start = line[1]
            w_end = line[2]
            w_size = int(w_end) - int(w_start)
            w_mid = (int(w_end) + int(w_start)) / 2
            w_mid_Mb = w_mid / 1000000
            w_identifier = scaff + "_" + str(w_mid_Mb)
            n_sites = line[3]
            coverage = int(n_sites) / w_size
            dxy = [float(line[i]) for i in dxy_inds]
            # **** need to skip regions which will give 0 height trees
            if sum(dxy) > 0 and coverage >= args.min_coverage:
                dxy_mat = gts.tab_to_df(dxy_h, dxy)
                dxy_mat_dict[w_identifier] = dxy_mat
                # UPGMA hierarchical clustering
                dxy_clusters = gts.cluster_df(dxy_mat)
                dxy_clus_dict[w_identifier] = dxy_clusters
                # get maximum height at which UPGMA merges populations i.e. the tree height
                tree_height = max(dxy_clusters[0:, 2])
                if tree_height == 0:
                    print(line)
                second_max = dxy_clusters[0:, 2][-2]
                # shortest root branch = tree height - penultimate cluster height
                srb = tree_height - second_max
                srb_dict[w_identifier] = srb

n_trees = len(dxy_clus_dict)

print(str(n_trees) + " trees above threshold coverage (" + str(args.min_coverage) + ")")

test_coph_cors = gts.test_coph_cor(100, dxy_clus_dict)

print("Mean cophenetic correlation coefficient, from comparisons of 100 random seed trees to all genomic trees, is " +
      str(sum(test_coph_cors) / len(test_coph_cors)))

gts_seeds = [i for i in dxy_clus_dict.keys()]
random.shuffle(gts_seeds)

gts_forests = gts.run_gts(gts_seeds, dxy_clus_dict, 0.8)

# tm = ((mat_dict_cp["Chr2B_383.702"] + mat_dict_cp["Chr2B_383.703"] + mat_dict_cp["Chr2B_383.704"] +
#        mat_dict_cp["Chr2B_383.705"] + mat_dict_cp["Chr2B_383.706"] + mat_dict_cp["Chr2B_383.707"] +
#        mat_dict_cp["Chr2B_383.708"] + mat_dict_cp["Chr2B_383.709"] + mat_dict_cp["Chr2B_760.273"] +
#        mat_dict_cp["Chr2B_760.274"] + mat_dict_cp["Chr2B_760.275"] + mat_dict_cp["Chr2B_786.418"] +
#        mat_dict_cp["Chr2B_786.419"] + mat_dict_cp["Chr5B_679.244"] + mat_dict_cp["Chr5B_679.245"] +
#        mat_dict_cp["Chr5B_679.246"] + mat_dict_cp["Chr5B_679.247"] + mat_dict_cp["Chr5B_679.248"] +
#        mat_dict_cp["Chr5B_679.249"] + mat_dict_cp["Chr5B_679.25" ] + mat_dict_cp["Chr5B_679.251"] +
#        mat_dict_cp["Chr5B_679.252"] + mat_dict_cp["Chr7B_662.065"] + mat_dict_cp["Chr7B_662.066"] +
#        mat_dict_cp["Chr7B_662.067"] + mat_dict_cp["Chr7B_662.068"]) / 26)

# clusters = shc.linkage(tm, method="average", metric="euclidean")
#
# shc.dendrogram(Z=clusters)
# plt.show()
# exit()

write_list = []
f_trees_list = []
mean_tree_list = []
for i, forest in enumerate(gts_forests):
    forest_id = str(i + 1)
    forest_size = str(len(forest))
    forest_trees = []
    forest_srb = []
    for tree in forest:
        tree_mat = dxy_mat_dict[tree]
        tree_srb = srb_dict[tree]
        forest_trees.append(tree_mat)
        forest_srb.append(tree_srb)

    mean_forest_tree = sum(forest_trees) / len(forest_trees)
    mean_forest_srb = sum(forest_srb) / len(forest_srb)
    write_list.append(forest_id + "," + forest_size + "," + str(mean_forest_srb) + "\n")
    f_trees_list.append(",".join(forest) + "\n")
    # convert mean_forest_tree to newick format
    mf_cluster = gts.cluster_df(mean_forest_tree)
    mf_node = shc.to_tree(mf_cluster, False)
    mf_newick = gts.get_newick(mf_node, mf_node.dist, gts.get_pop_names(dxy_h))
    mean_tree_list.append(mf_newick)

with(open(f_name_stats, "a")) as out_file:
    for line in write_list:
        out_file.write(line)

with(open(f_name_trees, "a")) as out_file:
    for line in f_trees_list:
        out_file.write(line)

with(open(f_name_mean_trees, "a")) as out_file:
    # this needs to be written as TSV, because Newick format contains commas
    out_file.write("\t".join(mean_tree_list))

# **** try writing info through logging e.g. n_trees
# **** write function to test mean coph cor
