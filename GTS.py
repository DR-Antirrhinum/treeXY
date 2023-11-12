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
f_name_stats = "GTS_forest_stats" + "_" + f_date + "_" + f_id + ".csv"
f_name_trees = "GTS_forests" + "_" + f_date + "_" + f_id + ".csv"

# open forest stats output file
with open(f_name_stats, "w") as out_file:
    # write header
    out_file.write("forest_ID" + "," + "forest_size" + "," + "mean_forest_SRB" + "\n")

# open forest members output file
open(f_name_trees, "w").close()

dxy_mat_dict = {}
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
            if coverage >= args.min_coverage:
                dxy_mat = gts.tab_to_matrix(dxy_h, dxy)
                dxy_mat_dict[w_identifier] = dxy_mat
                # get maximum height at which UPGMA merges populations i.e. the tree height
                tree_height = max(dxy_mat[0:, 2])
                second_max = dxy_mat[0:, 2][-2]
                # shortest root branch = tree height - penultimate cluster height
                srb = tree_height - second_max
                srb_dict[w_identifier] = (srb)

gts_seeds = [i for i in dxy_mat_dict.keys()]
random.shuffle(gts_seeds)

forest_list = []
while len(gts_seeds) > 0:
    curr_seed = gts_seeds[0]
    forest = gts.grouping_tree_scan(curr_seed, dxy_mat_dict, 0.5)
    forest_list.append(forest)
    for f_tree in forest:
        gts_seeds.remove(f_tree)
        del dxy_mat_dict[f_tree]

write_list = []
f_trees_list = []
for i, forest in enumerate(forest_list):
    forest_id = str(i + 1)
    forest_size = str(len(forest))
    forest_srb = []
    for tree in forest:
        tree_srb = srb_dict[tree]
        forest_srb.append(tree_srb)

    mean_forest_srb = sum(forest_srb) / len(forest_srb)
    write_list.append(forest_id + "," + forest_size + "," + str(mean_forest_srb) + "\n")
    f_trees_list.append(",".join(forest) + "\n")

with(open(f_name_stats, "a")) as out_file:
    for line in write_list:
        out_file.write(line)

with(open(f_name_trees, "a")) as out_file:
    for line in f_trees_list:
        out_file.write(line)
