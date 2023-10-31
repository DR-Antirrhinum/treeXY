import itertools
# import tracemalloc
from argparse import ArgumentParser
import treeXY_funcs as tf


##########################
# command line arguments #
##########################
parser = ArgumentParser(prog="treeXY",
                        description="Calculate window-averaged pi-statistics from mapped pool-seq data")

parser.add_argument("-f", "--file",
                    required=True,
                    help="Input SYNC file to be processed.",
                    metavar="sync_file")

# parser.add_argument("-p", "--pops",
#                     required=False,
#                     help="Input pop file describing the pools. "
#                          "Column 1 is pop names, column 2 is number in pool, and "
#                          "column 3 is 1/2 grouping (or 0 to exclude).",
#                     metavar="pop_file")

parser.add_argument("-m", "--min_depth",
                    type=int,
                    default=15,
                    help="Minimum depth, across all populations, for a site to be included")

parser.add_argument("-M", "--max_depth",
                    type=int,
                    default=200,
                    help="Maximum depth, across all populations, for a site to be included")

parser.add_argument("-a", "--min_allele_depth",
                    type=int,
                    default=2,
                    help="Minimum depth for an allele call")

parser.add_argument("-A", "--min_allele_pops",
                    type=int,
                    default=2,
                    help="Minimum number of populations required for an allele call")

parser.add_argument("-w", "--window_size",
                    type=int,
                    default=10000,
                    help="Size of sliding window")

parser.add_argument("-o", "--window_overlap",
                    type=int,
                    default=9000,
                    help="Overlap of consecutive sliding windows")

parser.add_argument("--ignore_multiallelic",
                    action="store_true",
                    help="Ignore sites with >2 alleles, rather than removing least common allele")

parser.add_argument("--write_sync",
                    action="store_true",
                    help="Write sites passing depth checks to new SYNC file")

parser.add_argument("--write_site_stats",
                    default=False,
                    help="Also write pi-within, pi-total, and dXY for individual sites")

parser.add_argument("--compute_trees",
                    action="store_true",
                    help="Run hierarchical clustering on site dXY values, returning tree stats and RD (see README)")

args = parser.parse_args()

##########
# treeXY #
##########

# tracemalloc.start()

# generate file names from args
sync_name = args.file
sync_name = sync_name.split(".")[0]
sync_name = sync_name.split("/")[-1]

h_args = args.min_depth, args.max_depth, args.min_allele_depth, args.min_allele_pops, \
         args.window_size, args.window_overlap
h_args = list(map(str, h_args))
file_name = sync_name + "_" + "_".join(h_args) + "_treeXY.csv"

if args.write_sync:
    f_sync_name = sync_name + "_treeXY_filtered.sync"
    open(f_sync_name, "w").close()

if args.compute_trees:
    site_args = args.min_depth, args.max_depth, args.min_allele_depth, args.min_allele_pops
    site_args = list(map(str, site_args))
    tree_file_name = sync_name + "_" + "_".join(site_args) + "_treeXY_topos.csv"
    open(tree_file_name, "w").close()

# initialise dict to store window averages
window_details = tf.initialise_windows(args.file, args.window_size, args.window_overlap)
window_dict = window_details[0]
max_pos = window_details[1]

with open(args.file) as file:
    for line in file:
        # strip trailing newline
        line = line.strip("\n")
        line = line.split("\t")
        # name of scaffold / chromosome
        scaff = line[0]
        # genomic position
        pos = line[1]
        # base call in reference
        ref = line[2]
        # all colon delimited count columns
        site_counts = line[3:]
        count_list = tf.get_sync_counts(site_counts)
        # **** eventually, want to read names from file
        pop_names = [str(i) for i in range(1, len(count_list) + 1)]

        # read depth checks
        pop_dpth = tf.check_read_depth(count_list, args.min_allele_depth, args.min_depth, args.max_depth)
        # in SW, all pops have to have > threshold depth for a site to be included
        # **** address this more elegantly
        if sum(pop_dpth) == len(pop_names):
            # get all possible alleles
            allele_stats = tf.check_allele_num(count_list, pop_dpth, args.min_allele_depth)
            pop_dpth = allele_stats[0]
            alleles = allele_stats[1]
            valid_allele_num = True

            # print("starting position " + pos)
            # print(tracemalloc.get_traced_memory())

            if len(alleles) > 2 and args.ignore_multiallelic:
                valid_allele_num = False
            else:
                while len(alleles) > 2:
                    line = tf.filter_triallelic(count_list, scaff, pos, ref)
                    # strip trailing newline
                    line = line.strip("\n")
                    line = line.split("\t")
                    # all colon delimited count columns
                    site_counts = line[3:]
                    count_list = tf.get_sync_counts(site_counts)
                    # get alleles
                    allele_stats = tf.check_allele_num(count_list, pop_dpth, args.min_allele_depth)
                    pop_dpth = allele_stats[0]
                    alleles = allele_stats[1]

            # print("completed multiallelic filtering for position " + pos)
            # print(tracemalloc.get_traced_memory())

            # proceed with piT / dXY / tree calculations for biallelic and monoallelic sites
            if valid_allele_num and len(alleles) > 0 and all([i > 0 for i in pop_dpth]):
                # optionally, write line to filtered SYNC file
                if args.write_sync:
                    with open(f_sync_name, "a") as f_sync_file:
                        f_sync_file.write("\t".join(line) + "\n")
                # get piw, piT, and dXY
                pos_piw_vals = tf.get_site_stats(alleles, count_list, pop_names, pop_dpth)[0]
                pos_pit_vals = tf.get_site_stats(alleles, count_list, pop_names, pop_dpth)[1]
                pos_dxy_vals = tf.get_site_stats(alleles, count_list, pop_names, pop_dpth)[2]
                pos_D_vals = tf.get_site_stats(alleles, count_list, pop_names, pop_dpth)[3]
                # stats to window(s)
                window_dict = tf.stats_to_windows(window_dict, pos, max_pos, pos_piw_vals, pos_pit_vals, pos_dxy_vals,
                                                  pos_D_vals, args.window_size, args.window_overlap)

            if args.compute_trees:
                if len(alleles) > 1 and all([i > 0 for i in pop_dpth]):
                    with open(tree_file_name, "a") as tree_file:
                        tree_stats = tf.get_site_trees(pop_dpth, pos_dxy_vals)
                        tree_stats = [str(i) for i in tree_stats]
                        tree_file.write(",".join([scaff, pos, ",".join(tree_stats)]) + "\n")

# print(tracemalloc.get_traced_memory())

# open file for writing
with open(file_name, "w") as out_file:
    piw_headers = ["piw_" + name for name in pop_names]
    pit_headers = ["_".join(pair) for pair in list(itertools.combinations(pop_names, 2))]
    pit_headers = ["piT_" + name for name in pit_headers]
    dxy_headers = ["_".join(pair) for pair in list(itertools.combinations(pop_names, 2))]
    dxy_headers = ["dXY_" + name for name in dxy_headers]
    D_headers = ["_".join(pair) for pair in list(itertools.combinations(pop_names, 2))]
    D_headers = ["D_" + name for name in D_headers]
    # write header
    out_file.write("scaffold" + "," + "window_start" + "," + "window_end" + "," + "n_window_sites" + "," +
                   ",".join(piw_headers) + "," + ",".join(pit_headers) + "," + ",".join(dxy_headers) + "," +
                   ",".join(D_headers) + "\n")

    for key in window_dict.keys():
        # need to take mean of each column for each window, for piw and dxy
        window_piw_vals = window_dict[key][0]
        window_pit_vals = window_dict[key][1]
        window_dxy_vals = window_dict[key][2]
        window_D_vals = window_dict[key][3]

        if len(window_piw_vals) > 0 and len(window_pit_vals) > 0 and len(window_dxy_vals) > 0:
            window_piw_means = tf.vals_to_pop_means(window_piw_vals)
            window_pit_means = tf.vals_to_pop_means(window_pit_vals)
            window_dxy_means = tf.vals_to_pop_means(window_dxy_vals)
            window_D_means = tf.vals_to_pop_means(window_D_vals)

            # convert to string and write to file
            window_piw_means = list(map(str, window_piw_means))
            window_pit_means = list(map(str, window_pit_means))
            window_dxy_means = list(map(str, window_dxy_means))
            window_D_means = list(map(str, window_D_means))
            # write to file
            out_file.write(scaff + "," + str(min(key)) + "," + str(max(key)) + "," +
                           str(len(window_piw_vals)) + "," + ",".join(window_piw_means) + "," +
                           ",".join(window_pit_means) + "," + ",".join(window_dxy_means) + "," +
                           ",".join(window_D_means) + "\n")
