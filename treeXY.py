import itertools
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

parser.add_argument("-m", "--min_depth",
                    type=tf.range_limited_int_type,
                    default=15,
                    help="Minimum depth, across all populations, for a site to be included")

parser.add_argument("-M", "--max_depth",
                    type=tf.range_limited_int_type,
                    default=200,
                    help="Maximum depth, across all populations, for a site to be included")

parser.add_argument("-a", "--min_allele_depth",
                    type=tf.range_limited_int_type,
                    default=2,
                    help="Minimum depth for an allele call")

parser.add_argument("-A", "--min_allele_pops",
                    type=tf.range_limited_int_type,
                    default=2,
                    help="Minimum number of populations required for an allele call")

parser.add_argument("--whole_scaffold",
                    action="store_true",
                    help="NOT CURRENTLY IMPLEMENTED"
                         "Take an average across all passing sites on the scaffold, rather than using windows")

parser.add_argument("-w", "--window_size",
                    type=tf.range_limited_int_type,
                    default=10000,
                    help="Size of sliding window")

parser.add_argument("-o", "--window_overlap",
                    type=tf.positive_int_type,
                    default=9000,
                    help="Overlap of consecutive sliding windows")

parser.add_argument("--binomial_adjustment",
                    action="store_true",
                    help="NOT CURRENTLY IMPLEMENTED"
                         "Adjust sample pi-within for binomial sampling")

parser.add_argument("--ignore_multiallelic",
                    action="store_true",
                    help="Ignore sites with >2 alleles, rather than removing least common allele")

parser.add_argument("--write_sync",
                    action="store_true",
                    help="Write sites passing depth checks to new SYNC file")

parser.add_argument("--write_site_stats",
                    default=False,
                    help="Also write pi-within, pi-total, and dXY for individual sites")

parser.add_argument("--dxy_trees",
                    action="store_true",
                    help="Generate UPGMA SNP trees using dXY, returning tree stats and RD (see README)")

parser.add_argument("--d_trees",
                    action="store_true",
                    help="Generate UPGMA SNP trees using Nei's D, returning tree stats and RD (see README)")

args = parser.parse_args()

##########
# treeXY #
##########

# write command line arguments to stdout
for arg in vars(args):
    print(str(arg) + ": " + str(getattr(args, arg)))

pop_names = tf.read_pop_names(args.file)
n_pops = len(pop_names)

# arg tests
tf.check_args(args.min_depth, args.max_depth, args.min_allele_pops, n_pops, args.window_size, args.window_overlap)

# generate file names from args
sync_name = args.file
sync_name = sync_name.split(".")[0]
sync_name = sync_name.split("/")[-1]

h_args = args.min_depth, args.max_depth, args.min_allele_depth, args.min_allele_pops, \
         args.window_size, args.window_overlap
h_args = list(map(str, h_args))
w_file_name = sync_name + "_" + "_".join(h_args) + "_treeXY.csv"

if args.write_sync:
    f_sync_name = sync_name + "_treeXY_filtered.sync"
    open(f_sync_name, "w").close()

if args.dxy_trees:
    site_args = args.min_depth, args.max_depth, args.min_allele_depth, args.min_allele_pops
    site_args = list(map(str, site_args))
    dxy_tree_file_name = sync_name + "_" + "_".join(site_args) + "_treeXY_topos_dXY.csv"
    open(dxy_tree_file_name, "w").close()

if args.d_trees:
    site_args = args.min_depth, args.max_depth, args.min_allele_depth, args.min_allele_pops
    site_args = list(map(str, site_args))
    da_tree_file_name = sync_name + "_" + "_".join(site_args) + "_treeXY_topos_da.csv"
    open(da_tree_file_name, "w").close()

# open window output file
with open(w_file_name, "w") as out_file:
    piw_headers = ["piw_" + name for name in pop_names]
    pit_headers = ["_".join(pair) for pair in list(itertools.combinations(pop_names, 2))]
    pit_headers = ["piT_" + name for name in pit_headers]
    dxy_headers = ["_".join(pair) for pair in list(itertools.combinations(pop_names, 2))]
    dxy_headers = ["dXY_" + name for name in dxy_headers]
    da_headers = ["_".join(pair) for pair in list(itertools.combinations(pop_names, 2))]
    da_headers = ["da_" + name for name in da_headers]
    # write header
    out_file.write("scaffold" + "," + "window_start" + "," + "window_end" + "," + "n_window_sites" + "," +
                   ",".join(piw_headers) + "," + ",".join(pit_headers) + "," + ",".join(dxy_headers) + "," +
                   ",".join(da_headers) + "\n")

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
        # filter alleles below AD threshold
        count_list = tf.filter_low_depth(count_list, args.min_allele_depth)
        pop_names = [str(i) for i in range(1, len(count_list) + 1)]

        # read depth checks
        pop_dpth = tf.check_read_depth(count_list, args.min_allele_depth, args.min_depth, args.max_depth)
        # proceed with calculations if threshold number of pops >= min_depth has been met
        if sum(pop_dpth) >= args.min_allele_pops:
            # get all possible alleles
            allele_stats = tf.check_allele_num(count_list, pop_dpth, args.min_allele_depth)
            pop_dpth = allele_stats[0]
            alleles = allele_stats[1]
            valid_allele_num = True

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

            # line with filtered counts
            f_count_list = []
            for pop_count in count_list:
                pop_count = [str(i) for i in pop_count]
                f_pop_count = ":".join(pop_count)
                f_count_list.append(f_pop_count)
            line = [scaff, pos, ref] + [str(i) for i in f_count_list]

            # proceed with piT / dXY / tree calculations for biallelic and monoallelic sites
            if valid_allele_num and len(alleles) > 0:
                # optionally, write line to filtered SYNC file
                if args.write_sync:
                    with open(f_sync_name, "a") as f_sync_file:
                        f_sync_file.write("\t".join(line) + "\n")
                # get piw, piT, and dXY
                pos_piw_vals = tf.get_site_stats(alleles, count_list, pop_names, pop_dpth)[0]
                pos_pit_vals = tf.get_site_stats(alleles, count_list, pop_names, pop_dpth)[1]
                pos_dxy_vals = tf.get_site_stats(alleles, count_list, pop_names, pop_dpth)[2]
                pos_da_vals = tf.get_site_stats(alleles, count_list, pop_names, pop_dpth)[3]
                # stats to window(s)
                window_dict = tf.stats_to_windows(window_dict, pos, max_pos, pos_piw_vals, pos_pit_vals, pos_dxy_vals,
                                                  pos_da_vals, args.window_size, args.window_overlap)

            lowest_w = next(iter(window_dict))
            keys_to_del = []

            if int(pos) >= max(lowest_w):
                w_average = tf.average_window(lowest_w, window_dict)
                if w_average:
                    tf.write_window(scaff, lowest_w, w_average, w_file_name)
                keys_to_del.append(lowest_w)

            # delete keys from window_dict if the window has already been written (to save memory)
            for key in keys_to_del:
                del window_dict[key]

            # write tree stats for biallelic sites
            n_dpth_passing_pops = len([i for i in pop_dpth if i != 0])
            if args.dxy_trees:
                if len(alleles) > 1 and n_dpth_passing_pops > 1:
                    with open(dxy_tree_file_name, "a") as tree_file:
                        tree_stats = tf.get_site_trees(pop_dpth, pos_dxy_vals)
                        tree_stats = [str(i) for i in tree_stats]
                        tree_file.write(",".join([scaff, pos, ",".join(tree_stats)]) + "\n")

            if args.d_trees:
                if len(alleles) > 1 and n_dpth_passing_pops > 1:
                    with open(da_tree_file_name, "a") as tree_file:
                        tree_stats = tf.get_site_trees(pop_dpth, pos_da_vals)
                        tree_stats = [str(i) for i in tree_stats]
                        tree_file.write(",".join([scaff, pos, ",".join(tree_stats)]) + "\n")

# write any remaining windows
for key in window_dict.keys():
    w_average = tf.average_window(key, window_dict)
    if w_average:
        tf.write_window(scaff, key, w_average, w_file_name)
