import random
import itertools
import tracemalloc
from argparse import ArgumentParser
import treeXY_funcs as tf
import math as maths

##########################
# command line arguments #
##########################
parser = ArgumentParser(prog="treeXY",
                        description="Calculate pi-statistics, and compute UPGMA trees, from mapped pool-seq data")

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
                    default=15,
                    help="Minimum depth, across all populations, for a site to be included")

parser.add_argument("-M", "--max_depth",
                    default=500,
                    help="Maximum depth, across all populations, for a site to be included")

parser.add_argument("-a", "--min_allele_depth",
                    default=2,
                    help="Minimum depth for an allele call")

parser.add_argument("-A", "--min_allele_pops",
                    default=2,
                    help="Minimum number of populations required for an allele call")

parser.add_argument("-w", "--window_size",
                    default=10000,
                    help="Size of sliding window")

parser.add_argument("-o", "--window_overlap",
                    default=9000,
                    help="Overlap of consecutive sliding windows")

args = parser.parse_args()

##########
# treeXY #
##########

# tracemalloc.start()

window_details = tf.initialise_windows()
window_dict = window_details[0]
max_pos = window_details[1]

# dict to store all stats for each position
# pos_stats_dict = {}

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
        # dpth_dict = {}
        pop_dpth = tf.check_read_depth(count_list)
        # retain for later window averaging
        # dpth_dict[pos] = pop_dpth

        # get all possible alleles
        # **** maybe instead of summing, I should do the per-pop depth measuring from the start?
        # **** otherwise, I have to redefine alleles later
        allele_stats = tf.check_allele_num(count_list, pop_dpth)
        pop_dpth = allele_stats[0]
        alleles = allele_stats[1]

        # print("starting position " + pos)
        # print(tracemalloc.get_traced_memory())

        while len(alleles) > 2:
            line = tf.filter_triallelic(count_list, scaff, pos, ref)
            # strip trailing newline
            line = line.strip("\n")
            line = line.split("\t")
            # all colon delimited count columns
            site_counts = line[3:]
            count_list = tf.get_sync_counts(site_counts)
            # get alleles
            allele_stats = tf.check_allele_num(count_list, pop_dpth)
            pop_dpth = allele_stats[0]
            alleles = allele_stats[1]

        # print("completed multiallelic filtering for position " + pos)
        # print(tracemalloc.get_traced_memory())

        # proceed with piT / dXY / tree calculations for biallelic and monoallelic sites
        # **** it is only here that I start to ignore lines below read depth threshold, does that make sense?
        # **** I could provide option to output filtered SYNC file here
        # **** if monoallelic, piT and dXY will always evaluate to zero, but need to include in output anyway
        # print(tf.check_pop_allele_depth(pop_dpth, site_counts, alleles))

        if len(alleles) > 0:
            # **** another check required here?

            # initialise dict key
            # print(pos, "PREDICT", tracemalloc.get_traced_memory())
            # **** creating this dict entry has massive memory cost in some instances
            # **** in test.sync, memory always spikes at position 975314 - why that position specifically?
            # pos_stats_dict[pos] = []
            # print(pos, "POSTDICT", tracemalloc.get_traced_memory())
            # get piw, piT, and dXY
            pop_piw_dict = tf.get_site_stats(alleles, count_list, pop_names, pop_dpth)[0]
            pop_pit_dict = tf.get_site_stats(alleles, count_list, pop_names, pop_dpth)[1]
            pop_dxy_dict = tf.get_site_stats(alleles, count_list, pop_names, pop_dpth)[2]
            # compile stats
            # pos_stats_dict[pos].append(pop_piw_dict)
            # pos_stats_dict[pos].append(pop_pit_dict)
            # pos_stats_dict[pos].append(pop_dxy_dict)
            # stats to window(s)
            # **** need to generate vals lists before populating window dict
            pos_piw_vals = tf.dict_to_vals(pop_piw_dict)
            pos_pit_vals = tf.dict_to_vals(pop_pit_dict)
            pos_dxy_vals = tf.dict_to_vals(pop_dxy_dict)

            # print(pos, "1", tracemalloc.get_traced_memory())
            window_dict = tf.stats_to_windows(window_dict, pos, max_pos, pos_piw_vals, pos_pit_vals, pos_dxy_vals)
            # print(pos, "2", tracemalloc.get_traced_memory())

            # print(pos_piw_vals)
            # print(pos_pit_vals)
            # print(pos_dxy_vals)

            print(pos, pos_piw_vals, pos_pit_vals)

            # print("stats loop completed")
            # print(tracemalloc.get_traced_memory())

            # print("pos loop " + pos + " completed. Moving to next position")

piw_headers = ["piw_" + name for name in pop_names]
pit_headers = ["_".join(pair) for pair in list(itertools.combinations(pop_names, 2))]
pit_headers = ["piT_" + name for name in pit_headers]
dxy_headers = ["_".join(pair) for pair in list(itertools.combinations(pop_names, 2))]
dxy_headers = ["dXY_" + name for name in dxy_headers]

# print("scaff" + "," + "window_start" + "," + "window_end" + "," + "n_window_sites" + "," +
#       ",".join(piw_headers) + "," + ",".join(pit_headers) + "," + ",".join(dxy_headers))

# print(tracemalloc.get_traced_memory())

for key in window_dict.keys():
    # need to take mean of each column for each window, for piw and dxy
    window_piw_vals = window_dict[key][0]
    window_pit_vals = window_dict[key][1]
    window_dxy_vals = window_dict[key][2]

    if len(window_piw_vals) > 0 and len(window_pit_vals) > 0 and len(window_dxy_vals) > 0:
        window_piw_means = tf.vals_to_pop_means(window_piw_vals)
        window_pit_means = tf.vals_to_pop_means(window_pit_vals)
        window_dxy_means = tf.vals_to_pop_means(window_dxy_vals)

        # convert to string and print
        window_piw_means = list(map(str, window_piw_means))
        window_pit_means = list(map(str, window_pit_means))
        window_dxy_means = list(map(str, window_dxy_means))
        # print(scaff + "," + str(min(key)) + "," + str(max(key)) + "," + str(len(window_piw_vals)) + "," +
        #       ",".join(window_piw_means) + "," + ",".join(window_pit_means) + "," + ",".join(window_dxy_means))
