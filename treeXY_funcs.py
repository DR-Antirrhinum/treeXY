import random
import itertools
import tracemalloc
from argparse import ArgumentParser
import math as maths

##########################
# command line arguments #
##########################
parser = ArgumentParser(prog="treeXY",
                        description="Calculate pi-statistics from mapped pool-seq data")

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


#############
# functions #
#############
def initialise_windows():
    # windows
    # check the total number of windows
    with open(args.file) as f:
        for i in f:
            # strip trailing newline
            i = i.strip("\n")
            i = i.split("\t")
            # genomic position
            p = i[1]
        max_pos = int(p)

    window_size = int(args.window_size)
    window_overlap = int(args.window_overlap)
    window_slide = window_size - window_overlap

    # initialise window_dict of correct length
    tot_windows = maths.floor(max_pos / window_slide)

    windows = {}
    for i in range(0, tot_windows):
        start_coord = i * window_slide
        end_coord = start_coord + window_size + 1
        if end_coord > max_pos:
            end_coord = max_pos + 1
        # three empty lists, to be populated with pi-w, piT, and dXY
        windows[range(start_coord, end_coord)] = [[], [], []]

    # print("initialise_windows", tracemalloc.get_traced_memory())

    return [windows, max_pos]


def get_sync_counts(sync_site_counts):
    # split all count columns and record counts as integers
    # order in SYNC files is A, T, C, G, N, del

    sync_count_list = []
    for count in sync_site_counts:
        count = count.split(":")
        count = list(map(int, count))
        # n_count = count[4]
        # non_n_count = count[0:4] + count[5:]
        sync_count_list.append(count)

    # print(pos, "get_sync_counts", tracemalloc.get_traced_memory())

    return sync_count_list


def remove_n_count(count):
    # remove N count from list of counts

    n_count = count[4]
    non_n_count = count[0:4] + count[5:]

    # print(pos, "remove_n_count", tracemalloc.get_traced_memory())

    return non_n_count


def remove_low_count(count):
    # remove counts below args threshold from list of counts

    new_counts = []
    for base in count:
        if base < args.min_allele_depth:
            new_counts.append(0)
        else:
            new_counts.append(base)

    # print(pos, "remove_n_count", tracemalloc.get_traced_memory())

    return new_counts


def check_read_depth(sync_count_list):
    # return list: 1 if pop >= threshold depth, 0 otherwise

    depth_list = []
    for count in sync_count_list:
        # remove bases below args threshold
        count = remove_low_count(count)
        count = remove_n_count(count)
        if args.min_depth <= sum(count) <= args.max_depth:
            depth_list.append(1)
        else:
            depth_list.append(0)

    # print(pos, "check_read_depth", tracemalloc.get_traced_memory())

    return depth_list


def check_allele_num(sync_count_list, dpth_pass_list):
    # return number of alleles
    # **** alleles should be returned in descending order of frequency (across all pops)

    # remove pops below threshold depth
    passing_dpth = [i for i, e in enumerate(dpth_pass_list) if e > 0]
    filtered_count_inds = [[i, sync_count_list[i]] for i in passing_dpth]
    filtered_count_list = [i[1] for i in filtered_count_inds]

    if len(filtered_count_inds) != 0:
        allele_inds = []
        for i, count in filtered_count_inds:
            # increment with each allele > min_allele_depth
            allele_incr = 0
            for base_i, base in enumerate(count):
                if base >= args.min_allele_depth:
                    allele_incr += 1
                    allele_inds.append(base_i)

            dpth_pass_list[i] = allele_incr
        allele_inds = list(set(allele_inds))

        # sort alleles so major allele is first in list
        tot_allele_counts = [sum(x) for x in zip(*filtered_count_list)]
        max_allele_ind = [i for i, j in enumerate(tot_allele_counts) if j == max(tot_allele_counts)]
        curr_allele_ind = [i for i, e in enumerate(allele_inds) if e == max_allele_ind[0]]
        allele_inds.insert(0, allele_inds.pop(curr_allele_ind[0]))

        return [dpth_pass_list, allele_inds]

    else:
        return [dpth_pass_list, []]


def filter_triallelic(sync_count_list, scaff, pos, ref):
    # remove least common allele from triallelic sites
    site_tots = [sum(x) for x in zip(*sync_count_list)]

    # set value of N count to 0, so never called as allele
    site_tots[4] = 0

    target = min([i for i in site_tots if i != 0])
    index = site_tots.index(target)
    edited_list = []
    for count in sync_count_list:
        count[index] = 0
        count = list(map(str, count))
        count = ":".join(count)
        edited_list.append(count)

    edited_list.insert(0, scaff)
    edited_list.insert(1, pos)
    edited_list.insert(2, ref)

    edited_list = "\t".join(edited_list)

    # print(pos, "filter_triallelic", tracemalloc.get_traced_memory())

    return edited_list


def get_allele_freqs(dpth_pass_list, allele_inds, sync_count_list, names):
    # *** p should be major allele i.e. most frequent

    pop_dict = {}
    for i, count in enumerate(sync_count_list):
        curr_pop = names[i]
        a1_ind = allele_inds[0]
        a1_count = count[a1_ind]

        if dpth_pass_list[i] == 2:
            a2_ind = allele_inds[1]
            a2_count = count[a2_ind]
        else:
            a2_count = 0

        if a1_count == 0 & a2_count == 0:
            p = 0
            q = 0
            pop_dict[curr_pop] = [p, q]
        else:
            p = a1_count / (a1_count + a2_count)
            q = a2_count / (a1_count + a2_count)
            pop_dict[curr_pop] = [p, q]

    return pop_dict


# calculate Nei's dXY between two populations
# corresponds to raw dXY in David's SlidingWindows program
def get_dxy(p1, p2):
    dxy = (p1 * (1 - p2)) + (p2 * (1 - p1))

    # print("get_dxy", tracemalloc.get_traced_memory())

    return dxy


# calculate pi-within for a given population
# **** implement adjustments for binomial sampling (See SW documentation)
def get_piw(p1):
    piw = 2 * p1 * (1 - p1)

    # print("get_piw", tracemalloc.get_traced_memory())

    return piw


# calculate piT for a given population pair
def get_pit(p1, p2, q1, q2):
    p_bar = (p1 + p2) / 2
    q_bar = (q1 + q2) / 2
    pit = 2 * p_bar * q_bar

    # print("get_pit", tracemalloc.get_traced_memory())

    return pit


# calculate piw across list of valid pops
def get_all_pop_piw(pop_names, dpth_pass_pops, freqs_dict):
    pop_piw_vals = []
    for pop in dpth_pass_pops:
        pop = pop_names[pop]
        pop_pq = freqs_dict[pop]
        pop_piw = get_piw(pop_pq[0])
        pop_piw_vals.append(pop_piw)

    # print(pos, "get_all_pop_piw", tracemalloc.get_traced_memory())

    return pop_piw_vals


# calculate stats for valid pops / comps
def get_all_pop_pit_dxy(pop_names, dpth_pass_comps, freqs_dict):
    pop_pit_vals = []
    pop_dxy_vals = []

    for comp in dpth_pass_comps:
        pop1 = pop_names[comp[0]]
        pop2 = pop_names[comp[1]]
        pop1_pq = freqs_dict[pop1]
        pop2_pq = freqs_dict[pop2]
        pops_pit = get_pit(pop1_pq[0], pop2_pq[0], pop1_pq[1], pop2_pq[1])
        pops_dxy = get_dxy(pop1_pq[0], pop2_pq[0])
        pop_pit_vals.append(pops_pit)
        pop_dxy_vals.append(pops_dxy)

    # print(pos, "get_all_pop_pit_dxy", tracemalloc.get_traced_memory())

    return [pop_pit_vals, pop_dxy_vals]


def get_site_stats(alleles, count_list, pop_names, pop_dpth):
    # calculate p and q for all pops
    freqs_dict = get_allele_freqs(pop_dpth, alleles, count_list, pop_names)
    # make list of valid pops based on indices of pop_dpth
    dpth_pass_pops = [i for i, e in enumerate(pop_dpth) if e != 0]
    # calculate piw for valid pops
    pop_piw_vals = get_all_pop_piw(pop_names, dpth_pass_pops, freqs_dict)

    # make list of valid comparisons based on indices of pop_dpth
    dpth_pass_comps = itertools.combinations(dpth_pass_pops, 2)
    # calculate stats for valid pops / comps
    pairwise_stats = get_all_pop_pit_dxy(pop_names, dpth_pass_comps, freqs_dict)
    pop_pit_vals = pairwise_stats[0]
    pop_dxy_vals = pairwise_stats[1]

    # print(pos, "get_site_stats2", tracemalloc.get_traced_memory())

    return [pop_piw_vals, pop_pit_vals, pop_dxy_vals]


def dict_to_vals(pop_dict):
    # print stats from dict
    val_list = []
    for val in pop_dict.values():
        val_list.append(val)

    # print(pos, "dict_to_vals", tracemalloc.get_traced_memory())

    return val_list


def stats_to_windows(curr_window_dict, curr_pos, w_max_pos, piw, pit, dxy):
    # which windows does this position fall into?
    min_w_start = int(curr_pos) - int(args.window_size)
    if min_w_start < 0:
        min_w_start = 0

    # how many bp does the window move by?
    w_slide = int(args.window_size) - int(args.window_overlap)

    min_w_index = maths.ceil(min_w_start / w_slide)
    # highest poss window start is the last window with a start value < pos
    max_w_index = maths.floor(int(curr_pos) / w_slide)
    for i in range(min_w_index, max_w_index + 1):
        w_start = i * w_slide
        w_end = w_start + int(args.window_size) + 1
        if w_end > w_max_pos:
            w_end = w_max_pos + 1
        if (w_end - w_start) > w_slide:
            range_key = range(w_start, w_end)
            # **** for now, I will append vals, but consider using dicts to retain pop_names
            curr_window_dict[range_key][0].append(piw)
            curr_window_dict[range_key][1].append(pit)
            curr_window_dict[range_key][2].append(dxy)

    # print(pos, "stats_to_windows", tracemalloc.get_traced_memory())

    return curr_window_dict


def vals_to_pop_means(val_list):
    # approach using filter - not very pythonic
    # window_sums = [sum(filter(None, i)) for i in zip(*val_list)]
    # other approach
    mean_list = []
    for col in zip(*val_list):
        # filter None
        col = [x for x in col if x is not None]
        col_len = len(col)
        if col_len > 0:
            col_sum = sum(col)
            col_mean = col_sum / col_len
            mean_list.append(col_mean)

    return mean_list


def get_site_trees():
    pass


def get_site_tree_stats():
    pass
