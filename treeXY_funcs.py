import random
import itertools
import tracemalloc
from argparse import ArgumentParser
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


def check_read_depth(sync_count_list):
    # return list: 1 if pop >= threshold depth, 0 otherwise

    depth_list = []
    for count in sync_count_list:
        count = remove_n_count(count)
        if args.min_depth <= sum(count) <= args.max_depth:
            depth_list.append(1)
        else:
            depth_list.append(0)

    # print(pos, "check_read_depth", tracemalloc.get_traced_memory())

    return depth_list


def check_allele_num(sync_count_list, dpth_pass_list):
    # return number of alleles
    # **** also return number of pops alleles present in?

    # remove pops below threshold depth
    passing_dpth = [i for i, e in enumerate(dpth_pass_list) if e == 1]
    filtered_count_list = [sync_count_list[i] for i in passing_dpth]

    if len(filtered_count_list) != 0:
        # sum depth across all remaining counts (excluding N)
        count_tots = [sum(x) for x in zip(*filtered_count_list)]
        # set value of N count to 0, so never called as allele
        count_tots[4] = 0
        # which alleles have >= threshold depth?
        alleles = [i for i, e in enumerate(count_tots) if e >= args.min_allele_depth]

    # **** consider adding 1 to dpth_pass value in pops which have the second allele above threshold depth?

        # print(pos, "check_allele_num", tracemalloc.get_traced_memory())

        return alleles

    else:
        # print(pos, "check_allele_num", tracemalloc.get_traced_memory())

        return []


def check_pop_allele_depth(dpth_pass_list, sync_count_list, called_alleles):
    passing_dpth = [i for i, e in enumerate(dpth_pass_list) if e == 1]
    print(sync_count_list[passing_dpth])

    return []


def check_pop_alleles(sync_count_list, dpth_pass_list, called_alleles):
    # check depth of alleles in individual pops
    # **** THIS IS REALLY BAD
    # print(pos, "check_pop_alleles1", tracemalloc.get_traced_memory())
    # remove pops below threshold depth
    passing_dpth = [i for i, e in enumerate(dpth_pass_list) if e == 1]
    # print(pos, "check_pop_alleles2", tracemalloc.get_traced_memory())
    filtered_count_list = [sync_count_list[i] for i in passing_dpth]
    # print(pos, "check_pop_alleles3", tracemalloc.get_traced_memory())

    if len(filtered_count_list) != 0:
        # print(pos, "check_pop_alleles4", tracemalloc.get_traced_memory())
        ad_list = []
        # print(pos, "check_pop_alleles5", tracemalloc.get_traced_memory())
        for count in filtered_count_list:
            # print(pos, "check_pop_alleles6", tracemalloc.get_traced_memory())
            ad_list.append([i for i in called_alleles if count[i] >= args.min_allele_depth])
            # print(pos, "check_pop_alleles7", tracemalloc.get_traced_memory())

        ad_len_list = []
        # print(pos, "check_pop_alleles8", tracemalloc.get_traced_memory())
        for i in ad_list:
            # print(pos, "check_pop_alleles9", tracemalloc.get_traced_memory())
            ad_len_list.append(len(i))
            # print(pos, "check_pop_alleles10", tracemalloc.get_traced_memory())

        expected_alleles = len(called_alleles)

        # does site pass as biallelic?
        if ad_len_list.count(expected_alleles) >= args.min_allele_pops:
            # print(pos, "check_pop_alleles12", tracemalloc.get_traced_memory())
            return True

        # does site pass as monoallelic?
        # if pop has good depth, but does not meet min_allele_pops, should PASS AS MONOALLELIC
        # *** BUT what should I do with the biallelic taxon?
        # **** PROBLEM HERE - need to remove minor allele prior to downstream calculations
        # elif ad_len_list.count(1) >= args.min_allele_pops:
        #     filter_biallelic...
            # sync_count_list = remove_minor_allele(sync_count_list)
            # print(pos, "check_pop_alleles12", tracemalloc.get_traced_memory())
            # return [True, sync_count_list]

        else:
            return False


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


def get_allele_freqs(allele_inds, counts, names):
    # print(pos, "get_allele_freqs", tracemalloc.get_traced_memory())
    pop_dict = {}
    # randomise which allele is which
    # conventionally, p is often designated as the major allele
    # however, we tend to assign it randomly
    # **** specify from command line whether p = major allele or random allele 1
    # random.shuffle(allele_inds)
    # print(pos, "get_allele_freqs1", tracemalloc.get_traced_memory())
    if len(allele_inds) == 2:
        # print(pos, "get_allele_freqs2", tracemalloc.get_traced_memory())
        # what is the index of p and q in the list of alleles?
        p_ind = allele_inds[0]
        # print(pos, "get_allele_freqs3", tracemalloc.get_traced_memory())
        q_ind = allele_inds[1]
        # print(pos, "get_allele_freqs4", tracemalloc.get_traced_memory())
        # calculate p and q for each population
        # only p is used in calculating dXY
        for ecounts in enumerate(counts):
            # print(pos, "get_allele_freqs5", tracemalloc.get_traced_memory())
            count_ind = ecounts[0]
            # print(pos, "get_allele_freqs6", tracemalloc.get_traced_memory())
            counts = ecounts[1]
            # print(pos, "get_allele_freqs7", tracemalloc.get_traced_memory())
            curr_pop = names[count_ind]
            # print(pos, "get_allele_freqs8", tracemalloc.get_traced_memory())
            p_count = counts[p_ind]
            # print(pos, "get_allele_freqs9", tracemalloc.get_traced_memory())
            q_count = counts[q_ind]
            # print(pos, "get_allele_freqs10", tracemalloc.get_traced_memory())
            if p_count == 0 & q_count == 0:
                # print(pos, "get_allele_freqs11", tracemalloc.get_traced_memory())
                p = 0
                # print(pos, "get_allele_freqs12", tracemalloc.get_traced_memory())
                q = 0
                # print(pos, "get_allele_freqs13", tracemalloc.get_traced_memory())
                pop_dict[curr_pop] = [p, q]
                # print(pos, "get_allele_freqs14", tracemalloc.get_traced_memory())
            else:
                p = p_count / (p_count + q_count)
                # print(pos, "get_allele_freqs15", tracemalloc.get_traced_memory())
                q = q_count / (p_count + q_count)
                # print(pos, "get_allele_freqs16", tracemalloc.get_traced_memory())
                pop_dict[curr_pop] = [p, q]
                # print(pos, "get_allele_freqs17", tracemalloc.get_traced_memory())

        return pop_dict

    if len(allele_inds) == 1:
        # what is the index of p in the list of alleles?
        p_ind = allele_inds[0]
        # calculate p and q for each population
        for ecounts in enumerate(counts):
            count_ind = ecounts[0]
            counts = ecounts[1]
            curr_pop = names[count_ind]
            p_count = counts[p_ind]
            if p_count == 0:
                p = 0
                q = 0
                pop_dict[curr_pop] = [p, q]
            else:
                p = 1
                q = 0
                pop_dict[curr_pop] = [p, q]

        # **** memory leak here?
        # print(pos, "get_allele_freqs", tracemalloc.get_traced_memory())

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
    pop_piw_dict = {}
    # create piw keys for each pop
    for pop in pop_names:
        pop_piw_dict["piw_" + pop] = None
    for pop in dpth_pass_pops:
        pop = pop_names[pop]
        pop_pq = freqs_dict[pop]
        pop_piw = get_piw(pop_pq[0])
        pop_key = "piw_" + pop
        pop_piw_dict[pop_key] = pop_piw

    # print(pos, "get_all_pop_piw", tracemalloc.get_traced_memory())

    return pop_piw_dict


# calculate stats for valid pops / comps
def get_all_pop_pit_dxy(pop_names, dpth_pass_comps, freqs_dict):
    pop_pit_dict = {}
    # create piT keys for each pop
    pit_keys = ["piT_" + i + "_" + j for i, j in itertools.combinations(pop_names, 2)]
    for i in pit_keys:
        pop_pit_dict[i] = None

    pop_dxy_dict = {}
    # create piT keys for each pop
    dxy_keys = ["dXY_" + i + "_" + j for i, j in itertools.combinations(pop_names, 2)]
    for i in dxy_keys:
        pop_dxy_dict[i] = None

    for comp in dpth_pass_comps:
        pop1 = pop_names[comp[0]]
        pop2 = pop_names[comp[1]]
        pop_1_2 = pop1 + "_" + pop2
        pop1_pq = freqs_dict[pop1]
        pop2_pq = freqs_dict[pop2]
        pops_pit = get_pit(pop1_pq[0], pop2_pq[0], pop1_pq[1], pop2_pq[1])
        pops_dxy = get_dxy(pop1_pq[0], pop2_pq[0])
        pop_pit_dict["piT_" + pop_1_2] = pops_pit
        pop_dxy_dict["dXY_" + pop_1_2] = pops_dxy

    # print(pos, "get_all_pop_pit_dxy", tracemalloc.get_traced_memory())

    return [pop_pit_dict, pop_dxy_dict]


def get_site_stats(alleles, count_list, pop_names, pop_dpth):
    # initialise dict key
    # pos_stats_dict[pos] = []
    # print(pos, "get_site_stats1", tracemalloc.get_traced_memory())
    # calculate p and q for all pops
    freqs_dict = get_allele_freqs(alleles, count_list, pop_names)
    # make list of valid pops based on indices of pop_dpth
    dpth_pass_pops = [i for i, e in enumerate(pop_dpth) if e == 1]
    # calculate piw for valid pops
    pop_piw_dict = get_all_pop_piw(pop_names, dpth_pass_pops, freqs_dict)

    # make list of valid comparisons based on indices of pop_dpth
    dpth_pass_comps = itertools.combinations(dpth_pass_pops, 2)
    # calculate stats for valid pops / comps
    pairwise_stats = get_all_pop_pit_dxy(pop_names, dpth_pass_comps, freqs_dict)
    pop_pit_dict = pairwise_stats[0]
    pop_dxy_dict = pairwise_stats[1]

    # print(pos, "get_site_stats2", tracemalloc.get_traced_memory())

    return [pop_piw_dict, pop_pit_dict, pop_dxy_dict]


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
