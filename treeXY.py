import random
import itertools
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

    window_dict = {}
    for i in range(0, tot_windows):
        start_coord = i * window_slide
        end_coord = start_coord + window_size + 1
        if end_coord > max_pos:
            end_coord = max_pos + 1
        window_dict[range(start_coord, end_coord)] = [[], [], []]

    return window_dict


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

    return sync_count_list


def remove_n_count(count):
    # remove N count from list of counts

    n_count = count[4]
    non_n_count = count[0:4] + count[5:]

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

        return alleles

    else:
        return []


def check_pop_alleles(sync_count_list, dpth_pass_list, called_alleles, expected_alleles):
    # check depth of alleles in individual pops
    # **** THIS IS REALLY BAD

    # remove pops below threshold depth
    passing_dpth = [i for i, e in enumerate(dpth_pass_list) if e == 1]
    filtered_count_list = [sync_count_list[i] for i in passing_dpth]

    if len(filtered_count_list) != 0:
        ad_list = []
        for count in filtered_count_list:
            ad_list.append([i for i in called_alleles if count[i] >= args.min_allele_depth])

        ad_len_list = []
        for i in ad_list:
            ad_len_list.append(len(i))

        if ad_len_list.count(expected_alleles) >= args.min_allele_pops:
            return True
        else:
            return False


def filter_triallelic(sync_count_list):
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

    return edited_list


def get_allele_freqs(allele_inds, counts, names):
    pop_dict = {}
    # randomise which allele is which
    # conventionally, p is often designated as the major allele
    # however, we tend to assign it randomly
    # **** specify from command line whether p = major allele or random allele 1
    random.shuffle(allele_inds)
    if len(allele_inds) == 2:
        # what is the index of p and q in the list of alleles?
        p_ind = allele_inds[0]
        q_ind = allele_inds[1]
        # calculate p and q for each population
        # only p is used in calculating dXY
        for ecounts in enumerate(counts):
            count_ind = ecounts[0]
            counts = ecounts[1]
            curr_pop = names[count_ind]
            p_count = counts[p_ind]
            q_count = counts[q_ind]
            if p_count == 0 & q_count == 0:
                p = 0
                q = 0
                pop_dict[curr_pop] = [p, q]
            else:
                p = p_count / (p_count + q_count)
                q = q_count / (p_count + q_count)
                pop_dict[curr_pop] = [p, q]

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

    return pop_dict


# calculate Nei's dXY between two populations
# corresponds to raw dXY in David's SlidingWindows program
def get_dxy(p1, p2):
    dxy = (p1 * (1 - p2)) + (p2 * (1 - p1))
    return dxy


# calculate pi-within for a given population
# **** implement adjustments for binomial sampling (See SW documentation)
def get_piw(p1):
    piw = 2 * p1 * (1 - p1)
    return piw


# calculate piT for a given population pair
def get_pit(p1, p2, q1, q2):
    p_bar = (p1 + p2) / 2
    q_bar = (q1 + q2) / 2
    pit = 2 * p_bar * q_bar
    return pit


def get_site_stats():
    pass


def get_site_trees():
    pass


def get_site_tree_stats():
    pass


##########
# treeXY #
##########
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
        count_list = get_sync_counts(site_counts)
        # **** eventually, want to read names from file
        pop_names = [str(i) for i in range(1, len(count_list) + 1)]

        # read depth checks
        dpth_dict = {}
        pop_dpth = check_read_depth(count_list)
        # retain for later window averaging
        dpth_dict[pos] = pop_dpth

        # get alleles
        alleles = check_allele_num(count_list, pop_dpth)

        while len(alleles) > 2:
            line = filter_triallelic(count_list)
            # strip trailing newline
            line = line.strip("\n")
            line = line.split("\t")
            # all colon delimited count columns
            site_counts = line[3:]
            count_list = get_sync_counts(site_counts)
            # get alleles
            alleles = check_allele_num(count_list, pop_dpth)

        # dict to store all stats for each position
        pos_stats_dict = {}

        # if biallelic, proceed with piT / dXY / tree calculations
        # **** it is only here that I start to ignore lines below read depth threshold, does that make sense?
        # **** I could provide option to output filtered SYNC file here
        if len(alleles) == 2:
            if check_pop_alleles(count_list, pop_dpth, alleles, 2):
                # initialise dict key
                pos_stats_dict[pos] = []
                # calculate p and q for all pops
                freqs_dict = get_allele_freqs(alleles, count_list, pop_names)
                # make list of valid pops based on indices of pop_dpth
                dpth_pass_pops = [i for i, e in enumerate(pop_dpth) if e == 1]
                # calculate piw for valid pops
                pop_piw_dict = {}
                for pop in dpth_pass_pops:
                    pop = pop_names[pop]
                    pop_pq = freqs_dict[pop]
                    pop_piw = get_piw(pop_pq[0])
                    pop_piw_dict["piw_" + pop] = pop_piw

                pos_stats_dict[pos].append(pop_piw_dict)

                # make list of valid comparisons based on indices of pop_dpth
                dpth_pass_comps = itertools.combinations(dpth_pass_pops, 2)
                # calculate stats for valid pops / comps
                pop_pit_dict = {}
                pop_dxy_dict = {}
                for comp in dpth_pass_comps:
                    pop1 = pop_names[comp[0]]
                    pop2 = pop_names[comp[1]]
                    pop1_pq = freqs_dict[pop1]
                    pop2_pq = freqs_dict[pop2]
                    pops_pit = get_pit(pop1_pq[0], pop2_pq[0], pop1_pq[1], pop2_pq[1])
                    pops_dxy = get_dxy(pop1_pq[0], pop2_pq[0])
                    pop_pit_dict["piT_" + pop] = pops_pit
                    pop_dxy_dict["dXY_" + pop] = pops_dxy

                pos_stats_dict[pos].append(pop_pit_dict)
                pos_stats_dict[pos].append(pop_dxy_dict)

        # if monoallelic, calculate piw
        # **** piT and dXY will always evaluate to zero, but need to include in output anyway
        if len(alleles) == 1:
            if check_pop_alleles(count_list, pop_dpth, alleles, 1):
                # initialise dict key
                pos_stats_dict[pos] = []
                # calculate p and q for all pops
                freqs_dict = get_allele_freqs(alleles, count_list, pop_names)
                # make list of valid pops based on indices of pop_dpth
                dpth_pass_pops = [i for i, e in enumerate(pop_dpth) if e == 1]
                # calculate piw for valid pops
                pop_piw_dict = {}
                for pop in dpth_pass_pops:
                    pop = pop_names[pop]
                    pop_pq = freqs_dict[pop]
                    pop_piw = get_piw(pop_pq[0])
                    pop_piw_dict["piw_" + pop] = pop_piw

                # make list of valid comparisons based on indices of pop_dpth
                dpth_pass_comps = itertools.combinations(dpth_pass_pops, 2)
                # calculate stats for valid pops / comps
                pop_pit_dict = {}
                pop_dxy_dict = {}
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

                # compile stats
                pos_stats_dict[pos].append(pop_piw_dict)
                pos_stats_dict[pos].append(pop_pit_dict)
                pos_stats_dict[pos].append(pop_dxy_dict)
