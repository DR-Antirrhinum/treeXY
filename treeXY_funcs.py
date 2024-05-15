import itertools
# import tracemalloc
import math as maths
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as shc
# from scipy.spatial import distance_matrix
from scipy.spatial.distance import squareform


#############
# functions #
#############
def read_pop_names(in_file):
    with open(in_file) as f:
        for i in f:
            # strip trailing newline
            i = i.strip("\n")
            i = i.split("\t")
            # genomic position
            sc = i[3:]
            pop_names = [str(i) for i in range(1, len(sc) + 1)]

            return pop_names


def initialise_windows(in_file, w_size, w_overlap):
    # windows
    # check the total number of windows
    with open(in_file) as f:
        for i, e in enumerate(f):
            print(i, e)
            # strip trailing newline
            e = e.strip("\n")
            e = e.split("\t")
            # genomic position
            p = e[1]
            if i == 0:
                min_pos = int(e[1])
        max_pos = int(p)

    window_size = int(w_size)
    window_overlap = int(w_overlap)
    window_slide = window_size - window_overlap

    # if window_size >= max_pos, return one window covering entire scaffold
    wc_window = {}
    if window_size >= max_pos:
        wc_window[range(1, max_pos + 1)] = [[], [], [], []]

        return [wc_window, max_pos]

    # initialise window_dict of correct length
    tot_windows = maths.floor(max_pos / window_slide)

    windows = {}
    for i in range(0, tot_windows):
        start_coord = i * window_slide
        end_coord = start_coord + window_size + 1
        if end_coord > max_pos:
            end_coord = max_pos + 1
        # three empty lists, to be populated with pi-w, piT, and dXY
        windows[range(start_coord, end_coord)] = [[], [], [], []]

    # remove windows where end coord is < min_pos
    # **** is it worth adding a loop to remove all redundant windows?
    # **** e.g. could record all positions in SYNC file and check against all window ranges
    s_keys = []
    for key in windows.keys():
        if max(key) < min_pos:
            s_keys.append(key)

    for key in s_keys:
        del windows[key]

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


def remove_low_count(count, min_allele_depth):
    # remove counts below args threshold from list of counts

    new_counts = []
    for base in count:
        if base < min_allele_depth:
            new_counts.append(0)
        else:
            new_counts.append(base)

    # print(pos, "remove_n_count", tracemalloc.get_traced_memory())

    return new_counts


def check_read_depth(sync_count_list, min_allele_depth, min_depth, max_depth):
    # return list: 1 if pop >= threshold depth, 0 otherwise

    depth_list = []
    for count in sync_count_list:
        # remove bases below args threshold
        count = remove_low_count(count, min_allele_depth)
        count = remove_n_count(count)
        if min_depth <= sum(count) <= max_depth:
            depth_list.append(1)
        else:
            depth_list.append(0)

    # print(pos, "check_read_depth", tracemalloc.get_traced_memory())

    return depth_list


def check_allele_num(sync_count_list, dpth_pass_list, min_allele_depth):
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
                if base >= min_allele_depth:
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


# def get_formatted_line(scaffold, position, ref_allele, count_list):
#     pass


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


def filter_low_depth(sync_count_list, ad_cutoff):
    new_counts = []
    for pop_counts in sync_count_list:
        for i, count in enumerate(pop_counts):
            if count < ad_cutoff:
                pop_counts[i] = 0
        new_counts.append(pop_counts)

    return new_counts


def get_allele_freqs(allele_inds, sync_count_list, names):
    # *** p should be major allele i.e. most frequent

    pop_dict = {}
    for i, count in enumerate(sync_count_list):
        curr_pop = names[i]
        if len(allele_inds) == 1:
            a1_ind = allele_inds[0]
            a1_count = count[a1_ind]
            a2_count = 0
        else:
            a1_ind = allele_inds[0]
            a1_count = count[a1_ind]
            a2_ind = allele_inds[1]
            a2_count = count[a2_ind]

        if a1_count == 0 and a2_count == 0:
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


# calculate Nei's D (da) for a given population pair
def get_da(p1, p2, dxy):
    piw1 = get_piw(p1)
    piw2 = get_piw(p2)
    piw_bar = (piw1 + piw2) / 2
    da = dxy - piw_bar

    # print("get_pit", tracemalloc.get_traced_memory())

    return da


# calculate FST for a given population pair
def get_fst(p1, p2, dxy):
    piw1 = get_piw(p1)
    piw2 = get_piw(p2)
    piw_bar = (piw1 + piw2) / 2
    da = dxy - piw_bar
    pi_tot = dxy + piw_bar
    fst = da / pi_tot

    # print("get_pit", tracemalloc.get_traced_memory())

    return fst


# calculate piw across list of valid pops
# **** I've removed None padding because I'm a pillock
def get_all_pop_piw(pop_names, dpth_pass_pops, freqs_dict):
    pop_piw_vals = [None] * len(pop_names)
    for pop in dpth_pass_pops:
        pop = pop_names[pop]
        pop_index = pop_names.index(pop)
        pop_pq = freqs_dict[pop]
        pop_piw = get_piw(pop_pq[0])
        pop_piw_vals[pop_index] = pop_piw

    # print(pos, "get_all_pop_piw", tracemalloc.get_traced_memory())

    return pop_piw_vals


# calculate stats for valid pops / comps
def get_all_pop_pit_dxy(pop_names, dpth_pass_comps, freqs_dict):
    # dpth_pass_comps = list(dpth_pass_comps)
    pop_pit_vals = [None] * len(dpth_pass_comps)
    pop_dxy_vals = [None] * len(dpth_pass_comps)
    pop_D_vals = [None] * len(dpth_pass_comps)
    pop_fst_vals = [None] * len(dpth_pass_comps)

    for comp in dpth_pass_comps:
        comp_index = dpth_pass_comps.index(comp)
        pop1 = pop_names[comp[0]]
        pop2 = pop_names[comp[1]]
        pop1_pq = freqs_dict[pop1]
        pop2_pq = freqs_dict[pop2]
        pops_pit = get_pit(pop1_pq[0], pop2_pq[0], pop1_pq[1], pop2_pq[1])
        pops_dxy = get_dxy(pop1_pq[0], pop2_pq[0])
        pops_D = get_da(pop1_pq[0], pop2_pq[0], pops_dxy)
        pop_pit_vals[comp_index] = pops_pit
        pop_dxy_vals[comp_index] = pops_dxy
        pop_D_vals[comp_index] = pops_D

    # print(pos, "get_all_pop_pit_dxy", tracemalloc.get_traced_memory())

    return [pop_pit_vals, pop_dxy_vals, pop_D_vals, pop_fst_vals]


def get_all_pop_fst(pop_names, dpth_pass_comps, freqs_dict):
    # dpth_pass_comps = list(dpth_pass_comps)
    pop_fst_vals = [None] * len(dpth_pass_comps)

    for comp in dpth_pass_comps:
        comp_index = dpth_pass_comps.index(comp)
        pop1 = pop_names[comp[0]]
        pop2 = pop_names[comp[1]]
        pop1_pq = freqs_dict[pop1]
        pop2_pq = freqs_dict[pop2]
        pops_dxy = get_dxy(pop1_pq[0], pop2_pq[0])
        if pops_dxy > 0:
            pops_fst = get_fst(pop1_pq[0], pop2_pq[0], pops_dxy)
            pop_fst_vals[comp_index] = pops_fst
        else:
            pop_fst_vals[comp_index] = None

    # print(pos, "get_all_pop_pit_dxy", tracemalloc.get_traced_memory())

    return pop_fst_vals


def get_site_stats(alleles, count_list, pop_names, pop_dpth):
    # calculate p and q for all pops
    freqs_dict = get_allele_freqs(alleles, count_list, pop_names)
    # make list of valid pops based on indices of pop_dpth
    dpth_pass_pops = [i for i, e in enumerate(pop_dpth) if e != 0]
    # calculate piw for valid pops
    pop_piw_vals = get_all_pop_piw(pop_names, dpth_pass_pops, freqs_dict)

    # make list of valid comparisons based on indices of pop_dpth
    dpth_pass_comps = list(itertools.combinations(dpth_pass_pops, 2))
    # calculate stats for valid pops / comps
    pairwise_stats = get_all_pop_pit_dxy(pop_names, dpth_pass_comps, freqs_dict)
    pop_pit_vals = pairwise_stats[0]
    pop_dxy_vals = pairwise_stats[1]
    pop_D_vals = pairwise_stats[2]

    pop_fst_vals = get_all_pop_fst(pop_names, dpth_pass_comps, freqs_dict)

    print(pop_fst_vals)

    # print(pos, "get_site_stats2", tracemalloc.get_traced_memory())

    return [pop_piw_vals, pop_pit_vals, pop_dxy_vals, pop_D_vals]


def dict_to_vals(pop_dict):
    # print stats from dict
    val_list = []
    for val in pop_dict.values():
        val_list.append(val)

    # print(pos, "dict_to_vals", tracemalloc.get_traced_memory())

    return val_list


def stats_to_windows(curr_window_dict, curr_pos, w_max_pos, piw, pit, dxy, D, window_size, window_overlap):
    # if window_size >= w_max_pos, there will only be one key encompassing the whole scaffold
    if window_size >= w_max_pos:
        curr_window_dict[range(1, w_max_pos + 1)][0].append(piw)
        curr_window_dict[range(1, w_max_pos + 1)][1].append(pit)
        curr_window_dict[range(1, w_max_pos + 1)][2].append(dxy)
        curr_window_dict[range(1, w_max_pos + 1)][3].append(D)

        return curr_window_dict

    else:
        # which windows does this position fall into?
        min_w_start = int(curr_pos) - int(window_size)
        if min_w_start < 0:
            min_w_start = 0

        # how many bp does the window move by?
        w_slide = int(window_size) - int(window_overlap)

        min_w_index = maths.ceil(min_w_start / w_slide)
        # highest poss window start is the last window with a start value < pos
        max_w_index = maths.floor(int(curr_pos) / w_slide)
        for i in range(min_w_index, max_w_index + 1):
            w_start = i * w_slide
            w_end = w_start + int(window_size) + 1
            if w_end > w_max_pos:
                w_end = w_max_pos + 1
            if (w_end - w_start) > w_slide:
                range_key = range(w_start, w_end)
                # **** for now, I will append vals, but consider using dicts to retain pop_names
                curr_window_dict[range_key][0].append(piw)
                curr_window_dict[range_key][1].append(pit)
                curr_window_dict[range_key][2].append(dxy)
                curr_window_dict[range_key][3].append(D)

        # print(pos, "stats_to_windows", tracemalloc.get_traced_memory())

        return curr_window_dict


def average_window(key_to_avg, curr_window_dict):
    # take mean of each column for each window, for piw and dxy
    window_piw_vals = curr_window_dict[key_to_avg][0]
    window_pit_vals = curr_window_dict[key_to_avg][1]
    window_dxy_vals = curr_window_dict[key_to_avg][2]
    window_da_vals = curr_window_dict[key_to_avg][3]

    n_sites = str(len(window_piw_vals))

    if len(window_piw_vals) > 0 and len(window_pit_vals) > 0 and len(window_dxy_vals) > 0:
        window_piw_means = vals_to_pop_means(window_piw_vals)
        window_pit_means = vals_to_pop_means(window_pit_vals)
        window_dxy_means = vals_to_pop_means(window_dxy_vals)
        window_da_means = vals_to_pop_means(window_da_vals)

        # convert to string and write to file
        window_piw_means = list(map(str, window_piw_means))
        window_pit_means = list(map(str, window_pit_means))
        window_dxy_means = list(map(str, window_dxy_means))
        window_da_means = list(map(str, window_da_means))

        return [n_sites, window_piw_means, window_pit_means, window_dxy_means, window_da_means]

    # return False if window has no associated values
    else:
        return False


def write_window(scaff, key, w_avg, w_out_file):
    with open(w_out_file, "a") as out_file:
        n_window_sites = w_avg[0]
        window_piw_means = w_avg[1]
        window_pit_means = w_avg[2]
        window_dxy_means = w_avg[3]
        window_da_means = w_avg[4]

        out_file.write(scaff + "," + str(min(key)) + "," + str(max(key)) + "," +
                       n_window_sites + "," + ",".join(window_piw_means) + "," +
                       ",".join(window_pit_means) + "," +
                       ",".join(window_dxy_means) + "," +
                       ",".join(window_da_means) + "\n")


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


def get_site_trees(pop_dpth, distance_vals):
    pop_pairs = []
    dpth_pass_pops = [i for i, e in enumerate(pop_dpth) if e != 0]
    dxy_headers = itertools.combinations(dpth_pass_pops, 2)
    for i, header in enumerate(dxy_headers):
        pop1 = header[0]
        pop2 = header[1]
        dxy = distance_vals[i]
        pop_pairs.append([pop1, pop2, dxy])

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

    # get maximum height at which UPGMA merges populations i.e. the tree height
    tree_height = max(clusters[0:, 2])
    second_max = clusters[0:, 2][-2]
    # shortest root branch = tree height - penultimate cluster height
    # **** implement piw and Nei's D
    srb = tree_height - second_max

    # split tree at root, yielding two clusters (root division analysis)
    # population membership to each cluster is represented as 0 or 1
    split = list(shc.cut_tree(clusters, n_clusters=2)[0:, 0])
    split = [n + 1 for n in split]

    # convert items in split list to strings, and join
    split = list(map(str, split))
    split = "".join(split)

    # return genomic position, root division summary, and tree height
    return [split, tree_height, srb]


def get_site_tree_stats():
    pass
