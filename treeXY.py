import sys
from argparse import ArgumentParser

##########################
# command line arguments #
##########################
parser = ArgumentParser(prog="get_tree_height",
                        description="Find all trees matching a root division pattern")

parser.add_argument("-f", "--file",
                    required=True,
                    help="Input biallelic SYNC file to be processed. ",
                    metavar="sync_file")

parser.add_argument("-p", "--pops",
                    required=False,
                    help="Input pop file describing the pools. "
                         "Column 1 is pop names, column 2 is number in pool, and "
                         "column 3 is 1/2 grouping (or 0 to exclude).",
                    metavar="pop_file")

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

args = parser.parse_args()


#############
# functions #
#############
# split all count columns and record counts as integers
# order in SYNC files is A, T, C, G, N, del
# here, I ignore Ns
def get_sync_counts(sync_site_counts):
    count_list = []
    for count in sync_site_counts:
        count = count.split(":")
        count = list(map(int, count))
        # n_count = count[4]
        # non_n_count = count[0:4] + count[5:]
        count_list.append(count)

    return count_list


def remove_n_count(count):
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


# def check_allelic_depth(sync_count_list, dpth_pass_list):
#     # check depth of alleles in individual pops
#     ad_list = []
#     for count in sync_count_list:
#         allele_depth = [count[i] for i in alleles]
#         alleles_present = [i for i in allele_depth if i > args.min_allele_depth]
#         ad_list.append(len(alleles_present))
#
#     # site is biallelic?
#     biallelic_count = [True for i in ad_list if i == 2]
#     biallelic = any(biallelic_count)
#
#     # site is multiallelic?
#     multiallelic_count = [True for i in ad_list if i > 2]
#     multiallelic = any(multiallelic_count)
#
#     return [ad_list, biallelic, multiallelic]


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


########
# main #
########
with open(args.file) as f:
    for line in f:
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

        # if biallelic, proceed with calculations
        if len(alleles) == 2:
            print(line)
