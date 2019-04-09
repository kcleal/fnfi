
from collections import OrderedDict
import click
import data_io
from sortedcontainers import SortedDict


def echo(*args):
    click.echo(args, err=True)


def merge_intervals(intervals, srt=True, pad=0):
    """
    >>> merge_intervals( [('chr1', 1, 4), ('chr1', 2, 5), ('chr2', 3, 5)] )
    >>> [['chr1', 1, 5], ['chr2', 3, 5]]
    """
    if srt:
        sorted_by_lower_bound = sorted(intervals, key=lambda tup: (tup[0], tup[1]))  # by chrom, start, end
    else:
        sorted_by_lower_bound = intervals

    if pad:
        sorted_by_lower_bound = [(c, 0 if i - pad < 0 else i - pad, j + pad) for c, i, j in intervals]

    merged = []
    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
            continue
        elif higher[0] != merged[-1][0]:  # Dont merge intervals on different chroms
            merged.append(higher)
        else:
            lower = merged[-1]  # Last item on merged (end of interval)
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[1] <= lower[2]:
                merged[-1] = (lower[0], lower[1], max(higher[2], lower[2]))
            else:
                merged.append(higher)
    return merged


def scan_whole_genome(inputbam, max_cov, tree):
    depth_d = {}
    approx_read_length = []
    for a in inputbam.fetch():

        # fails quality checks, PCR duplicate, unmapped, mate unmapped, not primary, proper pair, supplementary
        if a.flag & 3854:
            continue

        # Skip if both positions are non-canonical chromosomes
        rname = inputbam.get_reference_name(a.rname)

        bin_pos = (int(a.pos) / 100) * 100
        # Use an ordered dict, obviates need for sorting when merging intervals
        if rname not in depth_d:
            depth_d[rname] = OrderedDict()

        if bin_pos not in depth_d[rname]:
            depth_d[rname][bin_pos] = 0

        depth_d[rname][bin_pos] += 1
        if len(approx_read_length) < 1000:
            approx_read_length.append(len(a.seq))

    approx_read_length = float(sum(approx_read_length)) / len(approx_read_length)

    total_dropped, reads_dropped = 0, 0
    intervals = []
    for chrom, v in depth_d.iteritems():

        for k2, v2 in v.items():
            bin_cov = (v2 * approx_read_length) / 100
            if bin_cov < max_cov or data_io.intersecter(tree, chrom, k2, k2 + 100):
                intervals.append((chrom, k2, k2 + 201))  # Will cause a 1 bp overlap with adjacent window
            else:
                total_dropped += 200
                reads_dropped += v2

    click.echo("Skipping {} kb of reference, primary read-coverage >= {} bp - {} reads".format(
        round(total_dropped / 1e3, 3), max_cov, reads_dropped), err=True)
    merged = merge_intervals(intervals, srt=False, pad=1000)

    return depth_d, approx_read_length, merged


def scan_regions(inputbam, include, max_cov, tree):
    # Scan the input regions, then re-scan the locations of the mate-pairs around the genome
    depth_d = {}

    roi = data_io.get_include_reads(include, inputbam)
    approx_read_length = []

    target_regions = set([])
    intervals_to_check = set([])
    for a in roi:

        if a.flag & 1540:
            # Unmapped, fails Q, duplicate
            continue

        # Skip if both positions are non-canonical chromosomes
        rname, rnext = inputbam.get_reference_name(a.rname), inputbam.get_reference_name(a.rnext)

        bin_pos1 = (int(a.pos) / 100) * 100
        bin_pos2 = (int(a.pnext) / 100) * 100

        # Need to check 10kb around mate, find all intervals to check
        c1 = (rname, 0 if bin_pos1 - 10000 < 0 else bin_pos1 - 10000, bin_pos1 + 10000)
        c2 = (rnext, 0 if bin_pos2 - 10000 < 0 else bin_pos2 - 10000, bin_pos2 + 10000)
        if c1 not in intervals_to_check:
            intervals_to_check.add(c1)
        if c2 not in intervals_to_check:
            intervals_to_check.add(c2)

        if a.flag & 2314:
            # not primary, mate unmapped, supplementary, proper-pair
            continue

        sideA = (rname, bin_pos1)
        sideB = (rnext, bin_pos2)

        if sideA not in target_regions:
            target_regions.add(sideA)
        if sideB not in target_regions:
            target_regions.add(sideB)

        if len(approx_read_length) < 1000:
            if not a.flag & 2048:
                approx_read_length.append(len(a.seq))

    # Get coverage within 10kb
    merged_to_check = merge_intervals(intervals_to_check, srt=True)

    for item in merged_to_check:
        for a in inputbam.fetch(*item):

            if a.flag & 1540:
                # Unmapped, fails Q, duplicate
                continue

            rname = item[0]
            bin_pos1 = (int(a.pos) / 100) * 100

            # Use an ordered dict, obviates need for sorting when merging intervals later on
            if rname not in depth_d:
                depth_d[rname] = SortedDict()
            if bin_pos1 not in depth_d[rname]:
                depth_d[rname][bin_pos1] = 1
            else:
                depth_d[rname][bin_pos1] += 1

    approx_read_length = float(sum(approx_read_length)) / len(approx_read_length)
    click.echo("Read length {}".format(int(approx_read_length)), err=True)

    total_dropped, reads_dropped = 0, 0
    intervals = []
    for chrom, bin_start in target_regions:
        d = depth_d[chrom][bin_start]
        bin_cov = (d * approx_read_length) / 100

        if bin_cov < max_cov or data_io.intersecter(tree, chrom, bin_start, bin_start + 100):
            intervals.append((chrom, bin_start, bin_start + 201))  # Will cause a 1 bp overlap with adjacent window
        else:
            total_dropped += 200
            reads_dropped += d

    click.echo("Skipping {} kb of reference, read-coverage >= {} bp - {} reads".format(
        round(total_dropped / 1e3, 3), max_cov, reads_dropped), err=True)

    merged_targets = merge_intervals(intervals, srt=True)

    return depth_d, approx_read_length, merged_targets


def get_low_coverage_regions(inputbam, max_cov, tree, include):
    # depth_d contains read counts +- 10kb of regions of interest, merged is used to construct the main graph
    if tree is None:
        depth_d, approx_read_length, merged = scan_whole_genome(inputbam, max_cov, tree)

    else:
        depth_d, approx_read_length, merged = scan_regions(inputbam, include, max_cov, tree)

    return merged, approx_read_length, depth_d
