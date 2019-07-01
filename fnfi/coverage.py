from __future__ import absolute_import
from collections import OrderedDict
import click
from . import data_io
from sortedcontainers import SortedDict
import numpy as np


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

    approx_read_length = np.array(approx_read_length).mean()

    total_dropped, reads_dropped = 0, 0
    intervals = []
    for chrom, v in depth_d.items():

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
    inserts = []

    intervals_to_check = set([])  # A list of regions containing many bins
    target_regions = set([])  # A list of individual bins to work on
    for a in roi:

        if a.flag & 1540:
            # Unmapped, fails Q, duplicate
            continue

        if a.flag & 266:
            # not primary, mate unmapped, proper-pair
            continue

        # Skip if both positions are non-canonical chromosomes
        rname, rnext = inputbam.get_reference_name(a.rname), inputbam.get_reference_name(a.rnext)

        bin_pos1 = int((int(a.pos) / 100)) * 100
        bin_pos2 = int((int(a.pnext) / 100)) * 100

        # Need to check 10kb around pair and mate, find all intervals to check (limits where SVs will be called)
        c1 = (rname, 0 if bin_pos1 - 10000 < 0 else bin_pos1 - 10000, bin_pos1 + 10000)
        c2 = (rnext, 0 if bin_pos2 - 10000 < 0 else bin_pos2 - 10000, bin_pos2 + 10000)

        intervals_to_check.add(c1)
        intervals_to_check.add(c2)

        if len(approx_read_length) < 1000:
            if not a.flag & 2048:
                approx_read_length.append(len(a.seq))
                if abs(a.template_length) < 1000:
                    inserts.append(int(abs(a.template_length)))

        # Target regions include upstream and downstream of target - means supplementary reads are not missed later
        target_regions.add((rname, bin_pos1))
        target_regions.add((rnext, bin_pos2))

    # Get coverage within 10kb
    merged_to_check = merge_intervals(intervals_to_check, srt=True)

    # Get coverage in regions
    for item in merged_to_check:
        for a in inputbam.fetch(*item):

            if a.flag & 1540:
                # Unmapped, fails Q, duplicate
                continue

            rname = item[0]
            bin_pos1 = int(int(a.pos) / 100) * 100

            # Get depth at this locus
            # Use an ordered dict, obviates need for sorting when merging intervals later on
            if rname not in depth_d:
                depth_d[rname] = SortedDict()
            if bin_pos1 not in depth_d[rname]:
                depth_d[rname][bin_pos1] = 1
            else:
                depth_d[rname][bin_pos1] += 1

    approx_read_length = np.array(approx_read_length)
    if len(inserts) > 0:
        insert_std = np.array(inserts).std()
    else:
        insert_std = 600
    approx_read_length = approx_read_length.mean()

    # Increase target region space here. Currently target regions only point to primary alignments, increase to capture
    # supplementary mappings nearby
    pad = insert_std * 3

    new_targets = set([])
    for chrom, primary_site in target_regions:
        # Pad is determined by insert stdev
        lower_bin = int((primary_site - pad) / 100) * 100
        lower_bin = 0 if lower_bin < 0 else lower_bin
        upper_bin = (int((primary_site + pad) / 100) * 100) + 100
        for s in range(lower_bin, upper_bin, 100):
            new_targets.add((chrom, s))

    total_dropped, reads_dropped = 0, 0
    intervals = []
    for chrom, bin_start in new_targets:

        if bin_start in depth_d[chrom]:
            d = depth_d[chrom][bin_start]
        else:
            d = 0  # No reads found
        bin_cov = (d * approx_read_length) / 100

        if bin_cov < max_cov or data_io.intersecter(tree, chrom, bin_start, bin_start + 100):
            intervals.append((chrom, int(bin_start), int(bin_start + 201)))  # 1 bp overlap with adjacent window
        else:
            total_dropped += 200
            reads_dropped += d

    click.echo("Skipping {} kb of reference, read-coverage >= {} x - {} reads".format(
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
