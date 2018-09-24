import pysam
from collections import defaultdict, deque
import numpy as np
import networkx as nx
from pybedtools import BedTool
import quicksect
import itertools
import click
import sys
import pandas as pd
from src import assembler, caller
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


class Scoper:
    """Keeps track of which reads are in scope. Maximum distance depends on the template insert_median"""
    def __init__(self, max_dist, clip_length=20):
        self.max_dist = max_dist
        self.scope = deque([])
        self.clip_length = clip_length

    def update(self, current_read):

        if len(self.scope) == 0:
            self.scope.append(current_read)
            return

        elif len(self.scope) > 0 and current_read.rname != self.scope[-1].rname:
            self.scope = deque([current_read])  # Empty scope on new chromosome
            return

        current_pos = current_read.pos

        # Go through reads and check if they are in scope
        while True:
            if len(self.scope) > 0:
                p = self.scope[0].pos
                if current_pos - p < self.max_dist:
                    break
                self.scope.popleft()
                continue
            break
        self.scope.append(current_read)

    @staticmethod
    def overlap(start1, end1, start2, end2):
        return max(0, min(end1, end2) - max(start1, start2))

    def iterate(self):
        if len(self.scope) > 1:
            for i, t in enumerate(self.scope):
                if i == len(self.scope) - 1:
                    break  # Don't pair with self
                if t.qname == self.scope[-1].qname:
                    continue  # Don't pair with same template
                current = self.scope[-1]

                # If current is a supplementary, it MAY be hard-clipped so only look for alignment-overlaps
                if current.flag & 2048 and (current.cigar[0][0] == 5 or current.cigar[-1][0] == 5):

                    # Make sure hard-clip can be overlapped with soft-clip on t
                    if current.cigar[0][0] == 5 and current.cigar[0][1] >= self.clip_length:  # Left HARD-clip current
                        if t.cigar[0][0] == 4 and t.cigar[0][1] >= self.clip_length:  # Left soft-clip on t
                            # Extent of hard and soft-clip overlap
                            c_sc_start = current.pos - current.cigar[0][1]
                            t_sc_start = t.pos - t.cigar[0][1]
                            if self.overlap(t_sc_start, t.pos, c_sc_start, current.pos) >= self.clip_length:
                                # Read, soft-clip overlap primary edge, hard-clip overlap supplementary edge
                                yield t, False, True

                    elif current.cigar[-1][0] == 5 and current.cigar[-1][1] > self.clip_length:  # Right HARD-clip
                        if t.cigar[-1][0] == 4 and t.cigar[-1][1] > self.clip_length:  # Right soft-clip
                            c_sc_end = current.reference_end + current.cigar[-1][1]
                            t_sc_end = t.reference_end + t.cigar[-1][1]
                            if self.overlap(t.reference_end, t_sc_end, current.reference_end,
                                            c_sc_end) >= self.clip_length:
                                yield t, False, True

                #
                if current.cigar[0][0] == 4 and current.cigar[0][1] >= self.clip_length:  # Left soft-clip on current
                    if t.cigar[0][0] == 4:  # and t.cigar[0][1] >= self.clip_length:  # Left soft-clip on t
                        # Extent of sof-clip, make sure-soft-clip portion overlaps not the alignment part
                        c_sc_start = current.pos - current.cigar[0][1]
                        t_sc_start = t.pos - t.cigar[0][1]
                        if self.overlap(t_sc_start, t.pos, c_sc_start, current.pos) >= self.clip_length:
                            yield t, True, False

                    elif t.cigar[0][0] == 5 and t.flag & 2048 and t.cigar[0][1] >= self.clip_length:  # Hard clip on t
                        c_sc_start = current.pos - current.cigar[0][1]
                        t_sc_start = t.pos - t.cigar[0][1]
                        if self.overlap(t_sc_start, t.pos, c_sc_start, current.pos) >= self.clip_length:
                            yield t, False, True

                elif current.cigar[0][0] == 4:
                    if t.cigar[0][0] == 4 and t.cigar[0][1] >= self.clip_length:  # Left soft-clip on t
                        # Extent of sof-clip, make sure-soft-clip portion overlaps not the alignment part
                        c_sc_start = current.pos - current.cigar[0][1]
                        t_sc_start = t.pos - t.cigar[0][1]
                        if self.overlap(t_sc_start, t.pos, c_sc_start, current.pos) >= self.clip_length:
                            yield t, True, False

                elif current.cigar[-1][0] == 4 and current.cigar[-1][1] > self.clip_length:  # Right soft-clip
                    if t.cigar[-1][0] == 4: # and t.cigar[-1][1] > self.clip_length:
                        c_sc_end = current.reference_end + current.cigar[-1][1]
                        t_sc_end = t.reference_end + t.cigar[-1][1]
                        if self.overlap(t.reference_end, t_sc_end, current.reference_end, c_sc_end) >= self.clip_length:
                            yield t, True, False

                    elif t.cigar[-1][0] == 5 and t.flag & 2048 and t.cigar[-1][1] >= self.clip_length:
                        c_sc_end = current.reference_end + current.cigar[-1][1]
                        t_sc_end = t.reference_end + t.cigar[-1][1]
                        if self.overlap(t.reference_end, t_sc_end, current.reference_end, c_sc_end) >= self.clip_length:
                            yield t, False, True

                elif current.cigar[-1][0] == 4:
                    if t.cigar[-1][0] == 4 and t.cigar[-1][1] > self.clip_length:
                        c_sc_end = current.reference_end + current.cigar[-1][1]
                        t_sc_end = t.reference_end + t.cigar[-1][1]
                        if self.overlap(t.reference_end, t_sc_end, current.reference_end, c_sc_end) >= self.clip_length:
                            yield t, True, False

                else:
                    yield t, False, False

        yield None, None, None  # When not enough reads are in scope


def echo(*args):
    click.echo(args, err=True)


def get_regions(reg_path):
    if not reg_path:
        return None
    else:
        f = BedTool(reg_path)
        regions = defaultdict(lambda: quicksect.IntervalTree())
        for item in f:
            regions[item.chrom].add(item.start, item.end)
        return regions


def align_match(r, t):

    def extent(a):
        start = a.pos
        end = a.pos + len(a.seq)
        if a.cigartuples[0][0] == 4:
            start -= a.cigartuples[0][1]
            end -= a.cigartuples[0][1]
        return start, end

    r_start, r_end = extent(r)
    t_start, t_end = extent(t)

    r_offset = r_start - (min(r_start, t_start))
    t_offset = t_start - (min(r_start, t_start))

    size = max(len(t.seq) + t_offset, len(r.seq) + r_offset)
    arr = np.chararray((2, size), itemsize=1)
    arr[:] = ""

    # Convert string to character array
    r_s = np.array(list(r.seq))
    t_s = np.array(list(t.seq))

    arr[0, r_offset: len(r.seq) + r_offset] = r_s
    arr[1, t_offset: len(t.seq) + t_offset] = t_s

    istart = max(r_offset, t_offset)
    iend = min(len(r.seq) + r_offset, len(t.seq) + t_offset)

    mm = arr[0, istart:iend] != arr[1, istart:iend]
    mismatches = np.sum(mm)
    matches = len(mm) - mismatches
    identity = float(matches) / (mismatches + matches)

    p = 1.0
    # Find probability of mismatches
    if mismatches > 0:
        mm_idxs = np.where(mm)[0]
        r_idxs = [i - r_offset + istart for i in mm_idxs]
        t_idxs = [i - t_offset + istart for i in mm_idxs]

        rq = [r.query_qualities[i] for i in r_idxs]
        tq = [t.query_qualities[i] for i in t_idxs]

        minqual = np.min([rq, tq], axis=0).astype(float)

        p = np.multiply.reduce(10 ** (-minqual / 10))

    return identity, p


def get_reads(infile, sub_graph, max_dist, rl):
    nodes_u, nodes_v, edge_data = zip(*sub_graph.edges(data=True))
    nodes = set(nodes_u).union(set(nodes_v))
    coords = []
    for i in edge_data:
        coords += i['p']
    coords = itertools.groupby(sorted(set(coords), key=lambda x: (x[0], x[1])), key=lambda x: x[0])  # Groupby chrom

    # Make intervals of reads that are near one another
    # https://stackoverflow.com/questions/10016802/python-group-a-list-of-integer-with-nearest-values
    c = 0
    rd = defaultdict(lambda: defaultdict(dict))  # rname: flag: alignment
    for chrom, crds in coords:
        d = [i[1] for i in crds]
        chrom = infile.get_reference_name(chrom)
        m = [[chrom, d[0] - rl, d[0] + 2 * rl]]

        for x in d[1:]:
            if x - m[-1][-1] < max_dist:
                m[-1][-1] = x + 2 * rl
            else:
                m.append([chrom, x - rl, x + 2*rl])

        # Collect reads
        for chrom, start, end in m:
            for aln in infile.fetch(chrom, 1 if start <= 0 else start, end):
                if (aln.qname, aln.flag) in nodes:
                    rd[aln.qname][aln.flag] = aln
                    c += 1
    # assert(c == len(nodes))  #!
    if c != len(nodes):
        click.echo("WARNING: reads missing {} vs {}".format(c, len(nodes)), err=True)

    return rd, nodes_u, nodes_v, edge_data


def max_kmer(reads, k=27):
    kmers = defaultdict(int)
    for s in (i.seq for i in reads if not i.flag & 2048):  # Skip supplementary
        if s:
            for j in [s[i:i+k] for i in range(len(s) - k)]:
                kmers[j] += 1
    if len(kmers) > 0:
        return max(kmers.values())
    return 0


def merge_intervals(intervals):
    # thanks https://codereview.stackexchange.com/questions/69242/merging-overlapping-intervals
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []
    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = (lower[0], upper_bound)  # replace by merged interval
            else:
                merged.append(higher)
    return merged


def count_soft_clip_stacks(reads):
    blocks = []
    for r in reads:  # (i for i in reads if not i.flag & 2048):  # Skip supp
        cigar = r.cigartuples
        if cigar:
            if (cigar[0][0] == 4 and cigar[0][1] > 20) or (cigar[-1][0] == 4 and cigar[-1][1] > 20):
                blocks.append((r.pos, r.reference_end))
    mi = merge_intervals(blocks)
    return len(mi)


def get_mate_flag(flag, name, rds):
    if flag & 64:  # First in pair; get second-in-pair primary
        return [i for i in rds[name].keys() if not i & 2048 and i & 128][0]
    else:
        return [i for i in rds[name].keys() if not i & 2048 and i & 64][0]


def get_pri_flag(flag, name, rds):
    if flag & 64:
        return [i for i in rds[name].keys() if not i & 2048 and i & 64][0]
    else:
        return [i for i in rds[name].keys() if not i & 2048 and not i & 64][0]


def process_edge_set(edges, all_reads, bam, insert_size, insert_stdev, get_mate=True):
    # Grey edge means primary alignments share a similar rearrangement pattern. Therefore, the other mate-pair needs to
    # be collected for each node (read) along each edge, make a call for each read-mate pair

    break_points = []
    seen = set([])
    for u, v, d in edges:
        for n, f in (u, v):  # node name, flag
            if n in seen:
                continue
            seen.add(n)
            try:
                if f & 2048:
                    primary1 = all_reads[n][get_pri_flag(f, n, all_reads)]
                else:
                    primary1 = all_reads[n][f]
                if get_mate:
                    mate1 = all_reads[n][get_mate_flag(primary1.flag, n, all_reads)]
                else:
                    mate1 = all_reads[n][f]
            except ValueError:
                continue  # Badly formatted flag, or unpaired

            pair = sorted([primary1, mate1], key=lambda x: x.pos)
            join = []
            for read in pair:
                rid = {read.qname}
                # Sometimes a clip may be present, use this as a break point if available
                cl = []  # Get the biggest clipped portion of the read-pair to find the break-point
                if read.cigartuples[0][0] == 4 or read.cigartuples[0][0] == 5:  # Left clip
                    cl.append((5, rid, bam.get_reference_name(read.rname), read.pos))
                if read.cigartuples[-1][0] == 4 or read.cigartuples[-1][0] == 5:
                    cl.append((3, rid, bam.get_reference_name(read.rname), read.reference_end))
                if len(cl) > 0:
                    join.append(sorted(cl, key=lambda x: x[1])[-1])
                else:
                    # Breakpoint position is beyond the end of the last read
                    if read.flag & 16:  # Read on reverse strand guess to the left
                        p = read.pos - (insert_size / 2) + insert_stdev
                        t = 5
                    else:
                        p = read.reference_end + (insert_size / 2) - insert_stdev
                        t = 3
                    join.append((t, rid, bam.get_reference_name(read.rname), int(p)))
            break_points += join

    return break_points  # caller.call_break_points(break_points)


def score_reads(read_names, all_reads):
    if len(read_names) == 0:
        return {}

    pri, sup = [], []
    for name, flag in read_names:
        for flag in all_reads[name]:
            read = all_reads[name][flag]
            if not read.flag & 2048:
                pri.append(read)
            else:
                sup.append(read)

    data = {"Full_call": True}
    data["max_k_coverage"] = max_kmer(pri + sup)
    data["soct_clip_stacks"] = count_soft_clip_stacks(pri + sup)

    tags_primary = [dict(j.tags) for j in pri]
    tags_supp = [dict(j.tags) for j in sup]

    for item in ["DP", "DA", "NM", "DN", "AS"]:
        v = np.array([i[item] for i in tags_primary]).astype(float)
        if len(v) > 0:
            data[item + "pri"] = v.mean()
        else:
            echo("WARNING: {} tag missing", item)
    data["SP" + "pri"] = len([1 for i in tags_primary if int(i["SP"]) == 1]) / 2.  # Number of split pairs

    for item in ["EV", "DA", "NM", "NP", "AS"]:  # EV depends on mapper used
        if len(tags_supp) == 0:
            data[item] = None
        else:
            if item == "EV" and "EV" in tags_supp[0]:
                data[item + "sup"] = np.array([i[item] for i in tags_supp]).astype(float).min()   # Minimum E-value
            elif item != "EV":
                data[item + "sup"] = np.array([i[item] for i in tags_supp]).astype(float).mean()  # Mean
            else:
                data[item + "sup"] = None

    return data


def cluster_reads(args):

    # Read data into a graph
    # Link reads that share an overlapping soft-clip
    # Weak-link reads that are within ~500bb and have same rearrangement pattern
    # Link read-pairs
    # Link supplementary reads to primary reads that have a soft-clip for checking primary-primary overlap
    # Use a scope to weak-link reads

    infile = pysam.AlignmentFile(args["sv_bam"], "rb")
    regions = get_regions(args["include"])

    edges = []
    all_flags = defaultdict(lambda: defaultdict(list))  # Linking clusters together rname: (flag): (chrom, pos)

    max_dist = args["insert_median"]+(4*args["insert_stdev"])
    scope = Scoper(max_dist=max_dist)

    count = 0
    for r in infile:

        if count % 10000 == 0:
            click.echo("SV reads processed {}".format(count), err=True)
        count += 1

        n1 = (r.qname, r.flag)
        all_flags[r.qname][r.flag].append((r.rname, r.pos))

        # Limit to regions
        if regions:
            rname = infile.get_reference_name(r.rname)
            if rname not in regions or len(regions[rname].search(r.pos, r.pos+1)) == 0:
                continue

        scope.update(r)  # Add alignment to scope
        grey_added = 0

        for t, ol, sup_edge in scope.iterate():  # Other read, overlap current read, supplementary edge

            if not t:
                break

            n2 = (t.qname, t.flag)
            all_flags[t.qname][t.flag].append((t.rname, t.pos))

            # Four types of edges; black means definite match with overlapping soft-clips. Grey means similar
            # rearrangement with start and end co-ords on reference genome. Yellow means supplementary matches a primary
            # read; these edges need to be checked when both primary reads are available, change to black edge or delete
            # Finally white edge means a read1 to read2 edge
            if ol and not sup_edge:
                identity, prob_same = align_match(r, t)

                if prob_same > 0.01:  # Add a black edge
                    edges.append((n1, n2, {"p": [(r.rname, r.pos), (t.rname, t.pos)], "c": "b"}))

            elif not ol and not sup_edge:
                # Make sure both are discordant, mapped to same chromosome
                if (r.flag & 2) or (t.flag & 2) or r.rnext == -1 or r.rnext != t.rnext:
                    continue
                if abs(r.pnext - t.pnext) < max_dist:
                    # Add a grey edge
                    edges.append((n1, n2, {"p": [(r.rname, r.pos), (t.rname, t.pos)], "c": "g"}))
                    grey_added += 1
                    if grey_added > 60:
                        break  # Stop the graph getting unnecessarily dense

            elif sup_edge:
                # Add a yellow edge
                edges.append((n1, n2, {"p": [(r.rname, r.pos), (t.rname, t.pos)], "c": "y"}))

    # Add read-pair information to the graph, link regions together
    G = nx.Graph()
    G.add_edges_from(edges)
    new_edges = []
    for g in nx.connected_component_subgraphs(G):

        # Add white edges between read-pairs which are NOT in the subgraph
        nodes_to_check = [(n, all_flags[n].keys()) for n, f in g.nodes()]

        for n, flags in nodes_to_check:  # 2 reads, or 3 if supplementary read
            for f1, f2 in itertools.combinations_with_replacement(flags, 2):  # Check all edges within template
                u, v = (n, f1), (n, f2)
                has_v, has_u = g.has_node(v), g.has_node(u)

                # Only add an edge if either u or v is missing, not if both (or neither) missing
                if has_u != has_v:  # XOR
                   new_edges.append((u, v, {"p": list(set(all_flags[n][f1] + all_flags[n][f2])),  "c": "w"}))

    # Regroup based on white edges (link together read-pairs)
    G.add_edges_from(new_edges)

    if args["output"] == "-" or args["output"] is None:
        outfile = sys.stdout
    else:
        outfile = open(args["output"], "w")

    with outfile:
        head = ["chrA", "posA", "chrB", "posB", "svtype", "join_type", "nreads", "cipos95A", "cipos95B", "max_k_coverage",
         "Full_call", "soct_clip_stacks", "DPpri", "DApri", "DNpri", "NMpri", "SPpri", "ASpri", "EVsup", "DAsup",
         "NMsup", "ASsup", "best_sc", "contig"]

        outfile.write("\t".join(head) + "\n")
        c = 0

        for grp in nx.connected_component_subgraphs(G):  # Get large components, possibly multiple events
            reads, nodes_u, nodes_v, edge_data = get_reads(infile, grp, max_dist, rl=int(args['read_length']))

            events = assembler.merge_assemble(grp, reads, infile, args["clip_length"], args["insert_median"],
                                              args["insert_stdev"], args["read_length"])

            edges = grp.edges(data=True)
            edge_colors = {k: list(v) for k, v in itertools.groupby(edges, key=lambda x: x[2]['c'])}

            if len(events) == 0:
                # No assembly or linkage, look elsewhere
                echo("No assembly")
                continue

            for event_info, breakpoint_tuple in events:
                if len(breakpoint_tuple) == 1:
                    # Unlinked sub-cluster, partly assembled
                    echo("Unlinked")

                    echo(breakpoint_tuple)
                    echo(event_info)

                    echo("")
                elif len(breakpoint_tuple) == 2:
                    # Linked


                    # Make a call using the linked contigs as the base
                    call_info, contributing_reads = caller.call_break_points(breakpoint_tuple)
                    event_info.update(call_info)
                    # Get tag information for the supporting reads
                    event_info.update(score_reads(contributing_reads, reads))

                    outfile.write(caller.call_to_string(event_info))
                else:
                    raise ValueError("Too many sides to event")

                c += 1

        echo("Events called: {}".format(c))
