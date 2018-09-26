from collections import Counter
import itertools
from sklearn.cluster import KMeans
import numpy as np
import click
from src import assembler


def echo(*args):
    click.echo(args, err=True)


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


def call_break_points(break_points, thresh=500):
    """
    Makes a call from a list of break points. Can take a list of break points, or one merged cluster of breaks.
    Outliers are dropped. Breakpoints are clustered using kmeans into sets
    :param break_points: A 4 tuple (3 or 5 join, set([(read_name, flag)..]), chromosome, break point position)
    :param thresh: the distance threshold to determine if clustered
    :return: Info dict containing a summary of the call
    """
    break_points = sorted(break_points, key=lambda x: (x[2], x[3]))  # By chromosome and co-ordinate

    # If number chroms > 2, reduce
    chroms = Counter([i[2] for i in break_points])
    if len(chroms) > 2:
        c = [i[0] for i in sorted(chroms.items(), key=lambda x: x[1], reverse=True)][:2]
        chroms = {k: v for k, v in chroms.items() if k in c}
        break_points = [i for i in break_points if i[2] in chroms]

    # If minimum chrom counts < 0.1, reduce. Drop low coverage chromosomes
    if len(chroms) == 2:
        ci = list(chroms.items())
        total = float(sum(chroms.values()))
        c = None
        if ci[0][1] / total < 0.1:
            c = ci[1][0]  # Keep item 1
        elif ci[1][1] / total < 0.1:
            c = ci[0][0]  # Keep item 0
        if c:
            break_points = [i for i in break_points if i[2] == c]

    def separate_by_kmeans(break_points):
        # Cluster breakpoints by join-side and position
        X = [[i[0], np.log10(i[3])] for i in break_points]  # Log10 to scale the genome position down
        k_means = KMeans(init='k-means++', n_clusters=2, n_init=5, max_iter=20)
        labels = k_means.fit_predict(X)
        g = list(itertools.groupby(sorted(zip(labels, break_points)), key=lambda x: x[0]))

        c1 = [j[1] for j in g[0][1]]
        if len(g) > 1:
            c2 = [j[1] for j in g[1][1]]
        else:
            c2 = []
        if len(c2) > len(c1):
            c1p = c1
            c1 = c2
            c2 = c1p
        return c1, c2

    c1, c2 = [], []
    if len(break_points) == 1:
        c1 = break_points

    elif len(break_points) == 2:
        c1, c2 = [break_points[0]], [break_points[1]]

    # If more than 2 break points, drop breaks outside of template range that are on the same chromosome
    # No clear if this is necessary
    elif len(break_points) > 2:
        # Split into genomic clusters, keep biggest 2 only
        clst = []
        current = []
        for i in break_points:
            if len(current) == 0:
                current.append(i)
                continue
            chrom, pos = i[2], i[3]
            last_chrom, last_pos = current[-1][2], current[-1][3]
            if chrom == last_chrom and abs(last_pos - pos) < thresh:
                current.append(i)
            else:
                clst.append(current)
                current = [i]
        if len(current) > 0:
            clst.append(current)
        if len(clst) == 2:
            c1, c2 = clst
        if len(clst) > 2:  # Choose largest 2 clusters
            c1, c2 = [], []
            for item in clst:
                if len(item) > len(c1):
                    c1 = item
                elif len(item) > len(c2):
                    c2 = item
        elif len(clst) == 1:  # Use kmeans to seperate based on strand and position
            c1, c2 = separate_by_kmeans(break_points)

    info = {}
    count = 0
    contributing_reads = set([])
    for grp in [c1, c2]:
        if len(grp) == 0:
            continue  # When no c2 is found
        for i in grp:
            contributing_reads = contributing_reads.union(i[1])

        chrom = Counter([i[2] for i in grp]).most_common()[0][0]
        grp = [i for i in grp if i[2] == chrom]
        sc_side = Counter([i[0] for i in grp]).most_common()[0][0]
        bp = [i[3] for i in grp]
        mean_pos = int(np.mean(bp))
        mean_95 = abs(int(np.percentile(bp, [97.5])) - mean_pos)
        if count == 0:
            side = "A"
        else:
            side = "B"
        info["chr" + side] = chrom
        info["pos" + side] = mean_pos
        info["cipos95" + side] = mean_95
        info["join" + side] = sc_side
        count += 1
    info["nreads"] = len(contributing_reads)

    if "chrB" in info and info["chrA"] == info["chrB"]:
        if info["joinA"] == info["joinB"]:
            info["svtype"] = "INV"
            info["join_type"] = str(info["joinA"]) + "to" + str(info["joinB"])
        else:
            if info["posA"] <= info["posB"]:
                x = "joinA"
                y = "joinB"
            else:
                x = "joinB"
                y = "joinA"
            if info[x] == 3 and info[y] == 5:
                info["svtype"] = "DEL"
                info["join_type"] = "3to5"
            elif info[x] == 5 and info[y] == 3:
                info["svtype"] = "DUP"
                info["join_type"] = "5to3"

    elif "chrB" in info:
        info["svtype"] = "TRA"
        info["join_type"] = str(info["joinA"]) + "to" + str(info["joinB"])
    else:
        info["svtype"] = "INS"
        if "joinA" in info:
            info["join_type"] = str(info["joinA"]) + "to?"
        else:
            info["join_type"] = "?to?"

    return info, contributing_reads




def call_to_string(call_info):
    # Tab delimited string, in a bedpe style
    k = ["chrA", "posA", "chrB", "posB", "svtype", "join_type", "nreads", "cipos95A", "cipos95B", "max_k_coverage",
         "Full_call", "soct_clip_stacks", "DPpri", "DApri", "DNpri", "NMpri", "SPpri", "ASpri", "EVsup", "DAsup",
         "NMsup", "ASsup", "best_sc", "contig"]

    return "\t".join([str(call_info[ky]) if ky in call_info else "NA" for ky in k]) + "\n"



def single(bm, parent_graph, reads, bam, assemblies):

    # Get edges from parent graph
    echo(assemblies[list(bm.nodes())[0]])

    echo(parent_graph.nodes(data=True))
    cluster_pos = []
    for n in list(bm.nodes())[0]:
        echo(n)
        echo(parent_graph.node[n])
    echo("----")


def one_edge(bm, parent_graph, reads, bam, assemblies):

    for node_set in bm.nodes():
        echo(assemblies[node_set])

    pass


def multi(bm, parent_graph, reads, bam, assemblies):

    # Score edges
    edge_scores = []
    for u, v, data in bm.edges(data=True):
        score = data["weight"]
        if "result" in data:
            score = sum(data["result"].values())
        edge_scores.append((u, v, score))

    edge_scores = sorted(edge_scores, key=lambda k: k[2], reverse=True)
    seen = set([])
    for u, v, scr in edge_scores:
        if u in seen or v in seen:
            continue
        seen.add(u)
        seen.add(v)
        sub = bm.subgraph([u, v])
        yield one_edge(sub, parent_graph, reads, bam, assemblies)


def call_from_block_model(bm, parent_graph, reads, bam):

    # Try and assemble nodes in the block model
    assemblies = {}
    for node_set in bm.nodes():
        sub = parent_graph.subgraph(node_set)
        if any(i[2]["c"] == "b" for i in sub.edges(data=True)):
            # assemble reads if any overlapping soft-clips
            assemblies[node_set] = assembler.base_assemble(sub, reads, bam)
        else:
            assemblies[node_set] = {}

    if len(bm.edges) > 1:
        # Score edges
        echo('multi edges')

        for p in multi(bm, parent_graph, reads, bam, assemblies):
            echo("yielded", p)
        echo("----")

    if len(bm.nodes()) == 1:

        # Single isolated node
        echo('single node')
        return
        single(bm, parent_graph, reads, bam, assemblies)


    if len(bm.edges) == 1:

        # Easy case
        echo('one edge')
        return
        one_edge(bm, parent_graph, reads, bam, assemblies)
