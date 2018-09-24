from collections import Counter
import itertools
from sklearn.cluster import KMeans
import numpy as np


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

