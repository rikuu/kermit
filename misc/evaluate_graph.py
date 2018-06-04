import sys, pysam
from collections import defaultdict
from math import floor

if len(sys.argv) < 3:
    print('Usage: %s <alignments> <string graph>' % sys.argv[0])
    sys.exit(1)

ALIGNMENTS = sys.argv[1]
GRAPH = sys.argv[2]
INDEX_LEN = 5000
MIN_OVERLAP = 250

class Alignments:
    def __init__(self, filename):
        self.intervals = {}
        if filename[-3:] == 'sam' or filename[-3:] == 'bam':
            with pysam.AlignmentFile(filename, "rb") as f:
                for read in f:
                    if not read.is_unmapped and not read.query_name in self.intervals:
                        self.intervals[read.query_name] = (read.reference_name, read.reference_start, read.reference_end)
        else:
            with open(filename, 'r') as f:
                for line in f:
                    data = line.split()
                    self.intervals[data[0]] = (data[4], int(data[6]), int(data[7]))

        self.index = defaultdict(list)
        for read in self.intervals.keys():
            ref, start, end = self.intervals[read]
            si, ei = floor(start / INDEX_LEN), floor(end / INDEX_LEN) + 1
            while len(self.index[ref]) < ei: self.index[ref].append([])
            for i in range(si, ei):
                self.index[ref][i].append(read)

    def get_overlaps(self, v):
        if not v in self.intervals:
            return set()

        vr, v1, v2 = self.intervals[v]
        candidates = []
        for i in range(floor(v1 / INDEX_LEN), floor(v2 / INDEX_LEN) + 1):
            for w in self.index[vr][i]:
                if w != v and w not in candidates:
                    candidates.append(w)

        overlaps = []
        for read in candidates:
            wr, w1, w2 = self.intervals[read]
            overlap = min(v2, w2) - max(v1, w1)
            if wr == vr and overlap >= MIN_OVERLAP:
                overlaps.append(read)

        return set(overlaps)

    def distance(self, v, w):
        vr, v1, v2 = self.intervals[v]
        wr, w1, w2 = self.intervals[w]

        return (v1, v2, w1, w2, min(v2, w2) - max(v1, w1))

    def get_reads(self):
        return self.intervals.keys()

class Graph:
    def __init__(self, filename):
        self.overlaps = defaultdict(list)
        with open(filename, 'r') as f:
            for line in f:
                data = line.rstrip().split()
                assert data[0] == 'L'

                start, end = data[1].split(':')[0], data[3].split(':')[0]
                self.overlaps[start].append(end)
                self.overlaps[end].append(start)

    def get_overlaps(self, v):
        return set(self.overlaps[v])

alignments, graph = Alignments(ALIGNMENTS), Graph(GRAPH)
true_f, false_f, missing_f = open('true', 'w'), open('false', 'w'), open('missing', 'w')

all_alignment, all_graph = 0, 0
true_count, false_count, missing_count = 0, 0, 0
for read in alignments.get_reads():
    alignment_overlaps = alignments.get_overlaps(read)
    graph_overlaps = graph.get_overlaps(read)

    all_alignment += len(alignment_overlaps)
    all_graph += len(graph_overlaps)

    true = alignment_overlaps.intersection(graph_overlaps)
    false = graph_overlaps.difference(alignment_overlaps)
    missing = alignment_overlaps.difference(graph_overlaps)

    for edge in true:
        true_f.write(read + ' ' + edge + '\n')
    for edge in false:
        false_f.write(read + ' ' + edge + '\n')
    for edge in missing:
        missing_f.write(read + ' ' + edge + '\n')

    true_count += len(true)
    false_count += len(false)
    missing_count += len(missing)

print(all_alignment, all_graph)
print(true_count, missing_count, false_count)
print(100 * true_count / all_graph, 100 * false_count / all_graph)

true_f.close()
false_f.close()
missing_f.close()
