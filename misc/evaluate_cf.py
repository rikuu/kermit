import sys, pysam
from collections import defaultdict
from math import floor

if len(sys.argv) != 4:
    print('Usage: %s <markers> <alignments> <colors>' % sys.argv[0])
    sys.exit(1)

MARKERS = sys.argv[1]
ALIGNMENTS = sys.argv[2]
COLORS = sys.argv[3]

# TODO: These two should be made faster
def successor(v, x):
    for i in v:
        if i[0] >= x:
            return i[1]
    return None

def predecessor(v, x):
    for i1, i2 in zip(v[:-1], v[1:]):
        if i1[0] == x or i2[0] > x:
            return i1[1]
        elif i2[0] == x:
            return i2[1]
    return None

markers = defaultdict(list)
with open(MARKERS, 'r') as f:
    for line in f:
        data = line.split('\t')
        markers[data[0]].append((int(data[1]), int(data[2])))

for contig in markers.keys():
    markers[contig] = sorted(markers[contig], key = lambda x: x[0])

limits = {}
with pysam.AlignmentFile(ALIGNMENTS, "rb") as f:
    for read in f:
        if not read.is_unmapped:
            trunc = read.reference_name.split('_l')[0]

            lower = predecessor(markers[trunc], read.reference_start)
            if lower == None:
                lower = successor(markers[trunc], read.reference_start)

            upper = successor(markers[trunc], read.reference_end)
            if upper == None:
                upper = predecessor(markers[trunc], read.reference_end)

            if lower != None and upper != None:
                if not read.query_name in limits:
                    limits[read.query_name] = (lower, upper)
                else:
                    plower, pupper = limits[read.query_name]
                    limits[read.query_name] = (min(lower, plower), max(upper, pupper))

coloured, unmapped = 0, 0
inside, outside = [], []
with open(COLORS, 'r') as f:
    for line in f:
        data = line.split('\t')
        read = data[0]

        coloured += 1
        if not read in limits:
            unmapped += 1
            continue

        # 16 low bits tell color within a contig
        cl, cu = int(data[1]) & ((1 << 16)-1), int(data[2]) & ((1 << 16)-1)
        lower, upper = limits[read]

        overlap = max(0, min(cu+1, upper+1) - max(cl, lower))
        inside.append(overlap)
        outside.append(max(0, cu - cl - overlap))

all_reads = len(limits.keys())
print(all_reads, coloured, unmapped)

print(min(inside), max(inside), sum(inside) / len(inside))
print(min(outside), max(outside), sum(outside) / len(outside))

reads_inside = sum([1 for i in range(len(inside)) if inside[i] > 0 and outside[i] == 0])
reads_outside = sum([1 for i in range(len(inside)) if inside[i] == 0 and outside[i] > 0])
print(reads_inside, reads_outside)
print(100 * reads_inside / (coloured-unmapped), 100 * reads_outside / (coloured-unmapped))
