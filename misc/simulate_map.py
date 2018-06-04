import sys
from math import ceil
from random import randint

if len(sys.argv) != 4:
    print('Usage: %s <reference> <markers> <max length>' % sys.argv[0])
    sys.exit(1)

REFERENCE = sys.argv[1]
MARKERS = int(sys.argv[2])
MAX_LENGTH = int(sys.argv[3])

# Read reference genome to dict object
reference = {}
with open(REFERENCE, 'r') as f:
    name, seq = '', ''
    for line in f:
        if line[0] == '>':
            if seq != '':
                reference[name] = len(seq)
            name = line[1:].split()[0] #.replace(' ', '_')
            seq = ''
        else:
            seq += line[:-1]
    reference[name] = len(seq)

full_length = sum(reference.values())
contig_counts = {contig: ceil(MARKERS * (reference[contig] / full_length)) for contig in reference.keys()}

marker_positions = {contig: sorted([randint(0, reference[contig]) for i in range(contig_counts[contig])]) for contig in reference.keys()}

bins, distances, marker_count = 0, [], 0
for contig, positions in marker_positions.items():
    for current, next_position in zip(positions[:-1], positions[1:]):
        distance = next_position - current
        if distance > MAX_LENGTH:
            bins += 1

        distances.append(distance)
        print('%s\t%i\t%i' % (contig, current, bins))
        marker_count += 1
    print('%s\t%i\t%i' % (contig, positions[-1], bins))
    marker_count += 1
    bins += 1

# TODO: Generate fragmented (mis-)assembly

contig_count = len(reference.keys())
print('[simulate_map.py] contigs: %i, bases: %i' % (contig_count, full_length), file=sys.stderr)
print('[simulate_map.py] markers: %i, bins: %i' % (marker_count, bins), file=sys.stderr)
print('[simulate_map.py] base densities, marker: %.3f, bin: %.3f' % (marker_count / full_length, bins / full_length), file=sys.stderr)
print('[simulate_map.py] contig densities, marker: %.3f, bin: %.3f' % (marker_count / contig_count, bins / contig_count), file=sys.stderr)
print('[simulate_map.py] marker distances: (%i, %i, %.3f)' % (min(distances), max(distances), sum(distances) / len(distances)), file=sys.stderr)
