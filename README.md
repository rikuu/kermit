## Kermit
Kermit is a guided genome assembler using colored overlap graphs based on [miniasm].

## Getting Started
```sh
# Install kermit and minimap2 (requires zlib)
git clone --recursive https://github.com/rikuu/kermit && (cd kermit && make)
git clone https://github.com/lh3/minimap2 && (cd minimap2 && make)

# Color
minimap2/minimap2 -t8 -x map-pb reference.fa reads.fq > reference.paf
kermit/kermit-color reference.paf > reads.cf

# Overlap
minimap2/minimap2 -t8 -x ava-pb reads.fq reads.fq > reads.paf

# Layout
kermit/kermit -C reads.cf -f reads.fq reads.paf > reads.gfa
```

## Introduction
Kermit is heavily based on [miniasm] and as such shares most advantages and
disadvantages with it. [minimap2] is used to provide all-vs-all read
self-mappings to kermit. Kermit outputs an assembly graph in GFA format.

The main difference is an added coloring step in the pipeline. Kermit uses
minimap2 to map the reads to a reference genome (or draft assembly) and colors
reads based on the mappings.

Kermit supports assembly guiding using either a reference genome or a genetic
linkage map.

### Coloring
For reference-guided assembly, coloring can be done
```sh
kermit/kermit-color reference.paf > reads.cf
```

For linkage map-guided assembly, the markers 
```sh
kermit/kermit-color reference.paf map.txt > reads.cf
```

A linkage map is given as tab-separated files with each line defining a marker
by its position on a contig and a bin. [lepmap3] gives applicable linkage maps
and is recommended.

[miniasm]: https://github.com/lh3/miniasm
[minimap2]: https://github.com/lh3/minimap2
[lepmap3]: https://sourceforge.net/projects/lep-map3/
