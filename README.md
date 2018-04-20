## Kermit

Kermit is a guided genome assembler using coloured string graphs.

## Getting Started

```sh
# Color
minimap2/minimap2 -x map-pb reference.fa reads.fq > reference.paf
kermit/kermit-color reference.paf [map1.txt map2.txt ...] > reads.cf
# Overlap
minimap2/minimap2 -x ava-pb reads.fq reads.fq > reads.paf
# Layout
kermit/kermit -C reads.cf -f reads.fq reads.paf > reads.gfa
```

[minimap]: https://github.com/lh3/minimap
[miniasm]: https://github.com/lh3/miniasm
[paf]: https://github.com/lh3/miniasm/blob/master/PAF.md
[gfa]: https://github.com/pmelsted/GFA-spec/blob/master/GFA-spec.md
[sg]: http://bioinformatics.oxfordjournals.org/content/21/suppl_2/ii79.abstract
