## Kermit

Kermit is a guided genome assembler using colored string graphs.

## Getting Started

```sh
# Color
minimap/minimap -t8 reference.fa reads.fq | gzip -1 > reference.paf.gz
kermit/kermit-color reference.paf.gz | gzip -1 > reads.cf.gz
# Overlap
minimap/minimap -x ava10k -t8 reads.fq reads.fq | gzip -1 > reads.paf.gz
# Layout
kermit/kermit -f reads.fq -C reads.cf.gz reads.paf.gz > reads.gfa
```

[unitig]: http://wgs-assembler.sourceforge.net/wiki/index.php/Celera_Assembler_Terminology
[minimap]: https://github.com/lh3/minimap
[miniasm]: https://github.com/lh3/miniasm
[paf]: https://github.com/lh3/miniasm/blob/master/PAF.md
[gfa]: https://github.com/pmelsted/GFA-spec/blob/master/GFA-spec.md
[sg]: http://bioinformatics.oxfordjournals.org/content/21/suppl_2/ii79.abstract
