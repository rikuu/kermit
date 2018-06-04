#!/bin/bash
shopt -s extglob
set -eu
set -o pipefail

DIR=$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)
DATA=$DIR/elegans/

# Tools
TIME=/usr/bin/time
SIMLORD=~/simlord/simlord
MINIMAP=~/minimap2/minimap2
MINIASM=~/miniasm/miniasm
KERMIT=~/kermit/kermit
KERMITCOLOR=~/kermit/kermit-color
RACON=~/racon/build/bin/racon

# Simulation parameters
REFERENCE=$DATA/c_elegans.PRJNA275000.WS259.genomic.fa
READS=$DATA/reads.fastq
MAP=$DATA/map.txt
MARKERS=750000
BIN_DIST=200
COVERAGE=40

# Simulate reads
$SIMLORD --read-reference $REFERENCE --coverage $COVERAGE reads
ALN=reads.fastq.sam

# Simulate map
python3 $DIR/simulate_map.py $REFERENCE $MARKERS $BIN_DIST > $MAP

# Map reads
$TIME -v $MINIMAP -t16 -x ava-pb $READS $READS > ava.paf 2> ava.stderr
$TIME -v $MINIMAP -t16 -x map-pb $REFERENCE $READS > toref.paf 2> toref.stderr

# Miniasm
$MINIASM -S5 -p sg -f $READS ava.paf > miniasm.5.sg
$MINIASM -p sg -f $READS ava.paf > miniasm.sg
$TIME -v $MINIASM -f $READS ava.paf > miniasm.gfa 2> miniasm.stderr

# Kermit
$TIME -v $KERMITCOLOR toref.paf $MAP > colors.cf 2> color.stderr

$KERMIT -S5 -p sg -C colors.cf -f $READS ava.paf > kermit.5.sg
$KERMIT -p sg -C colors.cf -f $READS ava.paf > kermit.sg
$TIME -v $KERMIT -C colors.cf -f $READS ava.paf > kermit.gfa 2> kermit.stderr

# Evaluate coloring and cleaning
python3 $DIR/evaluate_cf.py $MAP $ALN colors.cf > eval-color.txt
python3 $DIR/evaluate_cf.py $MAP $ALN propagated.cf > eval-color-propagated.txt

python3 $DIR/evaluate_graph.py $ALN miniasm.5.sg > eval-graph-miniasm.5.txt
python3 $DIR/evaluate_graph.py $ALN kermit.5.sg > eval-graph-kermit.5.txt

python3 $DIR/evaluate_graph.py $ALN miniasm.sg > eval-graph-miniasm.txt
python3 $DIR/evaluate_graph.py $ALN kermit.sg > eval-graph-kermit.txt

# Evaluate assembly
awk '/^S/{print ">"$2"\n"$3}' kermit.gfa > kermit.fa
$TIME -v $MINIMAP -x map-pb -t16 kermit.fa $READS > kermit-racon.paf 2> kermit-racon-map.stderr
$TIME -v $RACON -t 16 $READS kermit-racon.paf kermit.fa > kermit-consensus.fa 2> kermit-racon.stderr

awk '/^S/{print ">"$2"\n"$3}' miniasm.gfa > miniasm.fa
$TIME -v $MINIMAP -x map-pb -t16 miniasm.fa $READS > miniasm-racon.paf 2> miniasm-racon-map.stderr
$TIME -v $RACON -t 16 $READS miniasm-racon.paf miniasm.fa > miniasm-consensus.fa 2> miniasm-racon.stderr

quast.py --threads 16 -R $REFERENCE miniasm-consensus.fa kermit-consensus.fa