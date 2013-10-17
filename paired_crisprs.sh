#!/bin/bash

#
# TODO:
#   stop using ah19 paths, make it configurable.
#   modify split bed to bail if a crispr has more than, say, 8k off targets.
#

function die () {
    if [ -n "$1" ]; then
        echo -e "$1"
    else
        echo "Die was called with no error message"
    fi
    exit 1
}

USAGE='Usage: paired_crisprs.sh <gene> <species> <exon ids>'

if [ -z "$LIMS2_DB" ]; then
    die "You must set the lims2 environment."
fi

if [ -z "$1" ]; then
    die "You must specify a gene:\n$USAGE"
fi

FILESTEM="$1"
shift #remove gene from the args

echo "Gene name is ${FILESTEM}"

if [ -z "$1" ]; then
    die "You must specify a species:\n$USAGE"
fi

SPECIES="$1"
shift

case "$SPECIES" in 
    [Hh]uman)
        GENOME=/lustre/scratch109/sanger/ah19/genomes/human/ncbi37/Homo_sapiens.NCBI36.48.dna.all.fa
        echo "Using human genome (${GENOME})"
        ;;
    [Mm]ouse)
        #i can no longer access the vrpipe one, so i've made a dir symlinking to srpipe
        #GENOME=/lustre/scratch105/vrpipe/refs/mouse/GRCm38/GRCm38_68.fa
        GENOME=/lustre/scratch109/sanger/ah19/genomes/mouse/GRCm38/Mus_musculus.GRCm38.68.dna.toplevel.fa
        echo "Using mouse genome (${GENOME})"
        ;;
    *)
        die "Species must be human or mouse:\n$USAGE"
        ;;
esac

if [ -z "$1" ]; then
    die "You must specify at least one exon id\n$USAGE"
fi

#make sure bedtools is in the path
bed_path=$(which bamToBed)
 if [ ! -x "$bed_path" ]; then
    die "Can't find bamToBed, make sure bedTools is in your path."
 fi

BASE_DIR=/lustre/scratch109/sanger/`whoami`/${FILESTEM}_paired_crisprs
echo "Working dir is $BASE_DIR"

if [ -d "$BASE_DIR" ]; then
    die "$BASE_DIR already exists, aborting"
fi

mkdir "$BASE_DIR" || die "Error creating $BASE_DIR"
cd "$BASE_DIR"

#write the genome and exons and stuff so we know what was run after the fact
echo -e "Gene name is ${FILESTEM}\nGenome used is ${GENOME}\nExons supplied:$@" > info.txt

#
# TOMORROW
#
# group off target output by mismatches so we can compare directly with zhang stuff. 
# remove --no-expand-seq
#
# important:
# do migration
# move code to farm3
# run on the exons we put in a comment in one of the files,
# insert crispr data
#
# 1: test the new find paired crisprs stuff. is the output correct? does it match hprt1 from before?
# 
# 3: how are we gonna put the off targets back in? i guess the same as before with ids
#

#
#what about UTR??? we need to label crisprs identified in the UTR. in the name somewhere so that later 
#it can be easily chosen
#
echo "Generating paired crisprs"
#perl ~ah19/work/paired_crisprs/find_paired_crisprs.pl "$@" | perl -MBio::Perl=revcom -we 'my $i = 1; my $current_exon; while( my $line = <> ) { next if $line =~ /^Exon/; chomp $line; my ($exon_id, $first, $spacer, $spacer_len, $second) = split ",", $line; if ( ! defined $current_exon || $current_exon ne $exon_id ) { $i = 1; $current_exon = $exon_id; } print "@" . $exon_id . "_" . $i . "A\n" . revcom($first)->seq . "\n"; print "@" . $exon_id . "_" . $i . "B\n" . $second . "\n"; $i++; }' > ${FILESTEM}_crisprs.fq || die "find_paired_crisprs.pl failed!"
perl ~ah19/work/paired_crisprs/find_paired_crisprs.pl --no-expand-seq --species "${SPECIES}" --exon-ids "$@" --fq-file "${FILESTEM}_crisprs.fq" --crispr-yaml-file "${FILESTEM}_crisprs.yaml" --pair-yaml-file "${FILESTEM}_pairs.yaml" || die "find_paired_crisprs.pl failed!"

echo 'Submitting bwa aln step'
bsub -K -o ${FILESTEM}_align.out -e ${FILESTEM}_align.err -G team87-grp -M 4000000 -R "select[mem>4000] rusage[mem=4000]" '/software/solexa/bin/bwa aln -n 6 -o 0 -l 21 -k 5 -N '"${GENOME}"' '"${FILESTEM}"'_crisprs.fq > '"${FILESTEM}"'.sai'

echo 'Submitting bwa samse step'
bsub -K -o ${FILESTEM}_samse.out -e ${FILESTEM}_samse.err -G team87-grp -M 4000000 -R "select[mem>4000] rusage[mem=4000]" '/software/solexa/bin/bwa samse -n 900000 '"${GENOME}"' '"${FILESTEM}"'.sai '"${FILESTEM}"'_crisprs.fq > '"${FILESTEM}"'.sam'

echo 'Getting sorted bam file'
/software/solexa/bin/aligners/bwa/current/xa2multi.pl $FILESTEM.sam | /software/solexa/bin/samtools view -bS - | /software/solexa/bin/samtools sort - ${FILESTEM}.sorted

echo 'Making bed file'
bamToBed -i ${FILESTEM}.sorted.bam > ${FILESTEM}.bed || die "bamToBed failed!"

echo 'Retrieving sequences'
bsub -K -o ${FILESTEM}_fasta.out -e ${FILESTEM}_fasta.err -G team87-grp -M 4000000 -R "select[mem>4000] rusage[mem=4000]" 'fastaFromBed -tab -fi '"${GENOME}"' -bed '"${FILESTEM}"'.bed -fo '"${FILESTEM}"'.with_seqs.tsv'

echo 'Merging sequences'
perl ~ah19/work/paired_crisprs/merge_fasta.pl ${FILESTEM}.bed ${FILESTEM}.with_seqs.tsv > ${FILESTEM}.with_seqs.bed || die "merge_fasta.pl failed!"

echo 'Removing invalid crisprs'
#perl ~ah19/work/paired_crisprs/remove_invalid_crisprs.pl ${FILESTEM}_crisprs.fq ${FILESTEM}.with_seqs.bed > "${FILESTEM}.valid.bed" || die "remove_invalid_crisprs.pl failed!"
perl ~ah19/work/paired_crisprs/remove_invalid_crisprs.pl --species "${SPECIES}" --fq-file "${FILESTEM}_crisprs.fq" --crispr-yaml-file "${FILESTEM}_crisprs.yaml" --bed-file "${FILESTEM}.with_seqs.bed" > "${FILESTEM}.valid.bed" || die "remove_invalid_crisprs.pl failed!"

echo 'Splitting bed file'
perl ~ah19/work/paired_crisprs/split_bed.pl "${FILESTEM}.valid.bed" yes > permutations.txt || die "split_bed.pl failed!"

#need to get pairs too
echo "Would running windowbed on" `wc -l permutations.txt | awk '{print $1}'` "files"

#run windowbed on all the permutations
#mkdir html || die "couldn't make html file"
#while read a b; do
#    echo "Finding valid pairs for ${a} vs ${b}"
#    windowBed -a "${FILESTEM}.valid-${a}.bed" -b "${FILESTEM}.valid-${b}.bed" -w 99999999999 | awk '!($2==$8 && #$1==$7)' | perl ~ah19/work/paired_crisprs/valid_pairs.pl ${FILESTEM}_crisprs.fq > html/${a}_vs_${b}_valid.html
#done < permutations.txt

#change permutations to just be the data from pairs.yaml somehow. parse the yaml with perl maybe?
#are we gonna throw away ones with more than x off-targets?
#mkdir paired_data || die "Couldn't make paired data folder"
#WINDOWSIZE=9000
#while read a b; do
#    echo "Finding valid pairs for ${a} vs ${b}"
#    #for each pair check all possible paired off-targets
#    windowBed -a "${FILESTEM}.valid-${a}.bed" -b "${FILESTEM}.valid-${b}.bed" -w "$WINDOWSIZE" > "paired_data/${a}_vs_${b}.txt"
#    windowBed -a "${FILESTEM}.valid-${a}.bed" -b "${FILESTEM}.valid-${a}.bed" -w "$WINDOWSIZE" >> "paired_data/${a}_vs_${b}.txt"
#    windowBed -a "${FILESTEM}.valid-${b}.bed" -b "${FILESTEM}.valid-${b}.bed" -w "$WINDOWSIZE" >> "paired_data/${a}_vs_${b}.txt"
    #now write some kind of summary to parse each file and find the closest match?
#done < permutations.txt

echo "Completed successfully."