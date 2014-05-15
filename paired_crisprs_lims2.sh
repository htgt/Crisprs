#!/bin/bash

function die () {
    if [ -n "$1" ]; then
        echo -e "$1"
    else
        echo "Die was called with no error message"
    fi

    #echo "job failed (${FILESTEM} ${SPECIES}): $@" | mail -s "Lims2 paired crispr search failed" ah19@sanger.ac.uk
    echo "job failed (${FILESTEM} ${SPECIES}): $@" | mail -s "Lims2 paired crispr search failed" tg6@sanger.ac.uk

    exit 1
}

USAGE='Usage: paired_crisprs.sh <gene> <species> <exon ids>'

#SCRIPT_PATH=~ah19/work/paired_crisprs

if [ -z "$SCRIPT_PATH" ]; then
    die "You must set the script path"
fi

if [ -z "$LIMS2_REST_CLIENT_CONFIG" ]; then
    die "You must set LIMS2_REST_CLIENT_CONFIG!"
fi

if [ -z "$1" ]; then
    die "You must specify a gene:\n$USAGE"
fi

PATH=$SCRIPT_PATH:$PATH
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
        GENOME=/lustre/scratch109/sanger/ah19/genomes/human/GRCh37/Homo_sapiens.GRCh37.dna.all.fa
        ASSEMBLY=GRCh37
        echo "Using human genome (${GENOME})"
        ;;
    [Mm]ouse)
        #i can no longer access the vrpipe one, so i've made a dir symlinking to srpipe
        #GENOME=/lustre/scratch105/vrpipe/refs/mouse/GRCm38/GRCm38_68.fa
        GENOME=/lustre/scratch109/sanger/ah19/genomes/mouse/GRCm38/Mus_musculus.GRCm38.68.dna.toplevel.fa
        ASSEMBLY=GRCm38
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

BASE_DIR=/lustre/scratch109/sanger/`whoami`/gibson/${FILESTEM}_paired_crisprs
echo "Working dir is $BASE_DIR"

if [ -d "$BASE_DIR" ]; then
    die "$BASE_DIR already exists, aborting"
fi

mkdir -p "$BASE_DIR" || die "Error creating $BASE_DIR"
cd "$BASE_DIR"

#write the genome and exons and stuff so we know what was run after the fact
echo -e "Gene name is ${FILESTEM}\nGenome used is ${GENOME}\nExons supplied:$@" > info.txt

export LIMS2_REST_CLIENT_CONFIG=/nfs/team87/farm3_lims2_vms/conf/lims2-live-rest-client.conf

PATH=/software/perl-5.16.2/bin:$PATH

PERL5LIB=/nfs/users/nfs_a/ah19/lib/perl5/x86_64-linux-thread-multi:/nfs/team87/farm3_lims2_vms/software/perl/lib/perl5:/nfs/team87/farm3_lims2_vms/software/perl/lib/perl5/x86_64-linux-thread-multi:/nfs/users/nfs_a/ah19/work/WGE/lib:/software/perl-5.16.2/lib/5.16.2/x86_64-linux-thread-multi:/nfs/users/nfs_a/ah19/lib/perl5:/software/pubseq/PerlModules/BioPerl/1_6_920:$PERL5LIB

#remove any brave new world crap, temporary until i have a new bashrc
PERL5LIB=$(perl -we 'print join ":", grep { $_ !~ /brave_new_world|perl-5\.8\.9/ } split ":", $ENV{PERL5LIB};')
PATH=$(perl -we 'print join ":", grep { $_ !~ /brave_new_world|perl-5\.8\.9/ } split ":", $ENV{PATH};')

#
#what about UTR??? we need to label crisprs identified in the UTR. in the name somewhere so that later 
#it can be easily chosen
#
echo "Generating paired crisprs"
#perl ${SCRIPT_PATH}/find_paired_crisprs.pl "$@" | perl -MBio::Perl=revcom -we 'my $i = 1; my $current_exon; while( my $line = <> ) { next if $line =~ /^Exon/; chomp $line; my ($exon_id, $first, $spacer, $spacer_len, $second) = split ",", $line; if ( ! defined $current_exon || $current_exon ne $exon_id ) { $i = 1; $current_exon = $exon_id; } print "@" . $exon_id . "_" . $i . "A\n" . revcom($first)->seq . "\n"; print "@" . $exon_id . "_" . $i . "B\n" . $second . "\n"; $i++; }' > ${FILESTEM}_crisprs.fq || die "find_paired_crisprs.pl failed!"
perl ${SCRIPT_PATH}/find_paired_crisprs.pl --species "${SPECIES}" --regions "$@" --fq-file "crisprs.fq" --crispr-yaml-file "crisprs.yaml" --pair-yaml-file "pairs.yaml" || die "find_paired_crisprs.pl failed!"

echo 'Aligning and outputting to bed file'
/software/solexa/bin/bwa samse -n 100000000 "${GENOME}" <(/software/solexa/bin/bwa aln -n 6 -o 0 -l 21 -k 5 -N -m 1000000000 "${GENOME}" crisprs.fq) crisprs.fq | /software/solexa/pkg/bwa/current/xa2multi.pl | /software/hpag/samtools/0.1.19/bin/samtools view -bS - | /software/hpag/samtools/0.1.19/bin/samtools sort - "${FILESTEM}.sorted" || die "aligning failed!"

bamToBed -i "${FILESTEM}.sorted.bam" > "${FILESTEM}.bed" || die "Couldn't make bed file"

echo 'Retrieving sequences'
fastaFromBed -tab -fi "${GENOME}" -bed "${FILESTEM}".bed -fo "${FILESTEM}".with_seqs.tsv

echo 'Merging sequences'
perl ${SCRIPT_PATH}/merge_fasta.pl "${FILESTEM}.bed" "${FILESTEM}.with_seqs.tsv" > ${FILESTEM}.with_seqs.bed || die "merge_fasta.pl failed!"

#we don't need these now we have a merged one
rm "${FILESTEM}.bed" "${FILESTEM}.with_seqs.tsv"

echo 'Removing invalid crisprs'
#perl ${SCRIPT_PATH}/remove_invalid_crisprs.pl ${FILESTEM}_crisprs.fq ${FILESTEM}.with_seqs.bed > "${FILESTEM}.valid.bed" || die "remove_invalid_crisprs.pl failed!"
perl ${SCRIPT_PATH}/remove_invalid_crisprs.pl --species "${SPECIES}" --fq-file "crisprs.fq" --crispr-yaml-file "crisprs.yaml" --bed-file "${FILESTEM}.with_seqs.bed" > "${FILESTEM}.valid.bed" || die "remove_invalid_crisprs.pl failed!"

#only keep valid bed file
rm "${FILESTEM}.with_seqs.bed"


echo 'Splitting bed file'
perl ${SCRIPT_PATH}/split_bed.pl "${FILESTEM}.valid.bed" || die "split_bed.pl failed!"

#
# we now have all the crispr information, get the paired data
#

#temporary fix to give us all the permutations in an easy to process format. replaces .yaml with .txt and outputs to that
#this should be modified to be a real script that filters out any crisprs with lots of off targets
perl -MYAML::Any=LoadFile -wE '(my $file = $ARGV[0]) =~ s/(.*)\.yaml/$1\.txt/; my $data = LoadFile($ARGV[0]) || die "yaml error"; open(my $fh, ">", $file) || die "Cant open $file"; while( my ($exon_id, $pairs) = each %{ $data } ) { say $fh "${exon_id}_" . $_->{left_crispr} . " ${exon_id}_" . $_->{right_crispr} for @{$pairs}  }' pairs.yaml || die "Creating pairs.txt filed"

#run windowbed on all the possible combinations
echo "Running windowbed on" `wc -l pairs.txt | awk '{print $1}'` "files"

mkdir paired_data || die "Couldn't make paired data folder"
WINDOWSIZE=9000
while read a b; do
    echo "Finding valid pairs for ${a} vs ${b}"
    #for each pair check all possible paired off-targets
    windowBed -a "${FILESTEM}.valid-${a}.bed" -b "${FILESTEM}.valid-${b}.bed" -w "$WINDOWSIZE" | awk '!($2==$8 && $1==$7)' > "paired_data/${a}_vs_${b}.txt"
    windowBed -a "${FILESTEM}.valid-${a}.bed" -b "${FILESTEM}.valid-${a}.bed" -w "$WINDOWSIZE" | awk '!($2==$8 && $1==$7)' >> "paired_data/${a}_vs_${b}.txt"
    windowBed -a "${FILESTEM}.valid-${b}.bed" -b "${FILESTEM}.valid-${b}.bed" -w "$WINDOWSIZE" | awk '!($2==$8 && $1==$7)' >> "paired_data/${a}_vs_${b}.txt"
done < pairs.txt

echo "Parsing windowbed output"
perl ${SCRIPT_PATH}/process_paired_data.pl --fq-file "crisprs.fq" --pair-yaml-file "pairs.yaml" --crispr-yaml-file "crisprs.yaml" --paired-output paired_data/*

echo "Persisting data"
perl ${SCRIPT_PATH}/persist_crisprs.pl --species "${SPECIES}" --assembly "${ASSEMBLY}" --crispr-yaml "crisprs.yaml" --commit || die "Crispr persist step failed"

perl ${SCRIPT_PATH}/persist_pairs.pl --species "${SPECIES}" --crispr-yaml "crisprs.yaml" --pair-yaml "pairs.yaml" --commit || die "Pairs persist step failed"

#cd ..
#rm -rf "$BASE_DIR"

echo "Completed successfully."
