#!/bin/bash

if [ -z $1 ]; then
    echo "You must provide a filestem"
    exit 0;
fi

FILESTEM="$1"
shift #remove gene from the args

echo "Gene name is ${FILESTEM}"

#make sure bedtools is in the path
bed_path=$(which bamToBed)
 if [ ! -x "$bed_path" ]; then
    die "Can't find bamToBed, make sure bedTools is in your path."
 fi

BASE_DIR=/lustre/scratch109/sanger/`whoami`/gibson/mouse/${FILESTEM}_paired_crisprs
echo "Working dir is $BASE_DIR"

if [ ! -d "$BASE_DIR" ]; then
    die "$BASE_DIR doesn't exist, aborting"
fi

cd "$BASE_DIR"

#temporary fix to give us all the permutations in an easy to process format. replaces .yaml with .txt and outputs to that
perl -MYAML::Any=LoadFile -wE '(my $file = $ARGV[0]) =~ s/(.*)\.yaml/$1\.txt/; my $data = LoadFile($ARGV[0]) || die "yaml error"; open(my $fh, ">", $file) || die "Cant open $file"; while( my ($exon_id, $pairs) = each %{ $data } ) { say $fh "${exon_id}_" . $_->{left_crispr} . " ${exon_id}_" . $_->{right_crispr} for @{$pairs}  }' ${FILESTEM}_pairs.yaml

#need to get pairs too
echo "Would running windowbed on" `wc -l ${FILESTEM}_pairs.txt | awk '{print $1}'` "files"

#run windowbed on all the permutations
#mkdir html || die "couldn't make html file"
#while read a b; do
#    echo "Finding valid pairs for ${a} vs ${b}"
#    windowBed -a "${FILESTEM}.valid-${a}.bed" -b "${FILESTEM}.valid-${b}.bed" -w 99999999999 | awk '!($2==$8 && #$1==$7)' | perl ~ah19/work/paired_crisprs/valid_pairs.pl ${FILESTEM}_crisprs.fq > html/${a}_vs_${b}_valid.html
#done < permutations.txt

#change permutations to just be the data from pairs.yaml somehow. parse the yaml with perl maybe?
#are we gonna throw away ones with more than x off-targets?
mkdir paired_data || die "Couldn't make paired data folder"
WINDOWSIZE=9000
while read a b; do
    echo "Finding valid pairs for ${a} vs ${b}"
    #for each pair check all possible paired off-targets
    windowBed -a "${FILESTEM}.valid-${a}.bed" -b "${FILESTEM}.valid-${b}.bed" -w "$WINDOWSIZE" > "paired_data/${a}_vs_${b}.txt"
    windowBed -a "${FILESTEM}.valid-${a}.bed" -b "${FILESTEM}.valid-${a}.bed" -w "$WINDOWSIZE" >> "paired_data/${a}_vs_${b}.txt"
    windowBed -a "${FILESTEM}.valid-${b}.bed" -b "${FILESTEM}.valid-${b}.bed" -w "$WINDOWSIZE" >> "paired_data/${a}_vs_${b}.txt"
    #now write some kind of summary to parse each file and find the closest match?
done < ${FILESTEM}_pairs.txt

perl ~ah19/work/paired_crisprs/process_paired_data.pl --fq-file "${FILESTEM}_crisprs.fq" --pair-yaml-file "${FILESTEM}_pairs.yaml" --crispr-yaml-file "${FILESTEM}_crisprs.yaml" --paired-output paired_data/*

echo "Completed successfully."