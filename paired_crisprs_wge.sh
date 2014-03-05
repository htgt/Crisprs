#!/bin/bash

#
#   the intention is to bsub this job, it will do no bsubbing 
#

function die () {
    if [ -n "$1" ]; then
        echo -e "$1"
    else
        echo "Die was called with no error message"
    fi

    #email me if it fails
    echo "job failed (${FILESTEM} ${SPECIES}): $@" | mail -s "Paired crispr search failed" ah19@sanger.ac.uk
    
    #if we have died update this pair accordingly
    update_crispr_pair.pl --species "${SPECIES}" --pair-ids "$FILESTEM" --status "-1"

    exit 1
}

USAGE='Usage: paired_crisprs.sh <gene> <species> <exon ids>'

SCRIPT_PATH=~ah19/work/paired_crisprs

if [ -z "$LIMS2_REST_CLIENT_CONFIG" ]; then
    die "You must set LIMS2_REST_CLIENT_CONFIG!"
fi

# if [ -z "$LIMS2_DB" ]; then
#     die "You must set the lims2 environment."
# fi

if [ -z "$1" ]; then
    die "You must specify a gene:\n$USAGE"
fi

#NOTE: for WGE this MUST be the pair id
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
        GENOME=/lustre/scratch110/sanger/ah19/genomes/human/GRCh37/Homo_sapiens.GRCh37.dna.all.fa
        ASSEMBLY=GRCh37
        echo "Using human genome (${GENOME})"
        ;;
    [Mm]ouse)
        #i can no longer access the vrpipe one, so i've made a dir symlinking to srpipe
        #GENOME=/lustre/scratch105/vrpipe/refs/mouse/GRCm38/GRCm38_68.fa
        GENOME=/lustre/scratch110/sanger/ah19/genomes/mouse/GRCm38/Mus_musculus.GRCm38.68.dna.toplevel.fa
        ASSEMBLY=GRCm38
        echo "Using mouse genome (${GENOME})"
        ;;
    *)
        die "Species must be human or mouse:\n$USAGE"
        ;;
esac

#we now require a pair id in the format LeftCrisprID_RightCrisprID
#so we can just update a single crispr if necessary, but still update pair data

if [ -z "$1" ]; then
    die "You must specify at least one crispr id\n$USAGE"
fi

#make sure bedtools is in the path
bed_path=$(which bamToBed)
 if [ ! -x "$bed_path" ]; then
    die "Can't find bamToBed, make sure bedTools is in your path."
 fi

#BASE_DIR=/lustre/scratch110/sanger/`whoami`/wge/${FILESTEM}_paired_crisprs
BASE_DIR=/lustre/scratch109/sanger/team87/crisprs_wge/${FILESTEM}_paired_crisprs
echo "Working dir is $BASE_DIR"

if [ -d "$BASE_DIR" ]; then
    die "$BASE_DIR already exists, aborting"
fi

mkdir -p "$BASE_DIR" || die "Error creating $BASE_DIR"
cd "$BASE_DIR"

#write the genome and exons and stuff so we know what was run after the fact
echo -e "Gene name is ${FILESTEM}\nGenome used is ${GENOME}\nExons supplied:$@" > info.txt

#PATH=/software/perl-5.16.2/bin:/nfs/users/nfs_a/ah19/work/WGE/bin:$PATH

#PERL5LIB=/nfs/team87/farm3_lims2_vms/software/perl/lib/perl5:/nfs/team87/farm3_lims2_vms/software/perl/lib/perl5/x86_64-linux-thread-multi:/nfs/users/nfs_a/ah19/work/WGE/lib:/software/perl-5.16.2/lib/5.16.2/x86_64-linux-thread-multi:/nfs/users/nfs_a/ah19/lib/perl5:/software/pubseq/PerlModules/BioPerl/1_6_920:/software/pubseq/PerlModules/Ensembl/www_73_1/ensembl/modules:$PERL5LIB

#remove any brave new world crap, temporary until i have a new bashrc
PERL5LIB=$(perl -we 'print join ":", grep { $_ !~ /brave_new_world|perl-5\.8\.9/ } split ":", $ENV{PERL5LIB};')
PATH=$(perl -we 'print join ":", grep { $_ !~ /brave_new_world|perl-5\.8\.9/ } split ":", $ENV{PATH};')

#export WGE_REST_CLIENT_CONFIG=/nfs/team87/farm3_lims2_vms/conf/wge-devel-rest-client.conf

echo "Generating paired crisprs"
# perl /nfs/users/nfs_a/ah19/work/WGE/bin/find_paired_crisprs.pl --environment development --species "${SPECIES}" --gene-ids "$@" --fq-file "crisprs.fq" --crispr-yaml-file "crisprs.yaml" --pair-yaml-file "pairs.yaml" || die "find_paired_crisprs.pl failed!"

#change back to this one once we know remove_invalid_crisprs is as expected
find_crisprs.pl --fq-file "crisprs.fq" --crispr-yaml-file "crisprs.yaml" --species "${SPECIES}" --ids "$@" || die "retrieving crisprs from db failed!"

update_crispr_pair.pl --species "${SPECIES}" --pair-ids "$FILESTEM" --status 2 || die "Couldn't update pair in WGE"

echo 'Aligning and outputting to bed file'
/software/solexa/bin/bwa samse -n 100000000 "${GENOME}" <(/software/solexa/bin/bwa aln -n 5 -o 0 -l 20 -k 4 -N -m 1000000000 "${GENOME}" crisprs.fq) crisprs.fq | /software/solexa/pkg/bwa/current/xa2multi.pl | /software/hpag/samtools/0.1.19/bin/samtools view -bS - | /software/hpag/samtools/0.1.19/bin/samtools sort - ${FILESTEM}.sorted || die "aligning failed!"

bamToBed -i "${FILESTEM}.sorted.bam" > "${FILESTEM}.bed" || die "Couldn't make bed file"

echo 'Retrieving sequences'
fastaFromBed -tab -fi "${GENOME}" -bed "${FILESTEM}".bed -fo "${FILESTEM}".with_seqs.tsv

echo 'Merging sequences'
perl ${SCRIPT_PATH}/merge_fasta.pl "${FILESTEM}.bed" "${FILESTEM}.with_seqs.tsv" > ${FILESTEM}.with_seqs.bed || die "merge_fasta.pl failed!"

#we don't need these now we have a merged one
rm "${FILESTEM}.bed" "${FILESTEM}.with_seqs.tsv"

echo 'Removing invalid crisprs'
perl ${SCRIPT_PATH}/remove_invalid_crisprs.pl --species "${SPECIES}" --fq-file "crisprs.fq" --crispr-yaml-file "crisprs.yaml" --bed-file "${FILESTEM}.with_seqs.bed" --max-edit-distance 5 > "${FILESTEM}.valid.bed" || die "remove_invalid_crisprs.pl failed!"

#only keep valid bed file
rm "${FILESTEM}.with_seqs.bed"

#perl ${SCRIPT_PATH}/tmp/wge_test.pl --species ${SPECIES} --bed-file "${FILESTEM}.valid.bed" --crispr-yaml-file crisprs.yaml --max-offs 5000

echo "Persisting crisprs"
update_crispr_pair.pl --species "${SPECIES}" --pair-ids "$FILESTEM" --status 3 || die "Couldn't update pair in WGE"

#this needs to be +x and in the path
persist_crisprs.pl --species "${SPECIES}" --crispr-yaml-file "crisprs.yaml" --bed-file "${FILESTEM}.valid.bed" --max-offs 2000 --commit || die "Crispr persist step failed!"

#this assumes we just got a single valid pair.
echo "Calculating paired off targets"
update_crispr_pair.pl --species "${SPECIES}" --pair-ids "$FILESTEM" --update-offs || die "Couldn't update pair in WGE"

echo "Completed successfully."

rm -rf $BASE_DIR
