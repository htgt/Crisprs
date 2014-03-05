#!/bin/bash

function die () {
    if [ -n "$1" ]; then
        echo -e "$1"
    else
        echo "Die was called with no error message"
    fi

    #email me if it fails
    echo "job failed (${BASE_DIR} ${SPECIES}): $@" | mail -s "Paired crispr search failed" ah19@sanger.ac.uk

    exit 1
}

USAGE='Usage: paired_crisprs.sh <folder> <species> <exon ids>'

SCRIPT_PATH=~ah19/work/paired_crisprs

if [ -z "$1" ]; then
    die "You must specify a folder:\n$USAGE"
fi

#NOTE: for WGE this MUST be the pair id
BASE_DIR="$1"
shift #remove dir from the args

echo "Working dir is ${BASE_DIR}"

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

#make sure bedtools is in the path
bed_path=$(which bamToBed)
if [ ! -x "$bed_path" ]; then
   die "Can't find bamToBed, make sure bedTools is in your path."
fi

if [ ! -d "$BASE_DIR" ]; then
    mkdir -p $BASE_DIR || die "Error creating $BASE_DIR"
fi

cd "$BASE_DIR"

#write the genome and exons and stuff so we know what was run after the fact
echo -e "Genome used is ${GENOME}\nExons supplied:$@" > info.txt

PATH=/software/perl-5.16.2/bin:/nfs/users/nfs_a/ah19/work/WGE/bin:$PATH

PERL5LIB=/nfs/team87/farm3_lims2_vms/software/perl/lib/perl5:/nfs/team87/farm3_lims2_vms/software/perl/lib/perl5/x86_64-linux-thread-multi:/nfs/users/nfs_a/ah19/work/WGE/lib:/software/perl-5.16.2/lib/5.16.2/x86_64-linux-thread-multi:/nfs/users/nfs_a/ah19/lib/perl5:/software/pubseq/PerlModules/BioPerl/1_6_920:/software/pubseq/PerlModules/Ensembl/www_73_1/ensembl/modules:$PERL5LIB

#remove any brave new world crap, temporary until i have a new bashrc
PERL5LIB=$(perl -we 'print join ":", grep { $_ !~ /brave_new_world|perl-5\.8\.9/ } split ":", $ENV{PERL5LIB};')
PATH=$(perl -we 'print join ":", grep { $_ !~ /brave_new_world|perl-5\.8\.9/ } split ":", $ENV{PATH};')

export WGE_REST_CLIENT_CONFIG=/nfs/team87/farm3_lims2_vms/conf/wge-devel-rest-client.conf

echo "Generating paired crisprs"
find_crisprs.pl --fq-file "crisprs.fq" --crispr-yaml-file "crisprs.yaml" --species "${SPECIES}" --ids "$@" || die "retrieving crisprs from db failed!"

echo 'Aligning and outputting to bed file'
/software/solexa/bin/bwa samse -n 100000000 "${GENOME}" <(/software/solexa/bin/bwa aln -n 5 -o 0 -l 20 -k 4 -N -m 1000000000 "${GENOME}" crisprs.fq) crisprs.fq | /software/solexa/pkg/bwa/current/xa2multi.pl | /software/hpag/samtools/0.1.19/bin/samtools view -bS - | /software/hpag/samtools/0.1.19/bin/samtools sort - ots.sorted || die "aligning failed!"

bamToBed -i "ots.sorted.bam" > "all_offs.bed" || die "Couldn't make bed file"

echo 'Retrieving sequences'
fastaFromBed -tab -fi "${GENOME}" -bed "all_offs.bed" -fo "all_offs.with_seqs.tsv"

echo 'Merging sequences'
perl ${SCRIPT_PATH}/merge_fasta.pl "all_offs.bed" "all_offs.with_seqs.tsv" > "all_offs.with_seqs.bed" || die "merge_fasta.pl failed!"

#we don't need these now we have a merged one
rm "all_offs.bed" "all_offs.with_seqs.tsv"

echo 'Removing invalid crisprs'
perl ${SCRIPT_PATH}/remove_invalid_crisprs.pl --species "${SPECIES}" --fq-file "crisprs.fq" --crispr-yaml-file "crisprs.yaml" --bed-file "all_offs.with_seqs.bed" --max-edit-distance 5 > "ots.valid.bed" || die "remove_invalid_crisprs.pl failed!"

#only keep valid bed file
rm "all_offs.with_seqs.bed"

echo "Persisting crisprs"
persist_crisprs.pl --species "${SPECIES}" --crispr-yaml-file "crisprs.yaml" --bed-file "ots.valid.bed" --max-offs 2000 --commit || die "Crispr persist step failed!"

echo "Completed successfully."
