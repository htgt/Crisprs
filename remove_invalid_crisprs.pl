#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use LIMS2::Util::EnsEMBL;
use Bio::Perl;
use Try::Tiny;
use Data::Dumper;
use JSON;
use YAML::Any qw( DumpFile LoadFile );

my ( $species, $crispr_fq_file, $crispr_yaml_file, $bed_file );
my $MAX_EDIT_DISTANCE = 6;

GetOptions(
    "help"                => sub { pod2usage( 1 ) },
    "man"                 => sub { pod2usage( 2 ) },
    "species=s"           => \$species,
    "fq-file=s"           => \$crispr_fq_file,
    "crispr-yaml-file=s"  => \$crispr_yaml_file,
    "bed-file=s"          => \$bed_file,
    'max-edit-distance=i' => \$MAX_EDIT_DISTANCE,
) or pod2usage( 2 );

#basically everything is required
die pod2usage( 2 ) unless $species and $crispr_yaml_file and $crispr_fq_file and $bed_file;

my $e = LIMS2::Util::EnsEMBL->new( species => $species );

#build a hash of all the crisprs so we can check hamming distances to be sure of orientation
my %crisprs;
parse_fq_file( \%crisprs );

#we store all the hamming distances so we can create a summary later.
#we also count the number of lines that have N sequence so we can say how many were thrown away
my ( $num_n_lines, %hd_data );

open(my $bed_fh, "<", $bed_file) or die "Couldn't open $bed_file: $!";
while ( my $line = <$bed_fh> ) {
    chomp $line;
    my ( $chr, $start, $end, $name, $unknown, $strand ) = split /\s+/, $line;

    my $seq;
    if ( $name =~ /(.+)-([CTAGNctagn]{23})/ ) {
        #print STDERR "Split $name into $1 $2\n";
        $seq = $2;
        $name = $1;
    }
    else {
        print STDERR "$name has no seq in it, fetching from ensembl\n";
        $seq = fetch_ensembl_seq( $chr, $start, $end );
    }

    die "No pam sites available for $name" unless exists $crisprs{$name};

    #skip all N sequences. we get a lot because bwa replace Ns in genome with random sequences
    if ( $seq =~ /N+/ ) {
        #push @n_lines, "$line\t$seq\n";
        $num_n_lines++;
        next;
    }

    #forward and reverse hamming distances so we don't have to keep computing them
    my $rev_hd = hamming_distance( $seq, $crisprs{$name}->{rev} ); 
    my $fwd_hd = hamming_distance( $seq, $crisprs{$name}->{fwd} );

    #names are in the format ENSE00003489858_7B. make sure the name matches or everything would break
    my ( $exon_id, $crispr_id ) = split '_', $name;
    die "fq header $name is not in the format EXONID_CRISPRID" unless $exon_id and $crispr_id;

    #if we havent looked at this crispr yet create an off targets hash with all the values set to 0,
    #so that when we take it out of the db no fields are missing.
    unless ( defined $hd_data{$exon_id}->{$crispr_id} ) {
        $hd_data{$exon_id}->{$crispr_id}{off_targets} = { map { $_ => 0 } (0 .. 5) };
    }

    #we're assuming reverse complemented pam, so we know orientation
    #check hamming distance to make sure the orientation is as we expect

    #
    # mismatch in the pam gets counted here. how can we get rid of it? 
    # will have to check N bp, if its a mismatch then $hd--
    # we also count the match that is the original location, so there's always 1 0mm entry
    #

    #we have to allow 6 mismatches total, as we don't count the N in the pam as a mismatch
    my $valid = 0;
    if ( $seq =~ /^CC/i && $rev_hd <= $MAX_EDIT_DISTANCE ) {
        $valid = 1;

        #compare the N in the pam sites; if they don't match then the edit distance
        #is reporting 1 more than we want it to (because mismatch in N doesn't count)
        if ( $rev_hd > 0 && substr($seq, 2, 1) ne substr($crisprs{$name}->{rev}, 2, 1) ) {
            #print STDERR "$seq\n" . $crisprs{$name}->{rev} . " has a mismatch in CCN.";
            $rev_hd--;
        }

        $hd_data{$exon_id}->{$crispr_id}{off_targets}{$rev_hd}++;
    }
    elsif ( $seq =~ /GG$/i && $fwd_hd <= $MAX_EDIT_DISTANCE ) {
        $valid = 1;

        #same as above, but the pam site is at the end.
        if ( $fwd_hd > 0 && substr($seq, 20, 1) ne substr($crisprs{$name}->{fwd}, 20, 1) ) {
            $fwd_hd--;
        }

        if ( $fwd_hd > 5 ) {
            print STDERR "$seq\n" . $crisprs{$name}->{fwd} . " has >5 mismatches??\n";
        }

        #we have to do this for both just in case
        $hd_data{$exon_id}->{$crispr_id}{off_targets}{$fwd_hd}++;
    }

    if ( $valid ) {
        #print STDERR "$seq looks like it has a valid pam\n";

        #print $line; #uncomment this to not put seqs in
        print "$chr\t$start\t$end\t$name-$seq\t$unknown\t$strand\n";
    }
    else {
        #print STDERR "Skipping $seq: invalid pam site.\n";
        #print STDERR " Hamming distances: ${rev_hd}, ${fwd_hd}\n";
    }
}

write_yaml_file();

print STDERR "Skipped a total of $num_n_lines N sequences:\n";
#print STDERR "$_\n" for @n_lines;

sub hamming_distance {
    #use string xor to get the number of mismatches between the two strings.
    #tr returns the number of changes it made, i.e. the number of letters that are different
    die "Strings passed to hamming distance differ" if length($_[0]) != length($_[1]);
    return (uc($_[0]) ^ uc($_[1])) =~ tr/\001-\255//;
}

sub write_yaml_file {
    #open the fq file (might need to use getopt to set $yaml_file)
    my $crisprs = LoadFile( $crispr_yaml_file );

    #create a json string of our mismatch data and add it to the crispr yaml
    for my $exon_id ( keys %hd_data ) {
        die "$exon_id isn't in $crispr_yaml_file file!" unless defined $crisprs->{$exon_id};
        while ( my ( $crispr_id, $crispr) = each %{ $hd_data{$exon_id} } ) {
            die "$crispr_id ($exon_id) isn't in $crispr_yaml_file!" unless defined $crisprs->{$exon_id}{$crispr_id};
            
            #display the off target information like: 0: 20, 1: 50, 2: 300, 3: 900
            my $ordered = join ", ", map { $_ . ": " . $crispr->{off_targets}{$_} } 
                                        sort { $a <=> $b } 
                                            keys %{ $crispr->{off_targets} };

            #make it the same format as before (not valid json, but it should be later)
            $crisprs->{$exon_id}{$crispr_id}{off_target_summary} = "{" . $ordered . "}";
        }
    }

    #write the yaml back out with the new off target information
    DumpFile( $crispr_yaml_file, $crisprs );
}

#this method takes a hashref and populates it with the fq data
sub parse_fq_file {
    my ( $crisprs ) = @_;

    open(my $fq_fh, "<", $crispr_fq_file) or die "Couldn't open $crispr_fq_file: $!";

    #parse the fq file. fq files can have a 3rd line which will break this but we generate them anyway
    my $current;
    while ( my $line = <$fq_fh> ) {
        chomp $line;

        if ( $line =~ /^@(.*)/ ) {
            $current = $1;
        }
        else {
            die "Error: there must only be one sequence per entry" if defined $crisprs->{ $current };
            die "Crisprs shouldn't be reverse complemented!" unless $line =~ /GG$/;

            #add fwd and reverse so we can compute hamming distances later.
            #we call it rev as we assume all crispr sequences have been reverse complemented.
            #cant automatically determine which is which cause CC-GG is a valid crispr
            $crisprs->{$current} = { rev => revcom( $line )->seq, fwd => $line };
        }
    }
}

sub fetch_ensembl_seq {
    my ( $chr, $start, $end ) = @_;

    #
    # we add 1 to start cause of the way bed files are. or maybe how bwa writes them idk
    #
    $start++;

    print STDERR "Fetching $chr:$start-$end\n";

    my $seq;

    try {
        $seq = $e->slice_adaptor->fetch_by_region("chromosome", $chr, $start, $end)->seq;
    }
    catch {
        print STDERR "Error fetching sequence, trying scaffold...\n";
    };

    #gross but whatever it works now and i don't use it anymore thanks to bed2fasta
    unless ( $seq ) {
        #must be scaffold.
        try {
            $seq = $e->slice_adaptor->fetch_by_region("scaffold", $chr, $start, $end)->seq;
        }
        catch {
            print STDERR "Couldn't find $chr:$start-$end on scaffold or chromosome\n";
        };
        
        unless ( $seq) {
            try {
                $seq = $e->slice_adaptor->fetch_by_region("supercontig", $chr, $start, $end)->seq;
            }
            catch {
                print STDERR "Couldn't find $chr:$start-$end on scaffold or chromosome or contig\n";
            };
        }

        #i give up
        next unless $seq;
    }

    return $seq;
}

1;

__END__

=head1 NAME

remove_invalid_crisprs.pl - remove potential off targets that don't have a valid pam site

=head1 SYNOPSIS

remove_invalid_crisprs.pl [options]

    --species            mouse or human
    --crispr-yaml-file   location of the crispr yaml data created by find_paired_crisprs.pl
    --fq-file            location of the crispr fq data created by find_paired_crisprs.pl
    --bed-file           location of bed file containing all potential off targets and their sequence
    --max-edit-distance  the maximum allowed edit distance for mismatches [optional, default 6]
    --help               show this dialog

Example usage:

remove_invalid_crisprs.pl --species mouse --fq-file /nfs/users/nfs_a/ah19/work/crisprs/hprt_crisprs.fq --crispr-yaml-file /nfs/users/nfs_a/ah19/work/crisprs/hprt_crisprs.yaml --bed-file /nfs/users/nfs_a/ah19/work/crisprs/HPRT1.with_seqs.bed > HPRT1.valid.bed

A bed file of only the valid crisprs from the source bed file is printed to STDOUT

=head1 DESCRIPTION

Find all possible crispr sites within exon sequences, and any possible pairs within them. 
The crispr-yaml file is created with the intention of being given to lims2-task to load into the database.
The fq file is intended to be handed to bwa.
It will also spit a csv out to stdout if you don't ask for a fq or yaml file

A hash is created of all the crisprs in the fq file in both directions, then every entry
in the bed file is checked against the sequence to: 
    1. check the hamming distance is below max-edit-distance
    2. see if there is a valid PAM site (i.e. does the off-target start CC or end GG)

If both of those conditions are true then it is a valid off-target, so will be printed to stdout
and the off-target information (i.e. the edit distance) will be stored in the crispr yaml file.

Again see paired_crisprs.sh for how I run this

=head AUTHOR

Alex Hodgkins

=cut