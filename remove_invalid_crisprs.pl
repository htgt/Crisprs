#!/usr/bin/perl

use strict;
use warnings;

use LIMS2::Util::EnsEMBL;
use Bio::Perl;
use Try::Tiny;
use Data::Dumper;

my $e = LIMS2::Util::EnsEMBL->new( species => "mouse" );

sub hamming_distance {
    #use string xor to get the number of mismatches between the two strings.
    #tr returns the number of changes it made, i.e. the number of letters that are different
    die "Strings passed to hamming distance differ" if length($_[0]) != length($_[1]);
    return (uc($_[0]) ^ uc($_[1])) =~ tr/\001-\255//;
}

#build a hash of all the crisprs so we can check hamming distances to be sure of orientation

die "Usage: remove_invalid_crisprs.pl <crisprs.fq> <bed_file.bed>" unless defined $ARGV[1];

#first file should be crisprs.fq, second the bed file
my $fq_file = shift;
open(my $fq_fh, "<", $fq_file) or die "Couldn't open $fq_file: $!";

my ( %crisprs, $current );
while ( my $line = <$fq_fh> ) {
    chomp $line;

    if ( $line =~ /^@(.*)/ ) {
        $current = $1;
    }
    else {
        die "Error: there must only be one sequence per entry" if defined $crisprs{ $current };
        #die "Crisprs should be reverse complemented!" unless $line =~ /^CC/;
        die "Crisprs shouldn't be reverse complemented!" unless $line =~ /GG$/;

        #add fwd and reverse so we can compute hamming distances later.
        #we call it rev as we assume all crispr sequences have been reverse complemented.
        #cant automatically determine which is which cause CC-GG is a valid crispr
        $crisprs{$current} = { rev => revcom( $line )->seq, fwd => $line };
    }
}

my $MAX_EDIT_DISTANCE = 6;



# die Dumper(%seqs);

# my %pams = (
#     "ARA"   => { pam => "CCA", revpam => "TGG" },
#     "ARB"   => { pam => "CCT", revpam => "AGG" },
#     "RAG11" => { pam => "CCA", revpam => "TGG" },
#     "RAG12" => { pam => "CCG", revpam => "CGG" },
#     "RAG1C" => { pam => "CCA", revpam => "TGG" },
#     "RAG1E" => { pam => "CCA", revpam => "TGG" },
#     "RAG1F" => { pam => "CCA", revpam => "TGG" },
#     "RAG1G" => { pam => "CCT", revpam => "AGG" },
# );

#take 2nd element as the pam (NGG), first is the bed file
#my $pam = pop @ARGV;
#my $revpam = revcom( $pam )->seq;

my @n_lines;

#print STDERR "Checking for $pam or $revpam\n";

while ( my $line = <> ) {
    chomp $line;
    #start is one too low so we +1 to it in every fetch_by_region
    my ( $chr, $start, $end, $name, $unknown, $strand ) = split /\s+/, $line;

    #die "No pam sites available for $name" unless exists $pams{$name};

    my $seq;
    if ( $name =~ /(.+)-([CTAGNctagn]{23})/ ) {
        #print STDERR "Split $name into $1 $2\n";
        $seq = $2;
        $name = $1;
    }
    else {
        print STDERR "$name has no seq in it, fetching from ensembl\n";
        $seq = fetch_ensembl_seq($chr, $start, $end);
    }


    if ( $seq =~ /N+/ ) {
        push @n_lines, "$line\t$seq\n";
        next;
    }

    #skip all the N sequences
    #next if $seq =~ /N+/;

    #if ( $seq !~ /^$revpam|$pam$/ ) {
    #    print "$seq didn't match $revpam or $pam\n";
    #}
    #else {
    #    print "$seq matches $revpam or $pam\n";
    #}

    #used the right way around here, unlike in the valid_pairs file. CONFUSING.
    #my ( $pam, $revpam ) = ( $pams{$name}->{revpam}, $pams{$name}->{pam} );

    #we're assuming reverse complemented pam, so we know orientation
    #check hamming distance to make sure the orientation is as we expect

    #we have to allow 6 mismatches total, as we don't count the N in the pam as a mismatch
    my $valid = 0;
    if ( $seq =~ /^CC/i && hamming_distance( $seq, $crisprs{$name}->{rev} ) <= $MAX_EDIT_DISTANCE ) {
        $valid = 1;
    }
    elsif ( $seq =~ /GG$/i && hamming_distance( $seq, $crisprs{$name}->{fwd} ) <= $MAX_EDIT_DISTANCE ) {
        $valid = 1;
    }

    #if ( $seq =~ /^$revpam|$pam$/ ) {
    if ( $valid ) {
        print STDERR "$seq looks like it has a valid pam\n";

        #print $line; #uncomment this to not put seqs in
        print "$chr\t$start\t$end\t$name-$seq\t$unknown\t$strand\n";
    }
    else {
        print STDERR "Skipping $seq: invalid crispr.";
        print STDERR " Hamming distances:" . hamming_distance( $seq, $crisprs{$name}->{rev} ) . ","
                    . hamming_distance( $seq, $crisprs{$name}->{fwd} ) . "\n";
    }

    #print $line if $seq =~ /^$revpam|$pam$/
    #die;
}

print STDERR "Skipped a total of " . scalar(@n_lines) . " N sequences:\n";
#print STDERR "$_\n" for @n_lines;

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
