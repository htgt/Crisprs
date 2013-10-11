#!/usr/bin/perl

use strict;
use warnings;

use Bio::Perl;
use Try::Tiny;

my %labels;
my $MAX_OFF_TARGETS = 8000;

die "Usage: split_bed.pl <bed_file> <print_windowbed_args?>" unless defined $ARGV[0];

#we'll append the label to the file. strip .bed extension
my $filename = substr( $ARGV[0], 0, -4 );

#permutations is for removing duplicates
my ( $print_windowbed, %permutations );
if ( defined $ARGV[1] ) {
    $print_windowbed = 1;
    pop @ARGV;
}

while ( my $line = <> ) {
    #chr, start, end, name, unknown, strand
    my $name = (split /\s+/, $line)[3];

    #if its an ARA-CGTCGACTAGTACG one then strip the seq
    if ( $name =~ /(.+)-[CGATcgat]+/ ) {
        $name = $1;
    }

    push @{ $labels{$name} }, $line;
}

#write everything out to separate files
for my $name ( keys %labels ) {
    my $out_filename = "$filename-$name.bed";
    open( my $fh, ">", $out_filename);
    print $fh $_ for @{ $labels{$name} };

    print STDERR "Total entries for $out_filename: " . scalar(@{ $labels{$name} }) . "\n";
    
    if ( $print_windowbed ) {
        print_windowbed_args( $name );
    }
}

sub print_windowbed_args {
    my ( $current ) = @_;

    if ( @{ $labels{$current} } > $MAX_OFF_TARGETS ) {
        print STDERR "Not adding $current to permutations, as it has too many off targets.";
        return;
    }

    #allow ENSE/ENSMUSE/whatever species
    my $re = qr/(ENS(?:\w+)?E\d+)_(\d+\w)/;

    my ( $exon_id, $crispr_id ) = $current =~ $re;

    for my $label ( keys %labels ) {
        #only print an entry if the exon ids match, no point checking pairs in different exons
        if ( $label =~ $re ) {
            next unless $1 eq $exon_id;
            next if defined $permutations{ "${exon_id}_${crispr_id}_$2" };

            #we store it backwards because that's how it will be found next time.
            $permutations{ "${exon_id}_${2}_${crispr_id}" } = 1;
            print "$current $label\n";
        }
    }    
}
