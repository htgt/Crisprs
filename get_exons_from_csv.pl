#!/usr/bin/env perl

use strict;
use warnings;
use 5.010; #needed for say

use feature qw( say );
use List::MoreUtils qw( first_index );

die "You must provide a csv file" unless @ARGV == 1;

my %genes;
while ( my $line = <> ) {
    chomp $line;
    my @cols = split ",", $line;

    my ( $marker_symbol, $exon_id ) = @cols[0, 2];

    if ( ! $exon_id ) {
        if ( $cols[3] ) {
            say STDERR "No exon id found, setting it to " . $cols[3];
            $exon_id = $cols[3]; #mt4_exon column, where mark manually chose an exon
        }
        else {
            say STDERR "$marker_symbol has no exon id ($exon_id) ($line)";
            next;
        }
    }

    push @{ $genes{$marker_symbol} }, $exon_id;
}

say STDERR "Extracted " . scalar( keys %genes ) . " genes.";

#print them in a format that's easy for bash to process
while ( my ( $marker_symbol, $exon_ids ) = each %genes ) {
    say $marker_symbol . " " . join " ", @{ $exon_ids };
}