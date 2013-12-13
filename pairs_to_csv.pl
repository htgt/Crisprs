#!/usr/bin/perl

use strict;
use warnings;

use feature qw( say );

use Getopt::Long;
use YAML::Any qw( LoadFile );

use Data::Dumper;

my $usage = "Usage: pairs_to_csv.pl --crispr-yaml-file crisprs.pl --pair-yaml-file pairs.pl > output.csv";

#should these all be caps? maybe. its really ugly though
my ( $crispr_yaml_file, $pair_yaml_file ); 

GetOptions(
    "crispr-yaml-file=s" => \$crispr_yaml_file,
    "pair-yaml-file=s"   => \$pair_yaml_file,
) or die $usage;

die $usage unless $crispr_yaml_file and $pair_yaml_file;

my $crispr_yaml = LoadFile( $crispr_yaml_file );
my $pair_yaml = LoadFile( $pair_yaml_file );

#use ensembl to get gene somehow? if we had real exon ids...

my @headers = qw(
    Exon
    Spacer
    Closest_Pair_OT
    L_Crispr_ID
    L_Crispr_Locus
    L_Crispr_Seq
    L_Crispr_OT
    R_Crispr_ID
    R_Crispr_Locus
    R_Crispr_Seq
    R_Crispr_OT
);

say join ",", @headers;

while ( my ( $exon_id, $pairs ) = each %{ $pair_yaml } ) {
    for my $pair ( @{ $pairs } ) {
        my ( $l_crispr, $r_crispr ) = get_crisprs( $exon_id, $pair );

        $pair->{off_target_summary} =~ tr/,/;/;

        say join ",", $exon_id,
                      $pair->{spacer},
                      $pair->{off_target_summary},
                      $l_crispr->{id},
                      $l_crispr->{locus_str},
                      $l_crispr->{seq},
                      $l_crispr->{off_target_summary},
                      $r_crispr->{id},
                      $r_crispr->{locus_str},
                      $r_crispr->{seq},
                      $r_crispr->{off_target_summary},
    }
}

sub get_crisprs {
    my ( $exon_id, $pair ) = @_;

    my ( $l, $r ) = ( $pair->{left_crispr}, $pair->{right_crispr} );

    for my $id ( $l, $r ) {
        my $crispr = $crispr_yaml->{$exon_id}{$id};
        die "Crispr $id isn't defined! Do you have the right files?"
            unless $crispr;

        $crispr->{id} = $id;

        $crispr->{locus_str} = $crispr->{locus}{chr_name} . ":" 
                             . $crispr->{locus}{chr_start} . "-" 
                             . $crispr->{locus}{chr_end};

        $crispr->{off_target_summary} =~ tr/,/;/;
    }

    return $crispr_yaml->{$exon_id}{$l}, $crispr_yaml->{$exon_id}{$r};

    #blah blab return the two crispr hashrefs
}