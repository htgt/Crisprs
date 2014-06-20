#!/usr/bin/env perl
use strict;
use warnings;

use WGE::Model::DB;
use Data::Dumper;
use Text::CSV;
use IO::Handle;
use IO::File;
use WGE::Util::FindPairs;
use Const::Fast;
use Log::Log4perl ':easy';
use feature qw( say );
use IPC::System::Simple qw( capture );

# hard coded species id of 2 - Mouse
const my $SPECIES_ID => 2;
Log::Log4perl->easy_init( { level => $DEBUG } );

LOGDIE "Usage: loxp_pairs.pl <loci.csv>" unless @ARGV;

my $w = WGE::Model::DB->new;
my $pair_finder = WGE::Util::FindPairs->new( schema => $w->schema );
$pair_finder->log->less_logging(2); #reduce output

my $csv = Text::CSV->new();
open ( my $fh, '<', $ARGV[0] ) or die( "Can not open $ARGV[0] " . $! );
# input csv file must have column headers
$csv->column_names( @{ $csv->getline( $fh ) } );

#
# run exon_coords from my .bash_functions to convert exon ids to 
# an appropriate format for this
#

my ( $line, $num ) = ( 0, 0 );
while ( my $data = $csv->getline_hr( $fh ) ) {
    say "Processing line " . ++$line;

    my $crisprs = crisprs_for_region( $data->{chr_name}, $data->{chr_start}-300, $data->{chr_end}+300 );
    my $pairs = crispr_pairs_for_region( $crisprs, $data->{chr_start}-300, $data->{chr_end}+300 );

    say "Found " . scalar( @{ $pairs } ) . " pairs";

    for my $pair ( @{ $pairs } ) {
        my $pair_id = $pair->{left_crispr}{id} . "_" . $pair->{right_crispr}{id};
        #say "Doing $pair_id";

        #pair crisprs is a resultset with just these two crisprs in,
        #which speeds up the data_missing call
        my ( $db_pair, $pair_crisprs ) = $w->find_or_create_crispr_pair( {
            left_id    => $pair->{left_crispr}{id}, 
            right_id   => $pair->{right_crispr}{id},
            species_id => $SPECIES_ID,
        } );

        #skip those with off target data
        if ( $db_pair->off_target_summary ) {
            WARN "Skipping $pair_id as it already has ots data";
            next;
        }

        if ( $db_pair->_data_missing( $pair_crisprs ) ) {
            WARN "$pair_id has missing individual data, skipping.";
            next;
        }
        else {
            $db_pair->calculate_off_targets;
        }

        say "Processed " . $num . " pairs" if ++$num % 100 == 0;
    }
}

=head crisprs_for_region

Find all the single crisprs in and around the target region.

=cut
sub crisprs_for_region {
    my ( $chr_name, $chr_start, $chr_end ) = @_;

    DEBUG("Getting crisprs for $chr_name:${chr_start}-${chr_end}");

    # we use 90 because the spaced between the crisprs in a pair can be 50 bases.
    # 50 + the size of 2 crisprs is around 90
    # that should bring back all the possible crisprs we want ( and some we do not want
    # which we must filter out )
    my @crisprs = $w->resultset('Crispr')->search(
        {
            'species_id'  => $SPECIES_ID,
            'chr_name'    => $chr_name,
            # need all the crisprs starting with values >= start_coord
            # and whose start values are <= end_coord
            'chr_start'   => {
                -between => [ $chr_start - 90, $chr_end + 90 ],
            },
        },
    )->all;

    return \@crisprs;
}

=head crispr_pairs_for_region

Identifies valid pairs within the list of crisprs for the region

=cut
sub crispr_pairs_for_region {
    my ( $crisprs, $loxp_start, $loxp_end ) = @_;

    # Find pairs amongst crisprs
    my $pairs = $pair_finder->find_pairs( $crisprs, $crisprs, { species_id => $SPECIES_ID, get_db_data => 0  } );

    return validate_crispr_pairs( $pairs, $loxp_start, $loxp_end );
}

=head validate_crispr_pairs

The crispr pair is valid if one or both of its crisprs lies within
the target region ( loxp site )

=cut
sub validate_crispr_pairs {
    my ( $pairs, $start, $end ) = @_;
    my @validated_pairs;

    for my $pair ( @{ $pairs } ) {
        if (   crispr_in_target_region( $pair->{left_crispr}, $start, $end )
            || crispr_in_target_region( $pair->{right_crispr}, $start, $end ) )
        {
            push @validated_pairs, $pair;
        }
    }

    return \@validated_pairs;
}

=head crispr_in_target_region

Crispr is in target region if the at least one base of the pam site is within
the target region.

=cut
sub crispr_in_target_region {
    my ( $crispr, $start, $end ) = @_;

    if ( $crispr->{pam_right} ) {
        if ( $crispr->{chr_end} > $start && $crispr->{chr_end} <= $end ) {
            return 1;
        }
    }
    else {
        if ( $crispr->{chr_start} >= $start && $crispr->{chr_start} < $end ) {
            return 1;
        }
    }

    return;
}

1;

=head1 NAME

loxp_pairs.pl - find groups of crisprs pairs in loxp regions

=head1 SYNOPSIS

promotor_crisprs_report.pl [input_file]

Find crispr pairs for loxp regions of conditional mouse designs

Example usage:

loxp_pairs.pl [input_file]

Input must be a csv file which has the following columns ( csv file must have column headers ):
loxp_start
loxp_end
chr_name

=head1 DESCRIPTION

Finds crispr pairs that cut in the wildtype sequence of the loxp region of conditional designs ( between D5 and D3 ).
Input csv file must have loxp_start, loxp_end and chr_name columns to work out what crisprs hit the region.

Two output files:
 -loxp_crispr_report.csv: copy of original csv file, plus extra column added with crispr / crispr pair ids and counts
 -loxp_crispr_ids.txt: A list of all the crispr ids

=cut

