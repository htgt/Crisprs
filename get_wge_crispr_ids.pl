#!/usr/bin/perl

use strict;
use warnings;

use Log::Log4perl qw( :easy );
use autodie;
use feature qw( say );
use Try::Tiny;
use Path::Class;

use WGE::Model::DB;

BEGIN { Log::Log4perl->easy_init( $DEBUG ) }

die "Usage: get_wge_crispr_ids.pl <exon_ids.txt>" unless @ARGV >= 1;

my $out_dir = dir( "/lustre/scratch110/sanger/ah19/mouse_crispr_ids/" );

my ( $out_file, $fh );
my ( $file_num, $num, $limit, $num_exons, $total ) = ( 1, 0, 8000, 0, 0 );

my $model = WGE::Model::DB->new;

#keep a species object so it goes quicker
my $species = $model->resultset('Species')->find( 1 );
WARN "Species is " . $species->id;

while ( my $line = <> ) {
    chomp $line;

    #either first time or new file
    if ( ! defined $fh || $num > $limit ) {
        #we have a file open so print the information
        if ( defined $fh ) {
            $total += $num;
            WARN "Wrote " . $num . " crisprs from " . $num_exons . " exons" if $num > 0;
            close( $fh ); #close old filehandle
        }

        $num = 0; #reset counter
        $num_exons = 0;

        #open new file

        my $filename = "crisprs_" . $file_num++;
        WARN "Opening $filename";
        $fh = $out_dir->file( $filename )->openw;
    }

    #get the associated crisprs for this exon
    #species id human with 200 base flanks
    my @ids = map { $_->id } $model->resultset('Exon')->find( $line )->crisprs( $species, 200 );
    $num_exons++;

    #add to counter
    $num += @ids;

    #write 1 id per line
    say $fh join "\n", @ids;
}

WARN "Found $total crisprs";


1;
