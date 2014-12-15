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

die "Usage: get_wge_crispr_ids.pl <species> <exon_ids.txt>" unless @ARGV >= 2;

#my $out_dir = dir( "/lustre/scratch109/sanger/ah19/human_flank_crispr_ids/" );
my $out_dir = dir( "/lustre/scratch109/sanger/ah19/grch38_crispr_ids/" );
#my $out_dir = dir( "/nfs/users/nfs_a/ah19/work/barry/" );

my ( $out_file, $fh );
my ( $file_num, $num, $limit, $num_exons, $total ) = ( 1, 0, 8000, 0, 0 );

my $model = WGE::Model::DB->new;

#keep a species object so it goes quicker
my $species = $model->resultset('Species')->find( { id => ucfirst( lc shift ) } );
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
    #only get crisprs without off targets
    my @ids = map { $_->id }
                grep { ! $_->off_target_summary }
                    $model->resultset('Exon')->find( $line )->crisprs( $species, 200 );
    $num_exons++;

    #add to counter
    $num += @ids;

    #write 1 id per line
    say $fh join "\n", @ids if @ids;
}

#if there arent enough to pass the limit total will be 0
$total = $num unless $total;

WARN "Found $total crisprs";


1;
