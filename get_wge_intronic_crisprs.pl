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

my $out_dir = dir( "/lustre/scratch109/sanger/ah19/intronic_crispr_ids/" );

my ( $out_file, $fh );
my ( $file_num, $num, $limit, $num_genes, $total ) = ( 1, 0, 8000, 0, 0 );

my $model = WGE::Model::DB->new;

#keep a species object so it goes quicker
my $species = $model->resultset('Species')->find( 1 );
WARN "Species is " . $species->id;

my @genes = $model->resultset('Gene')->search();

for my $gene ( @genes ) {

    #either first time or new file
    if ( ! defined $fh || $num > $limit ) {
        #we have a file open so print the information
        if ( defined $fh ) {
            $total += $num;
            WARN "Wrote " . $num . " crisprs from " . $num_genes . " genes" if $num > 0;
            close( $fh ); #close old filehandle
            die "Death";
        }

        $num = 0; #reset counter
        $num_genes = 0;

        #open new file

        my $filename = "crisprs_" . $file_num++;
        WARN "Opening $filename";
        $fh = $out_dir->file( $filename )->openw;
    }

    my $id = $gene->ensembl_gene_id;

    #get crisprs with no off target summary within the gene
    my @ids = map { $_->id } 
                $model->resultset('CrisprByGene')->search( 
                    { off_target_summary => undef  },
                    { bind => [ $id, $species->numerical_id ] }
                );

    $num_genes++;

    #add to counter
    $num += @ids;

    #write 1 id per line
    say $fh join "\n", @ids;
}

WARN "Found $total crisprs";


1;
