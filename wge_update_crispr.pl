#!/usr/bin/perl

use strict;
use warnings;

use Log::Log4perl qw( :easy );
use autodie;
use feature qw( say );
use Try::Tiny;
use Path::Class;
use List::MoreUtils qw( natatime );

use WGE::Model::DB;

BEGIN { Log::Log4perl->easy_init( $DEBUG ) }

die "Usage: run_batch_crisprs.pl" unless @ARGV >= 1;

my $out_dir = dir( "/lustre/scratch109/sanger/ah19/crispr_ids/" );

my ( $batch_size = 100 );
my @skip;

my $model = WGE::Model::DB->new;

#keep a species object so it goes quicker
my $species = $model->resultset('Species')->find( 1 );
WARN "Species is " . $species->id;

for my $file ( $out_dir->children ) {
    DEBUG "Processing $file";

    next if grep { $file->stringify =~ $_ } @skip;

    #we do the whole file in a transaction to avoid lots of commits
    try {
        $model->txn_do( sub {
            my @all_ids = $file->slurp( chomp => 1 );

            #use this terribly named method to fetch 100 crisprs at a time
            my $it = natatime $batch_size, @all_ids;

            #update the batch
            my $i = 0;
            while ( my @ids = $it->() ) {
                DEBUG "Batch " . ++$i;

                $model->resultset('Crispr')->search( 
                    { id => { -IN => \@ids } } 
                )->update( { genic => 1, exonic => 1 } );
            }
        } );
    }
    catch {
        WARN "Died commiting file $file!";
        die $_;
    };

    die "DEATH";
}