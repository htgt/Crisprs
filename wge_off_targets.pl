#!/usr/bin/perl

use strict;
use warnings;

use Log::Log4perl qw( :easy );
use UUID::Tiny ':std';
use autodie;

use IPC::System::Simple qw( run );
use Try::Tiny;
use File::Path qw( remove_tree );
use Path::Class;

use WGE::Model::DB;

BEGIN { Log::Log4perl->easy_init( $DEBUG ) }

die "Usage: run_batch_crisprs.pl <species> <ids.txt>" unless @ARGV >= 2;

my $batch_size = 100;
my @batch;
my $run_dir = dir( "/lustre/scratch110/sanger/ah19/wge_exons/" );
my $log_dir = dir( "/lustre/scratch109/sanger/ah19/crispr_logs/" );
my $index_file = file( "/lustre/scratch109/sanger/ah19/crisprs_human.bin" );

die "$index_file doesn't exist" unless -f $index_file;

my $model = WGE::Model::DB->new;

my $species = shift;
my $num_passed = 0;

DEBUG "Processing file " . $ARGV[0] . " for $species";

while( my $line = <> ) {
    chomp $line;

    die "ids must be numeric" unless $line =~ /\d+/;

    push @batch, $line;

    #wait until we have enough
    next if @batch < $batch_size;

    #use uuids because i'm lazy
    my $id = create_uuid_as_string();

    my $dir = $run_dir->subdir( $id );
    $dir->mkpath();

    #my $log_file = $log_dir->file( $id );

    #dump the ids into a text file
    {
        my $fh = $dir->file( "ids.txt" )->openw();
        print $fh join "\n", @batch;
    }

    DEBUG "Running batch with " . scalar( @batch ) . " ids, $id:";
    DEBUG join ", ", @batch;

    my $outfile = $dir->file('output.tsv');
    my $failed = 0;

    try {
        my $cmd = join " ", (
            "~/work/paired_crisprs/cpp/off_targets/find_off_targets",
            "align",
            "-i", $index_file->stringify,
            @batch,
            '>', $outfile->stringify,
        );

        DEBUG "Running $cmd";

        my $retval = run $cmd;

        DEBUG "returned $retval";

        #unlink $log_file->stringify;

        $num_passed++;
    }
    catch {
        WARN $_;
        $failed = 1;
    };

    #reset for next batch
    @batch = ();

    unless ( $failed ) {
        DEBUG "Verifying output against the db";

        my $fh = $outfile->openr();

        while ( my $line = <$fh> ) {
            chomp $line;

            my %fields;
            @fields{qw(id species_id off_target_ids off_target_summary)} = split /\t/, $line;

            if ( $fields{off_target_ids} eq 'NULL' ) {
                $fields{off_target_ids} = undef;
            }
            else {
                $fields{off_target_ids} = [ split /,/, substr($fields{off_target_ids}, 1, -1) ];
            }

            DEBUG "Processing crispr " . $fields{id};

            my $crispr = $model->resultset('Crispr')->find( $fields{id} );
            unless ( $crispr ) {
                WARN "Couldn't find crispr";
                next;
            }

            if ( $crispr->off_target_summary eq $fields{off_target_summary} ) {
                WARN "Summaries match";
            }
            else {
                WARN "Summaries don't match!";
                WARN "db: '" . $crispr->off_target_summary . "'";
                WARN "fi: '" . $fields{off_target_summary} . "'";
            }

            if ( ! $fields{off_target_ids} && ! $crispr->off_target_ids ) {
                WARN "Both have null off target arrays, so match.";
                next;
            }
            else {
                WARN "Old has: " . scalar( @{ $fields{off_target_ids} } );
                WARN "New has: " . scalar( @{ $crispr->off_target_ids } );
            }

            for my $new ( @{ $fields{off_target_ids} } ) {
                if ( ! grep { $new == $_ } @{ $crispr->off_target_ids } ) {
                    die "$new found in new and was NOT in old!";
                }
            }

            WARN "All new ids in old";

            #verify those missing from the db are CCGG
            my @missing;
            for my $old ( @{ $crispr->off_target_ids } ) {
                if ( ! grep { $old == $_ } @{ $fields{off_target_ids} } ) {
                    push @missing, $old;
                }
            }

            WARN "Fetching " . scalar( @missing ) . " off targets to verify";

            next unless @missing;

            #
            for my $c ( $model->resultset('Crispr')->search( { id => { -IN => \@missing} } ) ) {
                if ( $c->seq !~ /CC.*GG/ ) {
                    WARN $c->id . ": " . $c->seq;
                    die "Found off target that is NOT CC-GG!";
                }
            }
        }
    }
    
    DEBUG "Deleting $dir";
    #always delete the folder
    my $deleted = $dir->rmtree();

    WARN "Deleted $deleted files.";
}

if ( @batch ) {
    DEBUG "Stuff remaining in batch:";
    DEBUG "@batch";
}

WARN "All crisprs complete.";
WARN "$num_passed jobs passed";

1;
