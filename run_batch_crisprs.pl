#!/usr/bin/perl

use strict;
use warnings;

use Log::Log4perl qw( :easy );

use IPC::System::Simple qw( run );
use Try::Tiny;
use File::Path qw( remove_tree );

BEGIN { Log::Log4perl->easy_init( $DEBUG ) }

die "Usage: run_batch_crisprs.pl <species> <exons.txt>" unless @ARGV >= 2;

my $species = shift;
my $failed_file = "/lustre/scratch109/sanger/ah19/crispr_logs/failed.txt";
my ( $num_failed, $num_passed ) = ( 0, 0 );

DEBUG "Processing file " . $ARGV[0] . " for $species";

while( my $line = <> ) {
    chomp $line;

    my $filename;
    if ( $line =~ /^ENS/ ) {
        if ( $species eq "Human" ) {
            die "Non human exon id provided!" unless $line =~ /ENSE000/;
        }
        else {
            die "Non mouse exon id provided!" unless $line =~ /ENSMUSE000/;
        }
        $filename = $line;
    }
    elsif ( $line =~ /^\d+/ ) { #allow numbers
        DEBUG "Got " . scalar( split /\s+/, $line ) . " crispr ids:";
        #create a filename from truncated ids
        ( $filename = substr $line, 0, 40 ) =~ s/\s+/_/g;
        DEBUG "Filename is $filename";
    }
    else {
        die "Unrecognised id provided: $line\n";
    }

    #actual filenames get added depending on the type
    my $dir = "/lustre/scratch110/sanger/ah19/wge_exons/$filename";
    my $log_file = "/lustre/scratch109/sanger/ah19/crispr_logs/auto_$filename.txt";

    try {
        DEBUG "Creating directory $dir";
        mkdir $dir || die "Couldn't make $dir";

        DEBUG "Running crispr finding on $line";
        
        #my $cmd = "echo $line > $log_file";
        my $cmd = "crisprs_for_exon_wge.sh $dir $species $line > $log_file";
       
        DEBUG "Running $cmd";

        #this raises an exception on non 0 retval
        my $retval = run $cmd;

        DEBUG "Returned $retval";
        #if we got here everything was successful, 
        unlink $log_file;

        $num_passed++;
    }
    catch {
        WARN "$line failed, skipping";
        open my $fh, ">>", $failed_file;
        print $fh "$line\n";
        
        $num_failed++;
    };
    
    DEBUG "Deleting $dir";

    #always delete the folder
    my $deleted;
    remove_tree( $dir, { result => \$deleted } );

    WARN "Deleted " . scalar( @{ $deleted } ) . " files.";

    WARN "Processed " . ( $num_passed+$num_failed ) . " exons";
}

WARN "All exons complete:";
WARN "$num_passed jobs passed";
WARN "$num_failed jobs failed";

if ( $num_failed ) {
    WARN "Failed exons can be seen in this file: $failed_file";
}

1;
