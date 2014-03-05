#!/usr/bin/perl

use strict;
use warnings;

use Log::Log4perl qw( :easy );

use IPC::System::Simple qw( run );
use Try::Tiny;

BEGIN { Log::Log4perl->easy_init( $DEBUG ) }

die "Usage: run_batch_crisprs.pl <species> <exons.txt>" unless @ARGV >= 2;

my $species = shift;
my $failed_file = "/lustre/scratch109/sanger/ah19/crispr_logs/failed.txt";
my ( $num_failed, $num_passed ) = ( 0, 0 );

DEBUG "Processing file " . $ARGV[0] . " for $species";

while( my $line = <> ) {
    chomp $line;

    if ( $species eq "Human" ) {
        die "Non human exon id provided!" unless $line =~ /ENSE000/;
    }
    else {
        die "Non mouse exon id provided!" unless $line =~ /ENSMUSE000/;
    }

    my $dir = "/lustre/scratch110/sanger/ah19/wge_exons/$line";
    my $log_file = "/lustre/scratch109/sanger/ah19/crispr_logs/auto_$line.txt";

    try {
        DEBUG "Creating directory $dir";
        mkdir $dir || die "Couldn't make $dir";

        DEBUG "Running crispr finding on $line";

        #run("crisprs_for_exon_wge.sh auto_$line $species $line > $log_file");
        
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
    rmdir $dir;

    WARN "Processed " . ( $num_passed+$num_failed ) . " exons";
}

WARN "All exons complete:";
WARN "$num_passed jobs passed";
WARN "$num_failed jobs failed";

if ( $num_failed ) {
    WARN "Failed exons can be seen in this file: $failed_file";
}

1;
