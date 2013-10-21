#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use LIMS2::Model;
use YAML::Any qw( LoadFile DumpFile );
use Try::Tiny;

my ( $species, $assembly, @yaml_files );
my $commit = 0;
GetOptions(
    "help"               => sub { pod2usage( 1 ) },
    "man"                => sub { pod2usage( 2 ) },
    "species=s"          => sub { my ( $name, $val ) = @_; $species = ucfirst(lc $val); },
    "assembly=s"         => \$assembly, #e.g. GRCh37
    "crispr-yaml=s{,}"   => \@yaml_files,
    "commit!"            => \$commit,
) or pod2usage( 2 );

die pod2usage( 2 ) unless $species and $assembly and @yaml_files;

my $model = LIMS2::Model->new( user => 'tasks' );

for my $filename ( @yaml_files ) {
    die "$filename doesnt exist!" unless -f $filename;

    my $data = LoadFile( $filename ) || die "Couldn't load $filename";

    my $failed = 0;
    
    while ( my ( $exon_id, $crisprs ) = each %{ $data } ) {
        #we don't care about the crispr id so we just get the values
        for my $crispr ( values %{$crisprs} ) {
            my $crispr_id = delete $crispr->{id}; #db doesnt want this

            $crispr->{species} = $species;
            $crispr->{locus}{assembly} = $assembly;
            $crispr->{comment} = $crispr_id; #no harm in storing this

            try {
                $model->txn_do(
                    sub {
                        my $db_crispr = $model->create_crispr( $crispr );
                        if ( $commit ) {
                            $crispr->{db_id} = $db_crispr->id;
                        }
                        else {
                            $model->txn_rollback;
                        }
                    }
                );
            }
            catch {
                print STDERR "$_\n";
                $failed = 1;
            };

            $crispr->{id} = $crispr_id; #put the id back in
            delete $crispr->{comment}; #no point writing this to disk
        }
    }

    #now we've added all the ids write them back into the file. 
    try {
        return if $failed;
        #print "Would be writing data here. Num keys: " . scalar(keys %$data) . "\n";
        DumpFile( $filename, $data ) || die "Couldn't write to $filename!";
    }
    catch {
        #if it fails print to stduot
        print STDERR "$_\n";
        print STDERR "Couldn't write yaml to $filename, printing to stdout\n";
        print Dump( $data );
        $failed = 1;
    };

    if ( $failed ) {
        die "Something failed, skipping further yaml files.";
    }
}

1;

__END__

=head1 NAME

persist_crisprs.pl - insert crisprs into the db and write the database ids into the file

=head1 SYNOPSIS

persist_crisprs.pl [options]

    --species            mouse or human
    --assembly           the assembly used to generate the crisprs
    --crispr-yaml        yaml file(s) with crispr data in them
    --commit             whether or not to persist to the db. defaults to false
    --help               show this dialog

Example usage:

persist_crisprs.pl --species human --assembly GRCh37 --crispr-yaml HPRT1/HPRT1_crisprs.yaml CBX1/CBX1_crisprs.yaml --commit

=head1 DESCRIPTION

Find all possible crispr sites within exon sequences, and any possible pairs within them. 
The crispr-yaml file is created with the intention of being given to lims2-task to load into the database.
The fq file is intended to be handed to bwa.
It will also spit a csv out to stdout if you don't ask for a fq or yaml file

See paired_crisprs.sh for how I run this bad boy

Insert all crisprs in the yaml into the database, and write the inserted crispr id back into the yaml file.
If db_id is already set then the crispr will be skipped

=head AUTHOR

Alex Hodgkins

=cut