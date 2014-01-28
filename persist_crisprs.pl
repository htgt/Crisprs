#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use LIMS2::REST::Client;
use YAML::Any qw( LoadFile DumpFile Dump );
use Try::Tiny;
use Log::Log4perl qw(:easy);
use Data::Dumper;

Log::Log4perl->easy_init($DEBUG);

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

my $client = LIMS2::REST::Client->new_with_config(configfile => $ENV{LIMS2_REST_CLIENT_CONFIG});

for my $filename ( @yaml_files ) {
    unless ( -f $filename ) {
        print STDERR "$filename doesnt exist, skipping.";
        next;
    }

    my $data = LoadFile( $filename ) || (print "Couldn't load $filename, skippnig" && next);

    my $failed = 0;
    
    while ( my ( $exon_id, $crisprs ) = each %{ $data } ) {
        #we don't care about the crispr id so we just get the values
        print STDERR "Exon ID: $exon_id\n";
        for my $crispr ( values %{$crisprs} ) {
            #skip a crispr if it has already been entered (in case some fail and some don't)
            #next if defined $crispr->{db_id};

            #this is used for an existing crispr we just want to re-calculate OTS for
            if ( $exon_id eq "ENSExistingCrisprs" ) {
                print STDERR "Updating OTS instead of creating a new crispr\n";

                my ( $db_id ) = $crispr->{id} =~ /(\d+)/;

                print STDERR "Updating " . $db_id . " with " 
                             . $crispr->{off_target_summary} . "\n";

                if ( $commit ){
                    try{
                        my $summary = $client->POST("crispr_off_target_summary", 
                        {
                            id                 => $db_id,
                            algorithm          => 'bwa',
                            off_target_summary => $crispr->{off_target_summary},
                        });
                    }
                    catch{
                        print STDERR "Update failed: $_";
                        $failed = 1;
                    };

                }

                next;
            }

            delete $crispr->{db_id} if defined $crispr->{db_id} and $commit;

            my $crispr_id = delete $crispr->{id}; #db doesnt want this

            $crispr->{species} = $species;
            $crispr->{locus}{assembly} = $assembly;
            $crispr->{comment} = $crispr_id; #no harm in storing this

            #
            #TODO make transaction over the whole method, including the writing of the file
            # - if it fails we want to rollback
            #

            if ($commit){
                print STDERR "Creating crispr: ".Dumper($crispr);
                try{
                    my $new_crispr = $client->POST("single_crispr", $crispr);

                    print STDERR "Crispr created with id ".$new_crispr->{id}."\n";
                    $crispr->{db_id} = $new_crispr->{id};
                }
                catch{
                    print STDERR "Crispr creation failed: $_";
                    $failed = 1;
                };
            }

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
        #print Dump( $data );
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

Insert all crisprs in the yaml into the database, and write the inserted crispr id back into the yaml file.

=head AUTHOR

Alex Hodgkins

=cut
