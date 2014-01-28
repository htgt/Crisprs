#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use LIMS2::REST::Client;
use YAML::Any qw( LoadFile DumpFile );
use Try::Tiny;
use Log::Log4perl qw(:easy);
use Data::Dumper;

Log::Log4perl->easy_init($DEBUG);

my ( $species, @crisprs, @pairs );
my $commit = 0;
GetOptions(
    "help"               => sub { pod2usage( 1 ) },
    "man"                => sub { pod2usage( 2 ) },
    "species=s"          => sub { my ( $name, $val ) = @_; $species = ucfirst(lc $val); },
    "crispr-yaml=s{,}"   => \@crisprs,
    "pair-yaml=s{,}"     => \@pairs,
    "commit!"            => \$commit,
) or pod2usage( 2 );

die pod2usage( 2 ) unless $species and @crisprs and @pairs;

my $client = LIMS2::REST::Client->new_with_config(configfile => $ENV{LIMS2_REST_CLIENT_CONFIG});

for my $crispr_yaml ( @crisprs ) {
    my $pair_yaml = shift @pairs;

    my $crispr_data = LoadFile( $crispr_yaml );
    my $pair_data = LoadFile( $pair_yaml );

    while ( my ( $exon_id, $pairs ) = each %{ $pair_data } ) {
        for my $pair ( @{ $pairs } ) {
            my $l_crispr = $crispr_data->{$exon_id}{$pair->{left_crispr}};
            my $r_crispr = $crispr_data->{$exon_id}{$pair->{right_crispr}};

            my $failed = 0;

            if ( $exon_id eq "ENSExistingCrisprs" ) {
                #print STDERR "Updating OTS instead of creating a new pair\n";
                my ( $l_id ) = $l_crispr->{id} =~ /(\d+)/;
                my ( $r_id ) = $r_crispr->{id} =~ /(\d+)/;
                
                if($commit){
                    try{
                        my $updated_pair = $client->POST("crispr_pair_off_target_summary",
                        {
                            l_id               => $l_id,
                            r_id               => $r_id,
                            off_target_summary => $pair->{off_target_summary},
                        });
                

                        print STDERR "Updated " . $updated_pair->{id} . " with " 
                                       . $updated_pair->{off_target_summary} . "\n";
                    }
                    catch {
                        print STDERR "Update failed: $_\n";
                        $failed = 1;
                    };
                }

                next;
            }

            try {
                #make sure the ids match up
                if ( ! defined $l_crispr || ! defined $r_crispr ) {
                    die "Crispr and pair ids don't match -- have you put the right files in??";
                }

                #make sure this entry has db ids
                if ( ! defined $l_crispr->{db_id} || ! defined $r_crispr->{db_id} ) {
                    die "You must insert the crispr data into the db first!";
                }

                if ($commit){
                    my $new_pair = $client->POST("crispr_pair", 
                    {
                        l_id => $l_crispr->{db_id},
                        r_id => $r_crispr->{db_id},
                        spacer => $pair->{spacer},
                        off_target_summary => $pair->{off_target_summary},
                    });
                    $pair->{db_id} = $new_pair->{id};
                    print STDERR "Created or updated crispr pair with id ".$new_pair->{id}."\n";
                }
            }
            catch {
                print STDERR "$_ ($crispr_yaml & $pair_yaml)";
                print STDERR "Unsuccessful commit: $exon_id-" . $l_crispr->{id} . " and " . $r_crispr->{id} . "\n";
                $failed = 1;
            };

            die "Broken" if $failed;
        }
    }
}

1;

__END__

=head1 NAME

persist_pairs.pl - insert pairs into the db and write the database ids into the file

=head1 SYNOPSIS

persist_pairs.pl [options]

    --species            mouse or human
    --crispr-yaml        yaml file with crispr data (and db ids) in them
    --pair-yaml          yaml file containing the pair data
    --commit             whether or not to persist to the db. defaults to false
    --help               show this dialog

Example usage:

persist_pairs.pl --species human --crispr-yaml CBX1/CBX1_crisprs.yaml --pair-yaml CBX1/CBX1_pairs.yaml --commit

=head1 DESCRIPTION

Insert all pairs from the yaml file into the database, assuming all the crisprs referenced 
in the pairs already exist. The db id will then be added to the entry. 
If a db_id already exists it will be overwritten.

=head AUTHOR

Alex Hodgkins

=cut