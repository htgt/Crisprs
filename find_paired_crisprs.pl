#!/usr/bin/perl

use strict;
use warnings;

use feature qw( say );

use Getopt::Long;
use Pod::Usage;
use LIMS2::Util::EnsEMBL;
use Bio::Perl qw( revcom );
use YAML::Any qw( DumpFile );

use Data::Dumper;

#should these all be caps? maybe. its really ugly though
my ( $species, $fq_file, $yaml_file, @exon_ids );
my ( $expand_seq, $MIN_SPACER, $MAX_SPACER ) = ( 1, -10, 50 ); #set default values
GetOptions(
    "help"          => sub { pod2usage( 1 ) },
    "man"           => sub { pod2usage( 2 ) },
    "species=s"     => \$species,
    "exon-ids=s{,}" => \@exon_ids,
    "fq-file=s"     => \$fq_file,
    "yaml-file=s"   => \$yaml_file,
    "expand-seq!"   => \$expand_seq,
    "min-spacer=i"  => \$MIN_SPACER,
    "max-spacer=i"  => \$MAX_SPACER,
) or pod2usage( 2 );

die pod2usage( 2 ) unless $species and @exon_ids;

my $exon_re;
if ( $species =~ /mouse/i ) {
    $exon_re = qr/ENSMUSE/;
}
elsif ( $species =~ /human/i ) {
    $exon_re = qr/ENSE/;
}
else {
    die "Unknown species: $species";
}

my $e = LIMS2::Util::EnsEMBL->new( species => $species );

#
# exons with gibson designs:
# ENSE00003596810 ENSE00003663376 ENSE00003606976 ENSE00003682156 ENSE00003611316 ENSE00001427230 ENSE00001083354 ENSE00000683019 ENSE00000849850 ENSE00003701380 ENSE00003596529 
#

for my $exon_id ( @exon_ids ) {
    if ( $exon_id !~ $exon_re ) {
        say STDERR "$exon_id is not a valid exon for this species!";
        next;
    }

    say STDERR "Finding crisprs in $exon_id";

    #get the exon and 200bp surrounding if expand seq is set
    my $exon_slice = $e->slice_adaptor->fetch_by_exon_stable_id( $exon_id );

    #my ( $seq ) = $exon_slice->seq;

    if ( $expand_seq ) {
        my $orig_length = $exon_slice->length;
        say STDERR "Original length: " . $orig_length;
        #we need the gene so we can get the strand and take an asymmetrical slice:
        # 5' -----[ EXON ]----- 3'
        #      200        100
        my $strand = $e->gene_adaptor->fetch_by_exon_stable_id( $exon_id )->strand;

        #expand the slice considering the strand.
        if ( $strand == 1 ) {
            $exon_slice = $exon_slice->expand( 200, 100 );
        }
        elsif ( $strand == -1 ) {
            $exon_slice = $exon_slice->expand( 100, 200 );
        }
        else {
            die "Unexpected strand for gene associated with $exon_id";
        }
    }

    die "Couldn't get slice" unless $exon_slice;

    #say $exon_slice->seq, " (" . length($exon_slice->seq) . ")";
    #my @matches = get_matches( $exon_slice->seq );
    my @matches = get_matches_fixed( $exon_slice );
    #say "Found " . scalar( @matches ) . " paired crisprs:";

    write_fq( \@matches, $exon_id ) if defined $fq_file;
    write_yaml( \@matches, $exon_id ) if defined $yaml_file;

    #this is how it worked before the fq/yaml option so recreate that
    if ( ! defined $fq_file and ! defined $yaml_file ) {
        write_csv( \@matches, $exon_id );
    }

}

sub write_fq {
    my ( $pairs, $exon_id ) = @_;

    die "fq output file location not defined!" unless defined $fq_file;
    open( my $fh, ">", $fq_file ) || die "Couldn't open $fq_file for writing.";

    #print out all unique crisprs in fq format. if it ends in A it's a CC crispr and needs to be reverse
    #complemented for bwa

    my $unique_crisprs = _get_unique_crisprs( $pairs );
    #sort by integer part of id, ignoring A/B
    for my $crispr_id ( sort { ($a =~ /(\d+)[AB]$/)[0] <=> ($b =~ /(\d+)[AB]$/)[0] } keys %{ $unique_crisprs } ) {
        my $seq = $unique_crisprs->{$crispr_id}{seq};
        if ( $crispr_id =~ /A$/ ) {
            #should we insert it into the db reverse complemented??
            $seq = revcom( $seq )->seq;    
        }

        say $fh "\@" . $exon_id . "_" . $crispr_id;
        say $fh $seq;
    }
}

sub write_yaml {
    my ( $pairs, $exon_id ) = @_;

    die "yaml output file location not defined!" unless defined $yaml_file;

    YAML::Any::DumpFile( 
        $yaml_file, 
        {
            $exon_id => _get_unique_crisprs( $pairs )
        }
    );
}

#used to extract only unique crisprs from all the pairs
sub _get_unique_crisprs {
    my ( $pairs ) = @_;

    #create and return a hash of crispr ids pointing to crispr hashrefs
    my %unique;
    for my $pair ( @{ $pairs } ) {
        #first_crisprs have ids ending in A, second_crisprs end in B
        #some crisprs are in both lists (intentionally) -- for those we want
        #both to be in the hash, which is why they have different ids.
        for my $crispr ( qw( first_crispr second_crispr ) ) {
            my $id = $pair->{$crispr}{id};

            #add it if we don't already have it in the hash
            if( ! defined $unique{$id} ) {
                $unique{$id} = $pair->{$crispr};
            }
        }
    }

    return \%unique;
}

sub write_csv {
    my ( $pairs, $exon_id ) = @_;

    #chuck a csv onto stdout 
    #csv header
    say "Exon_ID, First_Crispr, Spacer, Spacer_Length, Second_Crispr";
    
    say $_ for map { join ",", $exon_id, 
                               $_->{first_crispr}->{seq}, 
                               $_->{spacer}, 
                               $_->{spacer_length}, 
                               $_->{second_crispr}->{seq} } @{ $pairs };
}

sub get_matches_fixed {
    my ( $slice ) = @_;

    my $seq = $slice->seq;

    my ( @pam_left, @pam_right );

    my $id = 1;

    while ( $seq =~ /(CC\S{21}|\S{21}GG)/g ) {
        #
        # NEED to check this is what we actually expect (the pos)
        #
        my $crispr_seq = $1;

        #
        #start/end are relative -- need to add $slice->start
        #

        my $data = {
            seq   => $crispr_seq,
            start => $slice->start + ((pos $seq) - length( $crispr_seq )), #is this right? i hope so
            end   => $slice->start + pos $seq,
            chr   => $slice->seq_region_name, #the same for all of them but we need it for the db
            id    => $id++,
            strand => $slice->strand,
        };

        my $type;

        #determine the direction, name appropriately and add to the correct list.
        if ( $crispr_seq =~ /^CC.*GG$/ ) {
            #its left AND right. what a joke
            $type = "both";
            my $right_data = { %$data }; #shallow copy data so they can be edited separately 
            $data->{id}       .= "A"; #this is left_data
            $right_data->{id} .= "B";

            push @pam_left, $data;
            push @pam_right, $right_data;
            
        }
        elsif ( $crispr_seq =~ /^CC/ ) {
            $type = "left";
            $data->{id} .= "A";

            push @pam_left, $data;
        }
        elsif ( $crispr_seq =~ /GG$/ ) {
            $type = "right";
            $data->{id} .= "B";

            push @pam_right, $data;
        }
        else {
            die "Crispr doesn't have a valid PAM site.";
        }

        say STDERR "Found $type crispr: $crispr_seq\@" . $slice->seq_region_name . ":"
                                                       . ($slice->start+$data->{start}) . "-" 
                                                       . ($slice->start+$data->{end}) . ":"
                                                       . $slice->strand;

        #go back to just after the pam site so we can find overlapping crisprs
        pos($seq) -= length( $crispr_seq ) - 2;
    }

    say STDERR "Found " . scalar(@pam_left) . " left crisprs and " . scalar(@pam_right) . " right crisprs";

    return get_pairs( \@pam_left, \@pam_right );
}

##
#
# COMPARE OUTPUT OF THIS WITH THE HPRT HUMAN STUFF I MADE THE OTHER DAY.
# IS THIS A SUPERSET. it had better be
#
##

sub get_pairs {
    my ( $pam_left, $pam_right ) = @_;

    my @pairs;

    #compare every left/right possibility, and see if they're a valid pair

    for my $l ( @{ $pam_left } ) {
        for my $r ( @{ $pam_right } ) {
            #
            # NEED TO SET MIN/MAX SPACER
            #
            my $distance = $l->{ end } - $r->{ start };
            if ( $distance <= $MAX_SPACER && $distance >= $MIN_SPACER ) {
                push @pairs, {
                    first_crispr  => $l,
                    spacer        => '-', #we can't easily get the spacer seq
                    spacer_length => $distance,
                    second_crispr => $r,
                }
            }

            #the lists are in order so if one is already smaller than the min spacer
            #all the remaining ones will be for this crispr
            #uncomment this when the rest is tested, check output is identical
            #last if $distance < $MIN_SPACER;
        }
    }

    say STDERR "Found " . scalar(@pairs) . " pairs.";

    return @pairs;
}

#this is the old method, still here while i check values
sub get_matches {
    my ( $seq ) = @_;

    my @matches;

    #while ( $seq =~ /CC(\S{21})(\S{0,30}?)(\S{21})GG/g ) {
    while ( $seq =~ /(CC\S{32,72}?GG)/g ) {
        #print "Match CC=>GG: \n CC$1 $2 $3GG (" . length($2) . ")\n";

        #extract the first and second (potentially overlapping) crisprs
        my ( $first_crispr, $second_crispr ) = ( substr( $1, 0, 23 ), substr( $1, -23 ) );

        #could just do - 46...
        my $spacer_length = length( $1 ) - length( $first_crispr ) - length( $second_crispr );
        #if the spacer length is 0 or less then its empty, but the substr would take sequence
        my $spacer = ( $spacer_length > 0 ) ? substr( $1, 23, -23 ) : "";

        push @matches, {
            first_crispr  => { seq => $first_crispr },
            spacer        => $spacer,
            spacer_length => $spacer_length,
            second_crispr => { seq => $second_crispr },
        };

        pos($seq) -= length($1) - 2;
    }

    #if ( @matches > 300 ) {
    #    print "Found " . scalar @matches . "\n";
    #}

    return @matches;
}

1;

__END__

=head1 NAME

find_paired_crisprs.pl - find paired crisprs for specific exons

=head1 SYNOPSIS

find_paired_crisprs.pl [options]

    --species     mouse or human
    --exon-ids    a list of ensembl exon ids
    --expand-seq  take an additional 200bp from the 5' and and an additional 
                  100bp from the 3' end. [optional, default true]
    --min-spacer  the minimum allowed spacer distance between pairs [optional, default -10]
    --max-spacer  the minimum allowed spacer distance between pairs [optional, default 50]
    --fq-file     where to output the fq data [optional]
    --yaml-file   where to output the yaml data [optional]
    --help        show this dialog

Example usage:

find_paired_crisprs.pl --species mouse --exon-ids ENSMUSE00001107660 ENSMUSE00001107880 --fq-file /nfs/users/nfs_a/ah19/work/crisprs/hprt_crisprs.fq --yaml-file /nfs/users/nfs_a/ah19/work/crisprs/hprt_crisprs.yaml
find_paired_crisprs.pl --species mouse --exon-ids ENSMUSE00001107660 ENSMUSE00001107880 --no-expand-seq --min-spacer 0 --max-spacer 30 > output.csv

=head1 DESCRIPTION

Find all possible crispr sites within exon sequences, and any possible pairs within them. 
The yaml file is created with the intention of being given to lims2-task to load into the database.
The fq file is intended to be handed to bwa.
It will also spit a csv out to stdout 

See paired_crisprs.sh for how I run this bad boy

=head AUTHOR

Alex Hodgkins

=cut