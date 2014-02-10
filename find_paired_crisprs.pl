#!/usr/bin/perl

use strict;
use warnings;

use feature qw( say );

use Getopt::Long;
use Pod::Usage;
use LIMS2::Util::EnsEMBL;
#use LIMS2::Model;
use Bio::Perl qw( revcom );
use YAML::Any qw( DumpFile );

use Data::Dumper;

#should these all be caps? maybe. its really ugly though
my ( $species, $crispr_fq_file, $crispr_yaml_file, $pair_txt_file ); 
my ( $pair_yaml_file, @exon_ids, @regions, @pair_ids, %specific_sites, $pair_gff_file );
my ( $expand_seq, $MIN_SPACER, $MAX_SPACER ) = ( 1, -10, 30 ); #set default values
GetOptions(
    "help"               => sub { pod2usage( 1 ) },
    "man"                => sub { pod2usage( 2 ) },
    "species=s"          => sub { my ( $name, $val ) = @_; $species = lc $val; },
    "exon-ids=s{,}"      => \@exon_ids,
    "regions=s{,}"       => \@regions,
    "pair-ids=s{,}"    => \@pair_ids,
    "specific-sites=s{,}" => sub { my ( $name, $val ) = @_; $specific_sites{$val} = 1; },
    "fq-file=s"          => \$crispr_fq_file,
    "crispr-yaml-file=s" => \$crispr_yaml_file,
    "pair-yaml-file=s"   => \$pair_yaml_file,
    "pair-text-file=s"   => \$pair_txt_file,
    "pair-gff-file=s"    => \$pair_gff_file,
    "expand-seq!"        => \$expand_seq,
    "min-spacer=i"       => \$MIN_SPACER,
    "max-spacer=i"       => \$MAX_SPACER,
) or pod2usage( 2 );

die pod2usage( 2 ) unless $species and (@exon_ids or @regions or @pair_ids);

if ( scalar keys %specific_sites ) {
    say "Specific sites are: " . join ",", keys %specific_sites;
}

#select the appropriate exon regex to validate any exon ids we get
my %exon_regexes = (
    'mouse' => qr/ENSMUSE/,
    'human' => qr/ENSE/,
);

die "Unknown species: $species" unless defined $exon_regexes{ $species };
my $exon_re = $exon_regexes{ $species };
my $region_re = qr/(\d+|X|Y):(\d+)-(\d+)/;

my $e = LIMS2::Util::EnsEMBL->new( species => $species );
#my $l = LIMS2::Model->new( user => "lims2" );
my %pairs_by_exon;

#why arent these functions. this whole file is a mess now and needs to be modularised 

for my $region ( @regions ) {
    my ( $chr, $start, $end ) = $region =~ $region_re;

    unless ( defined $chr and defined $start and defined $end ) {
        say STDERR "$region is not a valid region!";
        next;
    }

    my $slice = $e->slice_adaptor->fetch_by_region("chromosome", $chr, $start, $end);
    unless ( $slice ) {
        say STDERR "Couldn't get slice, skipping $region";
        next;
    }

    say STDERR "Region length: " . $slice->length;

    my $base = ( $species eq 'mouse' ) ? 'ENSMUSE9999' : 'ENSE9999';
    my $int_chr = ( $chr eq 'X' ) ? 22 : ( $chr eq 'Y' ) ? 23 : $chr; #change X/Y to 22/23
    my $fake_exon_id = $base . $int_chr . $start . $end;

    say STDERR "Fake exon id for $region: $fake_exon_id";

    my $matches = get_matches( $slice );
    say "Found " . scalar( @$matches ) . " paired crisprs:";

    $pairs_by_exon{ $fake_exon_id } = $matches;
}

for my $exon_id ( @exon_ids ) {
    if ( $exon_id !~ $exon_re ) {
        say STDERR "$exon_id is not a valid exon for this species!";
        next;
    }

    say STDERR "Finding crisprs in $exon_id";

    #get a slice of just the exon
    my $exon_slice = $e->slice_adaptor->fetch_by_exon_stable_id( $exon_id );
    say STDERR "Exon length: " .$exon_slice->length;

    #enable this if its a huge exon like AR
    #if ( $exon_slice->length > 400 ) {
    #    $exon_slice = $exon_slice->expand(0, -1100);
    #    say STDERR "New exon length: " .$exon_slice->length;
    #}

    #expand the slice if expand_seq is set (it is by default)
    if ( $expand_seq ) {        
        #we need the gene so we can get the strand and take an asymmetrical slice:
        # 5' -----[ EXON ]----- 3'
        #      200        100
        my $gene = $e->gene_adaptor->fetch_by_exon_stable_id( $exon_id );
        say STDERR "Gene identified as " . $gene->external_name;

        #expand the slice considering the strand.
        if ( $gene->strand == 1 ) {
            $exon_slice = $exon_slice->expand( 200, 100 );
        }
        elsif ( $gene->strand == -1 ) {
            $exon_slice = $exon_slice->expand( 100, 200 );
        }
        else {
            die "Unexpected strand for gene associated with $exon_id";
        }
    }

    die "Couldn't get slice" unless $exon_slice;

    #say $exon_slice->seq, " (" . length($exon_slice->seq) . ")";
    #my @matches = get_matches( $exon_slice->seq );

    #arrayref of pairs
    my $matches = get_matches( $exon_slice );
    #say "Found " . scalar( @$matches ) . " paired crisprs:";

    #add to the global pairs hash so we can output all data
    $pairs_by_exon{ $exon_id } = $matches;
}

if ( @pair_ids ) {
    #my @pairs = $l->schema->resultset( "CrisprPair" )->search(
    #    { 'me.id' => { -IN => \@pair_ids } },
    #    { prefetch => [ "left_crispr", "right_crispr" ] }
    #);
    my @pairs;

    #emulate the format used by everything else
    my @matches;
    for my $pair ( @pairs ) {
        push @matches, {
            first_crispr  => { 
                id    => $pair->left_crispr_id."A", 
                seq   => $pair->left_crispr->seq, 
                locus => $pair->left_crispr->loci->first->as_hash, 
            }, 
            second_crispr => { 
                id    => $pair->right_crispr_id."B", 
                seq   => $pair->right_crispr->seq,
                locus => $pair->right_crispr->loci->first->as_hash, 
            }
        };
    }

    #ENS needed to conform to regex
    $pairs_by_exon{ ENSExistingCrisprs } = \@matches;
}

#this should be a jump table or something
write_fq( \%pairs_by_exon ) if defined $crispr_fq_file;
write_crispr_yaml( \%pairs_by_exon ) if defined $crispr_yaml_file;
write_pair_yaml( \%pairs_by_exon ) if defined $pair_yaml_file;
write_pair_gff( \%pairs_by_exon ) if defined $pair_gff_file;

#this is how it worked before the fq/yaml option so recreate that
if ( ! defined $crispr_fq_file && ! defined $crispr_yaml_file ) {
    write_csv( \%pairs_by_exon );
}

sub write_fq {
    my ( $pairs_by_exon ) = @_;

    die "fq output file location not defined!" unless defined $crispr_fq_file;
    open( my $fh, ">", $crispr_fq_file ) || die "Couldn't open $crispr_fq_file for writing.";

    while ( my ( $exon_id, $pairs ) = each %{ $pairs_by_exon } ) {
        #print out all unique crisprs in fq format. if it ends in A it's a CC crispr and needs to be reverse
        #complemented for bwa
        my $unique_crisprs = _get_unique_crisprs( $pairs );
        #sort by integer part of id, ignoring A/B
        for my $crispr_id ( sort { ($a =~ /(\d+)[AB]$/)[0] <=> ($b =~ /(\d+)[AB]$/)[0] } keys %{ $unique_crisprs } ) {
            my $seq = $unique_crisprs->{$crispr_id}{seq};
            if ( $crispr_id =~ /A$/ ) {
                $seq = revcom( $seq )->seq;    
            }

            say $fh "\@" . $exon_id . "_" . $crispr_id;
            say $fh $seq;
        }
    }

    say STDERR "Crispr fq output written to $crispr_fq_file";
}

sub write_crispr_yaml {
    my ( $pairs_by_exon ) = @_;

    die "yaml output file location not defined!" unless defined $crispr_yaml_file;

    my %unique_crisprs;
    while ( my ( $exon_id, $pairs ) = each %{ $pairs_by_exon } ) {
        $unique_crisprs{ $exon_id } = _get_unique_crisprs( $pairs );
    }

    YAML::Any::DumpFile( $crispr_yaml_file, \%unique_crisprs );

    say STDERR "Crispr yaml output written to $crispr_yaml_file";
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
                $pair->{$crispr}{pam_right} = ( $id =~ /B$/ ) || 0; #set pam_right to true/false
                $unique{$id} = $pair->{$crispr};
            }
        }
    }

    return \%unique;
}

sub write_csv {
    my ( $pairs_by_exon ) = @_;

    #chuck a csv onto stdout 
    #csv header
    say "Exon_ID, First_Crispr, Spacer, Spacer_Length, Second_Crispr";

    while ( my ( $exon_id, $pairs ) = each %{ $pairs_by_exon } ) {
        say $_ for map { join ",", $exon_id, 
                                   $_->{first_crispr}->{seq}, 
                                   $_->{spacer}, 
                                   $_->{spacer_length}, 
                                   $_->{second_crispr}->{seq} } @{ $pairs };
    }
}

sub write_pair_yaml {
    my ( $pairs_by_exon ) = @_;

    die "pair yaml output file location not defined!" unless defined $pair_yaml_file;

    #we only want the pair information so we need to strip it out.
    #remember that pair_data is a hash whose keys point to an arrayref of hashrefs, e.g.:
    # ( ENSE0003596810 => [ {pair_hash_ref}, {pair_hash_ref}, ... ] )
    my %pair_data;
    while ( my ( $exon_id, $pairs ) = each %{ $pairs_by_exon } ) {
        #for every pair in this exon make a hashref of pair data.
        #store in an arrayref.
        my @db_pairs = map {
            { #each map iteration returns a hashref
                left_crispr  => $_->{first_crispr}{id},
                right_crispr => $_->{second_crispr}{id},
                spacer       => $_->{spacer_length},
            }
        } @{ $pairs };

        $pair_data{ $exon_id } = \@db_pairs;
    }

    YAML::Any::DumpFile( $pair_yaml_file, \%pair_data );

    say STDERR "Pair yaml output written to $pair_yaml_file";

    return;
}

sub write_pair_gff {
    my ( $pairs_by_exon ) = @_;

    die "pair yaml output file location not defined!" unless defined $pair_gff_file;

    open( my $fh, ">", $pair_gff_file ) || die "Couldn't open $pair_gff_file for writing: $!";

    while ( my ( $exon_id, $pairs ) = each %{ $pairs_by_exon } ) {
        my $id = 1;
        for my $pair ( @{ $pairs } ) {
            for my $crispr ( qw( first_crispr second_crispr ) ) {
                say $fh join "\t", 'chr'.$pair->{$crispr}{locus}{chr_name},
                                   'pair_finder',
                                   'pair',
                                   $pair->{$crispr}{locus}{chr_start},
                                   $pair->{$crispr}{locus}{chr_end},
                                   0,
                                   '.',
                                   '.',
                                   'group' . $id;
            }

            $id++;
        }
    }

    say STDERR "Pair gff output written to $pair_gff_file";

    return;
}

#we have it as text so it can be processed by paired_crisprs.sh easily
sub write_pair_txt {
    my ( $pairs_by_exon ) = @_;

    die "pair yaml output file location not defined!" unless defined $pair_txt_file;

    open(my $fh, "<", $pair_txt_file) || die "Couldn't open $pair_txt_file for reading.";

    #print out 1 pair per line separated by a space.
    while ( my ( $exon_id, $pairs ) = each %{ $pairs_by_exon } ) {
        for my $pair ( @{ $pairs } ) {
            say $fh "${exon_id}_" . $pair->{first_crispr}{id} . " "
                  . "${exon_id}_" . $pair->{second_crispr}{id};
        }
    }

    say STDERR "Pair txt output written to $pair_txt_file";
}

sub get_matches {
    my ( $slice ) = @_;

    my $seq = $slice->seq;

    my ( @pam_left, @pam_right );

    my $id = 1;

    say STDERR "Slice location:" . $slice->start . "-" . $slice->end . " (" . ($slice->length) . "bp)";

    while ( $seq =~ /(CC\S{21}|\S{21}GG)/g ) {
        my $crispr_seq = $1;

        if ( %specific_sites && ! defined $specific_sites{$crispr_seq} ) {
            say STDERR "Skipping $crispr_seq as it doesn't match specified sites.";
            pos($seq) -= length( $crispr_seq ) - 1;
            next;
        }

        my $data = {
            locus => {
                #chr_strand => $slice->strand, #this is specified below because it depends on the crispr orientation
                chr_start => $slice->start + ((pos $seq) - length( $crispr_seq )), #is this right? i hope so
                chr_end   => $slice->start + ((pos $seq) - 1), #need to subtract 1 or we get 24bp sequence back
                chr_name  => $slice->seq_region_name, #the same for all of them but we need it for the db
            },
            id        => $id++, #used to identify crisprs to add off targets later
            type      => 'Exonic', #for now -- some might be intronic
            seq       => $crispr_seq,
            #species   => $species,
            off_target_algorithm => 'bwa',
        };

        my $type;

        #determine the direction, name appropriately and add to the correct list.
        if ( $crispr_seq =~ /^CC.*GG$/ ) {
            #its left AND right. what a joke
            $type = "both ";
            #note that they still both point to the same locus, because why would you change that??
            my $right_data = { %$data }; #shallow copy data so they can be edited separately 
            $data->{id}       .= "A"; #pretend this is called left_data in this block
            $right_data->{id} .= "B";

            #the left crispr (beginning CC) is on the -ve, and the GG is on the +ve
            $data->{locus}{chr_strand} = -1;
            $right_data->{locus}{chr_strand} = 1;

            push @pam_left, $data;
            push @pam_right, $right_data;
        }
        elsif ( $crispr_seq =~ /^CC/ ) {
            $type = "left ";
            $data->{id} .= "A";

            #its CC so its on the global negative strand (all slices we have are on global +ve)
            $data->{locus}{chr_strand} = -1;

            push @pam_left, $data;
        }
        elsif ( $crispr_seq =~ /GG$/ ) {
            $type = "right";
            $data->{id} .= "B";

            #GG means it is on the global positive strand
            $data->{locus}{chr_strand} = 1;

            push @pam_right, $data;
        }
        else {
            die "Crispr doesn't have a valid PAM site.";
        }

        say STDERR "Found $type crispr " . $data->{id} . ": $crispr_seq\@" . $slice->seq_region_name . ":"
                                                       . $data->{locus}{chr_start} . "-" 
                                                       . $data->{locus}{chr_end} . ":"
                                                       . $slice->strand;

        #go back to just after the pam site so we can find overlapping crisprs
        pos($seq) -= length( $crispr_seq ) - 1;
    }

    say STDERR "Found " . scalar(@pam_left) . " left crisprs and " . scalar(@pam_right) . " right crisprs";

    return get_pairs( \@pam_left, \@pam_right );
}

sub get_pairs {
    my ( $pam_left, $pam_right ) = @_;

    my @pairs;

    #compare every left/right possibility, and see if they're a valid pair

    for my $l ( @{ $pam_left } ) {
        for my $r ( @{ $pam_right } ) {
            #
            # NEED TO SET MIN/MAX SPACER
            #
            my $distance = $r->{locus}{ chr_start } - $l->{locus}{ chr_end };
            if ( $distance <= $MAX_SPACER && $distance >= $MIN_SPACER ) {
                push @pairs, {
                    first_crispr  => $l,
                    spacer        => '-', #we can't easily get the spacer seq
                    spacer_length => $distance,
                    second_crispr => $r,
                }
            }

            #
            #TODO: 
            #   add else statement to keep crisprs that aren't in pairs
            #

            #the lists are in order so if one is already smaller than the min spacer
            #all the remaining ones will be for this crispr
            #uncomment this when the rest is tested, check output is identical
            #last if $distance < $MIN_SPACER;
        }
    }

    say STDERR "Found " . scalar( @pairs ) . " pairs.";

    return \@pairs;
}

#this is the old method, still here while i check values
sub get_matches_old {
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

    --species            mouse or human
    --exon-ids           a list of ensembl exon ids
    --expand-seq         take an additional 200bp from the 5' and and an additional 
                         100bp from the 3' end. [optional, default true]
    --min-spacer         the minimum allowed spacer distance between pairs [optional, default -10]
    --max-spacer         the minimum allowed spacer distance between pairs [optional, default 50]
    --fq-file            where to output the crispr fq data [optional]
    --crispr-yaml-file   where to output the crispr yaml data [optional]
    --pair-yaml-file     wgere to output the pair yaml data [optional]
    --help               show this dialog

Example usage:

find_paired_crisprs.pl --species mouse --exon-ids ENSMUSE00001107660 ENSMUSE00001107880 --fq-file /nfs/users/nfs_a/ah19/work/crisprs/hprt_crisprs.fq --crispr-yaml-file /nfs/users/nfs_a/ah19/work/crisprs/hprt_crisprs.yaml
find_paired_crisprs.pl --species mouse --exon-ids ENSMUSE00001107660 ENSMUSE00001107880 --no-expand-seq --min-spacer 0 --max-spacer 30 > output.csv

=head1 DESCRIPTION

Find all possible crispr sites within exon sequences, and any possible pairs within them. 
The crispr-yaml file is created with the intention of being given to lims2-task to load into the database.
The fq file is intended to be handed to bwa.
It will also spit a csv out to stdout if you don't ask for a fq or yaml file

See paired_crisprs.sh for how I run this bad boy

=head AUTHOR

Alex Hodgkins

=cut
