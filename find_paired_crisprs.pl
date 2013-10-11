#!/usr/bin/perl

use strict;
use warnings;

use feature qw( say );

use LIMS2::Util::EnsEMBL;

die "Usage: find_paired_crisprs.pl <species> <exon_ids>" unless @ARGV > 1;

#everything else should be ensembl ids
my $species = shift @ARGV;

#are these correct?
my ( $MIN_SPACER, $MAX_SPACER ) = ( -10, 50 );

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

#csv header
say "Exon_ID, First_Crispr, Spacer, Spacer_Length, Second_Crispr";

#
# exons with gibson designs:
# ENSE00003596810 ENSE00003663376 ENSE00003606976 ENSE00003682156 ENSE00003611316 ENSE00001427230 ENSE00001083354 ENSE00000683019 ENSE00000849850 ENSE00003701380 ENSE00003596529 
#

for my $exon_id ( @ARGV ) {
    if ( $exon_id !~ $exon_re ) {
        say STDERR "$exon_id is not a valid exon for this species!";
        next;
    }

    say STDERR "Finding crisprs in $exon_id";

    #get the exon and 200bp surrounding
    my $exon_slice = $e->slice_adaptor->fetch_by_exon_stable_id( $exon_id, 200 );

    my ( $seq ) = $exon_slice->seq;
    #we need the gene so we can get the strand and subsequently chop 100bp off the seq
    #because we want an asymmetrical slice
    # 5' -----[ EXON ]----- 3'
    #      200        100
    my $strand = $e->gene_adaptor->fetch_by_exon_stable_id( $exon_id )->strand;

    #strip away the excess 100bp we took because ensembl doesnt let us take asymmetrical slices
    if ( $strand == 1 ) {
        $seq = substr $seq, 0, -100;
    }
    elsif ( $strand == -1 ) {
        $seq = substr $seq, 100;
    }
    else {
        die "Unexpected strand for gene associated with $exon_id";
    }

    #say $exon_slice->seq, " (" . length($exon_slice->seq) . ")";
    #my @matches = get_matches( $exon_slice->seq );
    my @matches = get_matches_fixed( $exon_slice );
    #say "Found " . scalar( @matches ) . " paired crisprs:";
    
    say $_ for map { join ",", $exon_id, 
                               $_->{first_crispr}, 
                               $_->{spacer}, 
                               $_->{spacer_length}, 
                               $_->{second_crispr} } @matches;

}

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
            first_crispr  => $first_crispr,
            spacer        => $spacer,
            spacer_length => $spacer_length,
            second_crispr => $second_crispr,
        };

        pos($seq) -= length($1) - 2;
    }

    #if ( @matches > 300 ) {
    #    print "Found " . scalar @matches . "\n";
    #}

    return @matches;
}

sub get_matches_fixed {
    my ( $slice ) = @_;

    my $seq = $slice->seq;

    my ( @pam_left, @pam_right );

    while ( $seq =~ /(CC\S{21}|\S{21}GG)/g ) {
        #
        # NEED to check this is what we actually expect (the pos)
        #
        my $crispr_seq = $1;

        my $data = {
            seq   => $crispr_seq,
            start => (pos $seq) - length( $crispr_seq ), #is this right? i hope so
            end   => pos $seq, 
        };

        my $type;

        if ( $crispr_seq =~ /^CC.*GG$/ ) {
            #its left and right. what a joke
            push @pam_left, $data;
            push @pam_right, $data;
            $type = "both";
        }
        elsif ( $crispr_seq =~ /^CC/ ) {
            push @pam_left, $data;
            $type = "left";
        }
        elsif ( $crispr_seq =~ /GG$/ ) {
            push @pam_right, $data;
            $type = "right";
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
                    first_crispr  => $l->{ seq },
                    spacer        => '-', #we can't easily get the spacer seq
                    spacer_length => $distance,
                    second_crispr => $r->{ seq },
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