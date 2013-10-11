#!/usr/bin/perl

#usage:
#windowBed -a ar_with_n_valid_21-ARA.bed -b ar_with_n_valid_21-ARB.bed -w 999999 | awk '!($2==$8 && $1==$7)' | perl ~/work/paired_crisprs/valid_pairs.pl ~/work/paired_crisprs/ar_crisprs.fq > ~/work/paired_crisprs/ar-a_ar-b.html

use strict;
use warnings;

use Bio::Perl;
use Try::Tiny;
use Template;
use Data::Dumper;

my $tt = Template->new(PRE_CHOMP  => 1) || die Template->error();

#including pam mismatch
my $MAX_EDIT_DISTANCE = 6;
my $MAX_PAIRS = 11;

my %pams = (
    "ARA"   => { pam => "CCA", revpam => "TGG" },
    "ARB"   => { pam => "CCT", revpam => "AGG" },
    "RAG11" => { pam => "CCA", revpam => "TGG" },
    "RAG12" => { pam => "CCG", revpam => "CGG" },
    "RAG1C" => { pam => "CCA", revpam => "TGG" },
    "RAG1E" => { pam => "CCA", revpam => "TGG" },
    "RAG1F" => { pam => "CCA", revpam => "TGG" },
    "RAG1G" => { pam => "CCT", revpam => "AGG" },
);

sub hamming_distance {
    #use string xor to get the number of mismatches between the two strings.
    #the xor returns a string with the binary digits of each char xor'd, 
    #which will be an ascii char between 001 and 255. tr returns the number of characters replaced. 
    die "Strings passed to hamming distance differ" if length($_[0]) != length($_[1]);
    return (uc($_[0]) ^ uc($_[1])) =~ tr/\001-\255//;
}

die "Usage: cat ar-a_with_seqs.bed | valid_pairs.pl <crisprs.fq>" unless defined $ARGV[0];

#first file should be crisprs.fq, second the bed file
my $fq_file = shift;
open(my $fq_fh, "<", $fq_file) or die "Couldn't open $fq_file: $!";

my ( %crisprs, $current );
while ( my $line = <$fq_fh> ) {
    chomp $line;

    if ( $line =~ /^@(.*)/ ) {
        $current = $1;
    }
    else {
        die "Error: there must only be one sequence per entry" if defined $crisprs{ $current };
        #die "Crisprs should be reverse complemented!" unless $line =~ /^CC/;
        #pam must be at the end or the bwa wasn't the current one.
        die "Crisprs shouldn't be reverse complemented!" unless $line =~ /GG$/; 

        #add fwd and reverse so we can compute hamming distances later.
        #we call it rev as we assume all crispr sequences have been reverse complemented.
        #cant automatically determine which is which cause CC-GG is a valid crispr
        $crisprs{$current} = { rev => revcom( $line )->seq, fwd => $line };
    }
}

#my ( $pam_1, $pam_2 ) = ( "CCA", "AGG" );

#also need to try the reverse way. switch them around:
#my ( $revpam_1, $revpam_2 ) = ( revcom( $pam_2 )->seq, revcom( $pam_1 )->seq );

my %duplicates;
my $shortest = undef;

my %pairs;
my %names;
my $num_valid_pairs = 0;

while ( my $line = <STDIN> ) {
    chomp $line;
    #chr, start, end, name, unknown, strand x2
    my ( @data ) = split /\s+/, $line;

    die "Chromosomes don't match: $line" if $data[0] ne $data[6];

    #hashes would be better
    my ( @first_crispr, @second_crispr );

    #which crispr is earlier in the chromosome? that is the first crispr
    #if they're the regular way round, everything is good. we just search as normal,
    if ( $data[1] < $data[7] ) {
        @first_crispr  = @data[ 0 .. 5 ];
        @second_crispr = @data[ 6 .. 11 ];
    }
    else {
        #we need them the other way as they're reversed. we also need to switch them around,       

        @first_crispr  = @data[ 6 .. 11 ];
        @second_crispr =  @data[ 0 .. 5 ];
    }

    my $split = qr/(.+)-([CGATcgat]+)/;

    # if they are overlapping we need to check which side the crisprs are on or we'll get them the wrong way maybe
    if ( $second_crispr[1] - $first_crispr[2] < 0 ) {
        if ( get_direction( $first_crispr[3] =~ $split ) eq "fwd" && get_direction( $second_crispr[3] =~ $split ) eq "rev" ) {
            my @tmp = @first_crispr;
            @first_crispr = @second_crispr;
            @second_crispr = @tmp;
        }
    }

    my ( $lname, $lseq ) = $first_crispr[3] =~ /(.+)-([CGATcgat]+)/;
    my ( $rname, $rseq ) = $second_crispr[3] =~ /(.+)-([CGATcgat]+)/;

    #for labelling the output
    $names{$lname}++;
    $names{$rname}++;

    #die "'$lname' '$lseq' '$rname' '$rseq'";

    die "Unknown name: $lname" unless exists $crisprs{ $lname };
    die "Unknown name: $rname" unless exists $crisprs{ $rname };

    #get pam variables.

    #pam and revpam should have their names swapped. whatever, its confusing anyway. re-write
    #my ( $pam_1, $pam_2, $revpam_1, $revpam_2 ) = ( $pams{$lname}->{pam}, $pams{$rname}->{revpam}, $pams{$rname}->{pam}, $pams{$lname}->{revpam} );

    #see if the pam sites are end to end and if so is the sequence in the correct orientation.
    #we need this for sites that are CC-GG or GG-CC as it would appear to be end to end,
    #but the actual match

    #check we have a valid orientation

    my ( $tail_to_tail, $head_to_head, $type ) = ( 0, 0, "" );
    my ( $site_a, $site_b ); #these are the original crisprs. easier here

    if ( check_sequence( $lseq, $lname, "left" ) && check_sequence( $rseq, $rname, "right") ) {
    #if (   ( $lseq =~ /^CC/ and hamming_distance($lseq, $crisprs{$lname}->{rev}) <= $MAX_EDIT_DISTANCE )
    #    && ( $rseq =~ /GG$/ and hamming_distance($rseq, $crisprs{$rname}->{fwd}) <= $MAX_EDIT_DISTANCE ) ) {

        $tail_to_tail = 1;
        $type = "Tail to tail";

        $site_a = $crisprs{$lname}->{rev};
        $site_b = $crisprs{$rname}->{fwd};

        #colour mismatches. last param is pam_left
        $lseq = colour_mismatches( $lseq, $site_a, 1 );
        $rseq = colour_mismatches( $rseq, $site_b, 0 );

    }
    elsif ( check_sequence( $lseq, $lname, "right" ) && check_sequence( $rseq, $rname, "left" ) ) {
    #elsif (   ( $lseq =~ /GG$/ and hamming_distance($lseq, $crisprs{$rname}->{fwd}) <= $MAX_EDIT_DISTANCE )
    #       && ( $rseq =~ /^CC/ and hamming_distance($rseq, $crisprs{$rname}->{rev}) <= $MAX_EDIT_DISTANCE ) ) {

        $head_to_head = 1;
        $type = "Head to head";

        $site_a = $crisprs{$lname}->{fwd};
        $site_b = $crisprs{$rname}->{rev};

        #colour mismatches
        $lseq = colour_mismatches( $lseq, $site_a, 0 );
        $rseq = colour_mismatches( $rseq, $site_b, 1 );
    }

    #does the left sequence start CC AND have less than 6 mismatches with the crispr beginning CC?
    #does the second sequence end GG and have less than 6 mismatches with the crispr seq ending GG?
    if ( $tail_to_tail or $head_to_head ) {

        #store all matches by their start position and chromosome, which allows us to 
        #identify duplicates.
        my $dupe = 0;
        if ( exists $duplicates{$second_crispr[0]."_".$second_crispr[1]} ) {
            if ( $duplicates{$second_crispr[0]."_".$second_crispr[1]}->{$first_crispr[0]."_".$first_crispr[1]} ) {
                $dupe = 1;
            }
        }

        if ( $dupe ) {
            #print STDERR "$line is a duplicate!\n";
        }
        else {
            #print STDERR "$line is a valid crispr.\n";
            $num_valid_pairs++;
            #print "$line\n";

            my $distance = ($second_crispr[1] > $first_crispr[2]) ? ($second_crispr[1] - $first_crispr[2]) : ($first_crispr[1] - $second_crispr[2]);

            push @{ $pairs{$type} }, {
                #first_crispr  => \@first_crispr,
                #second_crispr => \@second_crispr,
                first_locus   => $first_crispr[0] . ":" . ($first_crispr[1]+1) . "-" . $first_crispr[2],
                second_locus  => $second_crispr[0] . ":" . ($second_crispr[1]+1) . "-" . $second_crispr[2],
                first_name    => $lname,
                second_name   => $rname,
                site_a        => $site_a,
                site_b        => $site_b,
                first_seq     => $lseq,
                second_seq    => $rseq,
                distance      => $distance,
                #type          => $type,
            };
            
            
            #set shortest if its not already set or distance is less. 
            $shortest = $distance if ! defined $shortest or $distance < $shortest;
        }

        $duplicates{$second_crispr[0]."_".$second_crispr[1]}->{$first_crispr[0]."_".$first_crispr[1]}++;
    }
    #else {
    #    die "$line didn't match! (^$pam_1 and $pam_2\$)";
    #}

}

sub get_direction {
    my ( $name, $seq ) = @_;

    my $fwd_distance = hamming_distance($seq, $crisprs{$name}->{fwd});
    my $rev_distance = hamming_distance($seq, $crisprs{$name}->{rev});

    die "Identical hamming distances found for both orientations of $seq" if $fwd_distance == $rev_distance; 

    return ( $fwd_distance < $rev_distance ) ? "fwd" : "rev";
}

sub check_sequence {
    my ( $seq, $name, $direction ) = @_;

    #print STDERR "running check seq with $seq $name $direction\n";

    my $fwd_distance = hamming_distance($seq, $crisprs{$name}->{fwd});
    my $rev_distance = hamming_distance($seq, $crisprs{$name}->{rev});

    #my $crispr_seq = ( $fwd_distance  < $rev_distance ) ? $crisprs{$name}->{fwd} : $crisprs{$name}->{rev};

    #die if $smallest > $MAX_EDIT_DISTANCE;

    #direction should actually be called location or something. its just the position
    if ( $fwd_distance < $rev_distance ) {
        return 1 if $seq =~ /GG$/ and $direction eq "right";
    }
    elsif ( $rev_distance < $fwd_distance ) {
        #valid pam left
        return 1 if $seq =~ /^CC/ and $direction eq "left";
    }
    else {
        #if they're equal i have no idea what to do.
        die "Identical hamming distances found for both orientations of $seq";
    }

    #its not valid for the direction specified
    return 0;
}

die "More than 2 crisprs being compared?? " . join(", ", keys %names) if scalar keys %names > 2;  

#make html for each type
for my $type ( keys %pairs ) {
    my @rows = sort { $a->{distance} <=> $b->{distance} } @{ $pairs{$type} };

    my $html = get_html( {
        rows        => [ splice( @rows, 0, $MAX_PAIRS ) ],
        first_name  => $rows[0]->{first_name}, #this is fine, right?? not really. should check same on all
        second_name => $rows[0]->{second_name},
        site_a      => $rows[0]->{site_a}, 
        site_b      => $rows[0]->{site_b},
        type        => $type
    } );

    #print join("\t", @{ $pair->{first_crispr} }, @{ $pair->{second_crispr} }, $pair->{distance} ) . "\n";
    #should be out of the loop. whatever
    print <<"END";
<style type="text/css">
table { font-family: 'Courier New', monospace; border-collapse:collapse; border:1px solid black; }
table td { border:1px solid black; }
</style>
END
    print $html, "\n";
}

sub get_html {
    my ( $data ) = @_;

    my $html = <<"END";
<h1>[% first_name %] spacer [% second_name %] ([% type %])</h1>
<table style="font-family:'Courier New',monospace;">
    <thead>
        <tr>
            <td>Site A</td>
            <td>Locus</td>
            <td>Spacer</td>
            <td>Locus</td>
            <td>Site B</td>
        </tr>
    </thead>
    <tbody>
        [% FOR c IN rows %]
        <tr>
            <td>[% c.site_a %]<br/>[% c.first_seq %] ([% c.first_name %])</td>
            <td><a href="http://www.ensembl.org/Mus_musculus/psychic?q=[% c.first_locus %]">[% c.first_locus %]</a></td>
            <td>[% c.distance %]</td>
            <td><a href="http://www.ensembl.org/Mus_musculus/psychic?q=[% c.second_locus %]">[% c.second_locus %]</a></td>
            <td>[% c.site_b %]<br/>[% c.second_seq %] ([% c.second_name %])</td>
        </tr>
        [% END %]
    </tbody>
</table>
END

    my $processed_html;
    $tt->process( \$html, $data, \$processed_html );

    return $processed_html;
}

sub colour_mismatches {
    my ( $seq, $ref, $pam_left ) = @_;

    #
    # gross. gross
    #

    my ( $new, $start, $end );

    if ( $pam_left ) {
        $start = 3;
        $end = length($seq);
        $new = "<span style='font-weight:bold;'>" . substr($seq, 0, 3) . "</span>"; #start with pam
    }
    else {
        $start = 0;
        $end = length($seq) - 3;
    }
    
    #need to know if pam is left or right side so we can start/end appropriately.
    #can't tell from just the two seqs

    for ( my $i=$start; $i<$end; $i++ ) {
        if ( substr($seq, $i, 1) ne substr($ref, $i, 1) ) {
            $new .= "<span style='color:#FF0000;'>" . substr($seq, $i, 1) . "</span>"
        }
        else {
            $new .= substr($seq, $i, 1);
        }
    }

    #we need to add the pam as it wasn't looped
    $new .= "<span style='font-weight:bold;'>" . substr($seq, $end, 3) . "</span>" if ! $pam_left;

    return $new;
}

#
# add head to head config as valid
#    $lseq =~ /GG$/ and hamming_distance($lseq, $crisprs{$rname}->{fwd}) <= 6 )
# && $rseq =~ /^CC/ and hamming_distance($rseq, $crisprs{$rname}->{rev}) <= 6 )
#
# instead of shortest store 10 best.
# output to html table with coloured differences
#

#stop concat error
$shortest = "None found." unless defined $shortest;

print STDERR "Closest crispr: $shortest\n";
print STDERR "Total valid pairs: $num_valid_pairs\n";
