#!/usr/bin/perl

use strict;
use warnings;

use feature qw( say );
use Bio::Perl;

use Tie::File;
use Fcntl 'O_RDONLY';

my @chromosome;
my $genomes_folder = "/lustre/scratch110/blastdb/Users/team87/Mouse/GRCm38/";

die "You must strip the pam." if length($ARGV[0]) != 20;

my ( $crispr, $rev_crispr ) = ( $ARGV[0], revcom( $ARGV[0] )->seq ); 
my $crispr_len = 20;

#just 1 for now to test
tie @chromosome, 'Tie::File', $genomes_folder . "1.fasta", mode => O_RDONLY, autochomp => 1 
    or die "Unable to open file: $!";
my $chr = "1"; # also need to set this per file

my ( @matches, $buf );

my $num_lines = scalar @chromosome;
my $line_start = 1;

my $all_n = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";

for my $line_num ( 0 .. $num_lines-1 ) {
    my $line = $chromosome[$line_num];
    $line_start += length($line); #so we can easily get the start location

    # hamming_distance($line, $all_n) or something

    next if length($line) < $crispr_len || $line =~ /^>/;

    #append 20 bp from the next line if there is one
    if ( $line_num != $num_lines-1 ) {
        #add 19 bp from the next line so we can check up to the very last bp of this line
        #should change this to take another 3bp so we can do pams
        $line .= substr $chromosome[$line_num+1], 0, $crispr_len-1;
    }

    #check every single bp on this line
    for my $index ( 0 .. ( length($line) - $crispr_len ) ) {
        #
        #check this is right
        #
        my $seq = substr $line, $index, $crispr_len;

        #say "$line(${line_num}): Index:$index" . ", size:" . ($index + $crispr_len) . ", len:" . length($seq);

        my $strand;
        #add bed file entry to @matches
        if ( hamming_distance( $seq, $crispr ) <= 5 ) {
            $strand = '+';
        }
        elsif ( hamming_distance( $seq, $rev_crispr ) <= 5 ) {
            $strand = '-';
        }
        else {
            next; #it isn't valid
        }

        #we could do with taking the next 3 base pairs too (annoying if we're at the end of a line,)
        #or maybe not? just take 22bp from next line instead of 19. and -3 from substr loop

        my $start = $line_start + $index;

        #should get a real ID
        push @matches, join( "\t", $chr, $start, $start+$crispr_len+3, "OT-$seq", 0, $strand );
        say join( "\t", $chr, $start, $start+$crispr_len+3, "OT-$seq", 0, $strand );
    }
}

say STDERR "found " . scalar( @matches ) . " off targets.";
say $_ for @matches;

untie @chromosome;

sub hamming_distance {
    #use string xor to get the number of mismatches between the two strings.
    #the xor returns a string with the binary digits of each char xor'd, 
    #which will be an ascii char between 001 and 255. tr returns the number of characters replaced. 
    die "Strings passed to hamming distance differ: " . $_[0] . ", " . $_[1]
        if length($_[0]) != length($_[1]);
    return (uc($_[0]) ^ uc($_[1])) =~ tr/\001-\255//;
}