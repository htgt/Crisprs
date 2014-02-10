#!/usr/bin/perl

#usage:
#windowBed -a ar_with_n_valid_21-ARA.bed -b ar_with_n_valid_21-ARB.bed -w 999999 | awk '!($2==$8 && $1==$7)' | perl ~/work/paired_crisprs/valid_pairs.pl ~/work/paired_crisprs/ar_crisprs.fq > ~/work/paired_crisprs/ar-a_ar-b.html

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Bio::Perl;
use Try::Tiny;
use Data::Dumper;
use feature qw( say );
use YAML::Any qw( DumpFile LoadFile );

use List::MoreUtils;
use List::Util qw(sum);

my ( $crispr_fq_file, $crispr_yaml_file, @paired_files, $pair_yaml_file );
GetOptions(
    "help"                => sub { pod2usage( 1 ) },
    "man"                 => sub { pod2usage( 2 ) },
    "fq-file=s"           => \$crispr_fq_file,
    "pair-yaml-file=s"    => \$pair_yaml_file,
    "crispr-yaml-file=s"  => \$crispr_yaml_file,
    "paired-output=s{,}"  => \@paired_files,
) or pod2usage( 2 );

die pod2usage( 2 ) unless $pair_yaml_file and $crispr_fq_file and $crispr_yaml_file and @paired_files;

my %CRISPRS = process_fq_file( $crispr_fq_file );
my %OFF_TARGET_DATA;
my $CRISPR_NAME_RE = qr/^(ENS.*)_.*$/;

#we need this data to make sure the off target we're adding isn't the same as the crispr
my $CRISPR_DATA = LoadFile( $crispr_yaml_file );

for my $filename ( @paired_files ) {
    my ( %duplicates, %crispr_names );
    my ( %num_valid_pairs, $shortest );

    my @all;

    open( my $fh, "<", $filename ) || die "Couldn't open $filename";

    #say STDERR "Processing $filename";

    while ( my $line = <$fh> ) {
        chomp $line;

        #split the line into two hashrefs of bed data
        my ( $l_crispr, $r_crispr ) = process_windowbed_line( $line );

        next unless $l_crispr and $r_crispr;

        #have we already processed this pair but the other way?
        #windowbed gives you each pair twice.
        if ( exists $duplicates{ $r_crispr->{locus_str} }->{ $l_crispr->{locus_str} } ) {
            next;
        }
        else {
            #we've done this one now so add it to our duplicates hash
            $duplicates{ $r_crispr->{locus_str} }->{ $l_crispr->{locus_str} }++;
        }

        my $pair_key = $l_crispr->{name}."-".$r_crispr->{name};

        if ( ! defined $num_valid_pairs{$pair_key} ) {
            $num_valid_pairs{$pair_key} = 0; #so we can accurately report how many for each
        }

        #for labelling the output
        $crispr_names{ $l_crispr->{name} }++;
        $crispr_names{ $r_crispr->{name} }++;

        #check we have a valid orientation:
        #see if the pam sites are end to end and if so is the sequence in the correct orientation.
        #we need this for sites that are CC-GG or GG-CC as it would appear to be end to end,
        #but the actual match

        my ( $valid_pair, $type ) = ( 0, "" );

        #we have already validated the seqs so we just need to check if its CC-GG or GG-CC
        if ( $l_crispr->{direction} eq "rev" && $r_crispr->{direction} eq "fwd" ) {
            $valid_pair = 1;
            $type = "Tail to tail";

        }
        elsif ( $l_crispr->{direction} eq "fwd" && $r_crispr->{direction} eq "rev" ) {
            $valid_pair = 1;
            $type = "Head to head";
        }

        if ( $valid_pair ) {
            #the entry wasn't a duplicate so lets process it

            #find the entry from the crisprs.yaml to see if this is the original pair
            my $l_crispr_yaml = $CRISPR_DATA->{ $l_crispr->{exon_id} }{ $l_crispr->{crispr_id} } || die;
            my $r_crispr_yaml = $CRISPR_DATA->{ $r_crispr->{exon_id} }{ $r_crispr->{crispr_id} } || die;

            my $left_identical = ( $l_crispr->{chr} eq $l_crispr_yaml->{locus}{chr_name} && 
                                   $l_crispr->{start} eq $l_crispr_yaml->{locus}{chr_start} );

            my $right_identical = ( $r_crispr->{chr} eq $r_crispr_yaml->{locus}{chr_name} && 
                                    $r_crispr->{start} eq $r_crispr_yaml->{locus}{chr_start} );

            #if its the original pair, skip it
            if ( $left_identical and $right_identical ) {
                say "Skipping original pair entry.";
                next;
            }

            $num_valid_pairs{$pair_key}++;

            my $distance;
            if ( $r_crispr->{start} > $l_crispr->{start} ) {
                $distance = $r_crispr->{start} - $l_crispr->{end};
            }
            else {
                $distance = $l_crispr->{start} - $r_crispr->{end};
            }

            my $data = {
                exon_id       => $l_crispr->{exon_id}, #they will both have the same
                first_locus   => $l_crispr->{locus_str},
                second_locus  => $r_crispr->{locus_str},
                first_name    => $l_crispr->{name},
                second_name   => $r_crispr->{name},
                first_seq     => $l_crispr->{seq},
                second_seq    => $r_crispr->{seq},
                distance      => $distance,
                type          => $type,
            };

            #push @all, $data;

            #if we already have a shortest and this one is longer we can ignore it
            #we use abs so that -10 isnt "better" than 2
            next if defined $shortest && $distance > abs( $shortest->{distance} );
            
            #we'll set shortest if its not already set or distance is less. 
            $shortest = $data;
        }
    }

    #uncomment this and the above push to display summary data.

    # my @rows = sort { abs($a->{distance}) <=> abs($b->{distance}) } @all;
    
    # my $total_in_150 = 0;
    # say "l_name,l_locus,distance,r_name,r_locus";

    # for ( @rows ) {
    #     last if $_->{distance} > 9000;
    #     ++$total_in_150;

    #     say join ",", $_->{first_name}, $_->{first_locus}, $_->{distance}, $_->{second_name}, $_->{second_locus};
    # }

    # say "$filename";
    # say "\tTotal within 9000: " . scalar(@all);
    # say "\tTotal within 150: " . $total_in_150;
    # say "\tClosest: " . $rows[0]->{distance};

    #make sure we didnt somehow get more than 2 crispr names    
    die "More than 2 crisprs being compared?? ($filename):" . join(", ", keys %crispr_names) 
        if scalar keys %crispr_names > 2;

    #add the shortest one we found in this file to the global hash
    if ( defined $shortest ) {
        say STDERR "Closest crispr for " . join(" vs ", keys %crispr_names) . ": " 
                 . $shortest->{distance}; 
        say STDERR "Total valid: " . sum(values %num_valid_pairs) . " - " .
                   join(", ", map { $_ . ": " . $num_valid_pairs{$_} } keys %num_valid_pairs);

        my ( $f_name ) = $shortest->{first_name} =~ /_(.*)$/;
        my ( $s_name ) = $shortest->{second_name} =~ /_(.*)$/;

        $OFF_TARGET_DATA{ $shortest->{exon_id} }->{$f_name}{$s_name} =
            "{" . 
                'distance: "' . $shortest->{distance}     . '", ' . 
                'left: "'     . $shortest->{first_locus}  . '", ' . 
                'right: "'    . $shortest->{second_locus} . '", ' . 
                'type: "'     . $shortest->{type}         . '", ' .
                'id: "'       . "${f_name}_${s_name}"     . '", ' .
                'total_in_9k: "'.sum(values %num_valid_pairs). '", ' .
                'algorithm: "bwa"' .
            "}";
    }
    else {
        say STDERR "Couldn't find a valid paired off target for crisprs:" . join ", ", keys %crispr_names;
    }
}

#loop through pairs yaml and insert any off target data that we found
my $pair_data = LoadFile( $pair_yaml_file );
add_off_targets( $pair_data ); #edits the hash in place
DumpFile( $pair_yaml_file, $pair_data );



#
# METHODS
#

sub process_fq_file {
    my ( $crispr_fq_file ) = @_;

    my ( %crisprs, $current );

    open( my $fq_fh, "<", $crispr_fq_file ) or die "Couldn't open $crispr_fq_file: $!";
    while ( my $line = <$fq_fh> ) {
        chomp $line;

        #fq files have an @ line to signify a label, with the sequence on the next line
        if ( $line =~ /^@(.*)/ ) {
            $current = $1;
        }
        else {
            die "Error: there must only be one sequence per entry" if defined $crisprs{ $current };
            #pam must be at the end for bwa
            die "Crisprs shouldn't be reverse complemented!" unless $line =~ /GG$/; 

            #add fwd and reverse so we can compute hamming distances later.
            #we call it rev as we assume all crispr sequences have been reverse complemented.
            #cant automatically determine which is which cause CC-GG is a valid crispr
            $crisprs{$current} = { rev => revcom( $line )->seq, fwd => $line };
        }
    }

    return %crisprs;
}

sub process_windowbed_line {
    my ( $line ) = @_;

    #chr, start, end, name, unknown, strand x2
    my ( @data ) = split /\s+/, $line;
    die "Chromosomes don't match: $line" if $data[0] ne $data[6];

    my @bed_cols = qw( chr start end name score strand );

    #not first and second crispr cause we need to check the order first
    my ( %crispr_one, %crispr_two );
    @crispr_one{@bed_cols} = @data[ 0 .. 5 ]; #half the line is crispr one,
    @crispr_two{@bed_cols} = @data[ 6 .. 11 ]; #half is crispr two

    #we need to do some additional processing on each crispr entry and get the direction.
    #its the same for both crisprs so we do it in a loop
    for my $crispr ( \%crispr_one, \%crispr_two ) {
        $crispr->{start} += 1; #fix coords to match ensembl api, bed file number is different
        $crispr->{chr} =~ s/Chr//; #strip Chr prefix that human genome has on chromosome names
        $crispr->{locus_str} = $crispr->{chr} . ":" . $crispr->{start} . "-" . $crispr->{end}; #to detect dupes

        #split the name field into name and seq (name looks like ARA-CCGCTAGCTAGCTGATG) with a hash slice
        @{$crispr}{qw(name seq)} = delete($crispr->{name}) =~ /(.+)-([CGATcgat]+)/;

        #now we have the name without the seq extract the exon_id and crispr_id
        #the name is formatted like ENSMUSE000064698257_25A
        @{$crispr}{qw(exon_id crispr_id)} = $crispr->{name} =~ /^(ENS.*)_(.*)$/;

        die "Unknown name: " . $crispr->{name} unless exists $CRISPRS{ $crispr->{name} };

        $crispr->{direction} = get_direction( $crispr->{name}, $crispr->{seq} );
    }

    #this should never happen
    die "Pair is not in the same exon!" unless $crispr_one{exon_id} eq $crispr_two{exon_id};

    #make sure we got a direction for both
    return unless $crispr_one{direction} and $crispr_two{direction};

    #determine which is the left crispr and which is the right 

    my ( $left_crispr, $right_crispr );
    #which crispr is earlier in the chromosome? that is the first crispr
    #if they're the regular way round, everything is good.
    if ( $crispr_one{ start } < $crispr_two{ start } ) {
        $left_crispr = \%crispr_one;
        $right_crispr = \%crispr_two;
    }
    else {
        #the latter crispr comes first genomically
        #so we need them the other way as they're reversed
        $left_crispr = \%crispr_two;
        $right_crispr = \%crispr_one;
    }

    #do they overlap? if so we need to make sure that we have them the right way round:
    #in the fwd direction the pam is at the very end, so the pam site
    #will actually be further along than the rev pam site
    if ( $right_crispr->{ start } - $left_crispr->{ end } < 0 ) {
        if ( $left_crispr->{direction} eq "fwd" && $right_crispr->{direction} eq "rev" ) {
            #they're the wrong way around so swap them
            ( $left_crispr, $right_crispr ) = ( $right_crispr, $left_crispr );
        }
    }

    return $left_crispr, $right_crispr;
}

sub hamming_distance {
    #use string xor to get the number of mismatches between the two strings.
    #the xor returns a string with the binary digits of each char xor'd, 
    #which will be an ascii char between 001 and 255. tr returns the number of characters replaced. 
    die "Strings passed to hamming distance differ" if length($_[0]) != length($_[1]);
    return (uc($_[0]) ^ uc($_[1])) =~ tr/\001-\255//;
}

sub get_pair_exon_id {
    my ( $pair_data ) = @_;

    my ( $l_exon_id ) = $pair_data->{first_name} =~ $CRISPR_NAME_RE;
    my ( $r_exon_id ) = $pair_data->{second_name} =~ $CRISPR_NAME_RE;

    die "Pair is not in the same exon!" unless $l_exon_id eq $r_exon_id;

    return $l_exon_id;
}

sub get_direction {
    my ( $name, $seq ) = @_;

    my $fwd_distance = hamming_distance($seq, $CRISPRS{$name}->{fwd});
    my $rev_distance = hamming_distance($seq, $CRISPRS{$name}->{rev});

    #i hope this doesn't happen but its possible
    #die "Identical hamming distances found for both orientations of $seq" if $fwd_distance == $rev_distance;

    #make sure that we disregard ones we definitely cant identify
    if ( $fwd_distance == $rev_distance && $seq =~ /^CC.*GG$/ ) {
        say STDERR "Dual PAM and identical hamming distances for both orientations of $seq, skipping.";
        return;
    }   
    
    #try and get a direction
    if ( $fwd_distance <= $rev_distance && $seq =~ /GG$/ ) {
        #all seqs at this point should have a valid pam site, if not we complain.
        #die "$name does not have a valid pam: $seq" unless $seq =~ /GG$/;
        return if hamming_distance(substr($seq, 0, 20), substr($CRISPRS{$name}->{fwd},0,20)) > 5;
        return "fwd";
    }
    elsif ( $fwd_distance >= $rev_distance && $seq =~ /^CC/ ) {
        #die "$name does not have a valid pam: $seq" unless $seq =~ /^CC/;
        return if hamming_distance(substr($seq, 3), substr($CRISPRS{$name}->{rev},3)) > 5;
        return "rev";
    }
    else {
        say STDERR "Crispr doesn't look valid, skipping ($fwd_distance $rev_distance $seq)";
        return;
    }
}

sub add_off_targets {
    my ( $pair_data ) = @_; #should be hashref to yaml data
    
    #if we didn't get an off target we notify the user by saying there were non within the search bounds
    my $no_off_target = '{distance: "9000+"}';

    while ( my ( $exon_id, $pairs ) = each %{ $pair_data } ) {
        
        for my $pair ( @{ $pairs } ) {
            #shorthands because the full name is too long
            my ( $l_name, $r_name ) = ( $pair->{left_crispr}, $pair->{right_crispr} );

            $pair->{off_target_summary} = ''; #clear any old off target summary data.

            #see if there's any off target data for this pair
            if ( defined $OFF_TARGET_DATA{$exon_id} ) {
                #see if we have any off target data for this entry (we have to try both ways just in case)

                #loop through all possibilities, e.g for 2A and 5B: 2A 5B, 5B 2A, 2A 2A, 5B 5B
                my $found = 0;
                for my $first ( $l_name, $r_name ) {
                    next unless defined $OFF_TARGET_DATA{$exon_id}->{$first};

                    for my $second ( $l_name, $r_name ) {
                        if ( defined $OFF_TARGET_DATA{$exon_id}->{$first}{$second} ) {
                            #if there's already an off target summary see if we have a shorter one
                            if ( $pair->{off_target_summary} && $pair->{off_target_summary} !~ /9000\+/ ) {
                                $found = 1;
                                my ( $distance1 ) = $pair->{off_target_summary} =~ /distance:\s+"(-?\d+)",/;
                                my ( $distance2 ) = $OFF_TARGET_DATA{$exon_id}->{$first}{$second} =~ /distance:\s+"(-?\d+)",/;

                                next if abs($distance1) < abs($distance2);

                                #we found a shorter one so add it
                                $pair->{off_target_summary} = $OFF_TARGET_DATA{$exon_id}->{$first}{$second};
                                say STDERR "Found shorter distance.";
                            }
                            else {
                                $pair->{off_target_summary} = $OFF_TARGET_DATA{$exon_id}->{$first}{$second};
                                $found = 1;

                                #say STDERR "found $first $second, set to " . $OFF_TARGET_DATA{$exon_id}->{$first}{$second};
                            }
                        }
                    }
                }

                next if $found;

                unless ( $found ) {
                    say STDERR "Couldn't find $l_name vs $r_name pair in off target data.";
                }

                # if ( defined $OFF_TARGET_DATA{$exon_id}->{$l_name}{$r_name} ) {
                #     #we do have data so add it
                #     $pair->{off_target_summary} = $OFF_TARGET_DATA{$exon_id}->{$l_name}{$r_name};
                #     next;
                # }
                # elsif ( defined $OFF_TARGET_DATA{$exon_id}->{$r_name}{$l_name} ) {
                #     $pair->{off_target_summary} = $OFF_TARGET_DATA{$exon_id}->{$r_name}{$l_name};
                #     next;
                # }
                # else {
                #     say STDERR "Couldn't find $l_name vs $r_name pair in off target data.";
                #     say STDERR Dumper( \%OFF_TARGET_DATA );
                # }
            }
            else {
                say STDERR "$exon_id not defined in off target data hash";
            }

            #if we get here there's no data
            $pair->{off_target_summary} = $no_off_target;
        }
    }

    #say Dumper( \%OFF_TARGET_DATA );
}

1;

__END__

=head1 NAME

process_paired_data.pl - process windowbed output and merge the summary into the pairs file

=head1 SYNOPSIS

process_paired_data.pl [options]

    --fq-file            location of the crispr fq data created by find_paired_crisprs.pl
    --pair-yaml-file     location of the pair yaml data created by find_paired_crisprs.pl
    --crispr-yaml-file   location of crispr yaml data
    --paired-output      location of windowbed output files (can be multiple)
    --help               show this dialog

Example usage:

process_paired_data.pl --fq-file /nfs/users/nfs_a/ah19/work/crisprs/hprt_crisprs.fq --crispr-yaml-file /nfs/users/nfs_a/ah19/work/crisprs/hprt_crisprs.yaml --pair-yaml-file /nfs/users/nfs_a/ah19/work/crisprs/hprt_pairs.yaml --paired-output /nfs/users/nfs_a/ah19/work/crisprs/paired_output/*

=head1 DESCRIPTION

Iterate through all the paired output files to determine the closest off targets for the
given crispr pairings. The shortest for each pair is then written as a string into the original
pairs yaml file, ready to be persisted by persist_pairs.pl

=head AUTHOR

Alex Hodgkins

=cut
