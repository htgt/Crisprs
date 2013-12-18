#!/usr/bin/perl

use strict;
use warnings;

use DBI;
use feature qw( say );

die "Usage: wge_test.pl <bed file>" unless @ARGV;

my ( $db, $host, $port ) = ( 'wge', 'wge-db', 5448 );
my ( $user, $password ) = ( "user", "pass" ); # YOU NEED TO ADD THESE IN
my $dbh = DBI->connect( "dbi:Pg:dbname=$db;host=$host;port=$port", $user, $pass )
            or die "Couldn't connect to db.";

my $base_query = <<EOT;
WITH bed(chr_name, chr_start) as ( 
    VALUES %s
) 
SELECT c.id, c.seq, c.pam_right FROM bed 
JOIN crisprs c ON c.chr_name=bed.chr_name AND c.chr_start=bed.chr_start
WHERE c.species_id=1;
EOT

# my $where_query = <<EOT;
# SELECT id, seq, pam_right FROM crisprs 
# WHERE species_id=1 AND (%s);
# EOT

# my $values = <<EOT;
# ('1'::text, 15595::int ), 
# ('1', 43060 )
# EOT

#for wensheng_crisprs.bed 268.8s (1s on re-entry) -- with temp table/join
#for 52k rows of CBX1, only took 1 minute to insert?? half of the row were valid though

 my $values = join ",", 
                 map { '(' . $_->[0] . ',' . $_->[1] . ')'  } 
                     @{ process_bed( ) };


my $res = $dbh->selectall_arrayref( sprintf($base_query, $values) );

say "Found " . scalar( @{$res} ) . " rows.";
say $_ for @{ $res->[0] };

sub process_bed {
    my @data;
    my $count = 0;
    while ( <> ) {
        my @cols = split /\s+/, $_;

        $cols[0] =~ s/^Chr//;
        #first column is chromosome, second is start
        push @data, [ "'" . $cols[0] . "'", $cols[1] ];

        last if $count++ > 52103;
    }

    die "No data found in bed file." unless @data;

    say "Total lines: " . scalar( @data );

    #we have to force the types on the first entry
    #so postgres knows what they should be
    $data[0]->[0] .= '::text';
    $data[0]->[1] .= '::int';

    return \@data;
}

1;
