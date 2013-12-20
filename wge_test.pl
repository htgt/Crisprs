#!/usr/bin/perl

use strict;
use warnings;

use DBI;
use feature qw( say );

die "Usage: wge_test.pl <bed file>" unless @ARGV;

my ( $db, $host, $port ) = ( 'wge', 'wge-db', 5448 );
my ( $user, $password ) = ( "user", "pass" ); # YOU NEED TO ADD THESE IN
my $dbh = DBI->connect( "dbi:Pg:dbname=$db;host=$host;port=$port", $user, $pass, { AutoCommit => 0 } )
            or die "Couldn't connect to db.";

#format the data as tsv
my $data = join "\n",
              map { $_->[0] . "\t" . $_->[1] }
                  @{ process_bed() };

#make temp table, give it the tsv data
$dbh->do( "CREATE TEMP TABLE bed (chr_name TEXT, chr_start INTEGER) ON COMMIT DROP;" );
$dbh->do( "COPY bed (chr_name, chr_start) FROM STDIN" );
$dbh->pg_putcopydata( $data );
$dbh->pg_putcopyend();

#join temporary table to the crispr table with explicit species for now
my $res = $dbh->selectall_arrayref( <<'EOT' );
SELECT c.id, c.seq, c.pam_right, c.chr_start, c.chr_name FROM bed 
JOIN crisprs c ON c.chr_name=bed.chr_name AND c.chr_start=bed.chr_start
WHERE c.species_id=1;
EOT

$dbh->rollback;
$dbh->disconnect;

say "Found " . scalar( @{$res} ) . " rows.";
say $_ for @{ $res->[0] }; #show the first entry

sub process_bed {
    my @data;
    my $count = 0;
    while ( <> ) {
        my @cols = split /\s+/, $_;

        $cols[0] =~ s/^Chr//;
        #first column is chromosome, second is start
        push @data, [ "'" . $cols[0] . "'", $cols[1] ];

        #last if $count++ > 52103;
    }

    die "No data found in bed file." unless @data;

    say "Total lines: " . scalar( @data );

    #we have to force the types on the first entry
    #so postgres knows what they should be
    #$data[0]->[0] .= '::text';
    #$data[0]->[1] .= '::int';

    return \@data;
}

1;
