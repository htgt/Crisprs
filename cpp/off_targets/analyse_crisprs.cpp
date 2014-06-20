#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <vector>
#include <array>
#include <stdexcept>
#include <ctime>
#include <climits>
#include <unistd.h>
#include <map>

#include "utils.h"
#include "crisprutil.h"

using namespace std;

/*
build with:
g++ -std=c++0x -O3 -W -Wall crisprutil.o utils.o analyse_crisprs.cpp -o analyse_crisprs
*/

/*
should add finding of all crisprs from genome into these options
maybe a method to find an id given a sequence?
*/
int usage() {
    fprintf(stderr, "\n");
    fprintf(stderr, "Program: find_off_targets\n");
    fprintf(stderr, "Contact: Alex Hodgkins <ah19@sanger.ac.uk>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: find_off_targets <command> [options]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Command: index      Create binary index of all CRISPRs\n");
    fprintf(stderr, "         align      Find potential off targets for CRISPRs\n");
    fprintf(stderr, "\n\n");
    fprintf(stderr, "Example end to end usage:\n");
    fprintf(stderr, "Create index:\n");
    fprintf(stderr, "\tfind_off_targets index -i human_chr1-11.csv -i human_chr12_on.csv -o index.bin\n");
    fprintf(stderr, "Calculate off targets for single CRISPR:\n");
    fprintf(stderr, "\tfind_off_targets align 873245 > crispr_data.tsv\n");

    return 1;
}

int align_usage() {
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: find_off_targets align [options] <ids>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options: -i FILE    The file containing the CRISPR index, (from the index step)\n");
    fprintf(stderr, "         -s int     The CRISPR range start\n");
    fprintf(stderr, "         -n int     How many CRISPRs to compute after start)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Note:\n");
    fprintf(stderr, "You can specify a list of ids, OR specify a range of CRISPRs. ");
    fprintf(stderr, "For example:\n\n");
    fprintf(stderr, "./find_off_targets align -s 43275 -n 1000\n");
    fprintf(stderr, "  This will calculate off targets for 1000 crisprs,\n");
    fprintf(stderr, "  the ids of which will be 43275-44725\n");
    fprintf(stderr, "\n\n");
    fprintf(stderr, "find_off_targets align 873245 923577 237587 109583\n");
    fprintf(stderr, "  This will calculate off targets for the given 4 CRISPRs\n");

    return 1;
}

int main(int argc, char * argv[]) {
    int c = -1;

    map<uint64_t, int> ids;
    string index = "";

    while ( (c = getopt(argc, argv, "i:")) != -1 ) {
        switch ( c ) {
            case 'i': index = optarg; break;
            case '?': return align_usage();
        }
    }

    if ( index == "" ) {
        cerr << "An index file must be specified with the -i option\n";
        return align_usage();
    }

    CrisprUtil finder = CrisprUtil();
    //finder.load_binary( index );

    finder.analyse_binary( index );

    // //loop all sequences incrementing our map
    // for ( uint64_t i = 1; i < finder.num_seqs(); i++ ) {
    //     if ( i % 50000000 == 0 ) {
    //         cerr << "Loaded " << i << " sequences" << endl;
    //     }

    //     ids[finder.get_crispr_int(i)]++;
    // }

    // uint64_t seq_length = finder.seq_length();

    // for ( auto it = ids.begin(); it != ids.end(); ++it ) {
    //     if ( it->second > 10 )
    //         continue;

    //     cout << util::bits_to_string( it->first, seq_length )
    //          << " " << it->second;
    // }

    return 0;
}


/*
    Copyright (C) 2014 Genome Research Limited

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 
    Written by Alex Hodgkins (ah19@sanger.ac.uk) in 2014 
    Some code taken from scanham written by German Tischler
*/
