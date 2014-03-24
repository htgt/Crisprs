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
#include <bitset>

using namespace std;

uint8_t cmap[256];
const unsigned int max_mismatches = 4;
unsigned int max_offs = 2000;
uint32_t version = 2;

typedef struct {
    uint64_t num_seqs;
    uint64_t seq_length;
    uint8_t species_id;
    char species[30]; //use fixed char arrays so we don't have to store size
    char assembly[30];
} metadata_t;

//crisprs we want to compare are inserted into this
typedef struct {
    uint64_t id;
    uint64_t seq;
    uint64_t rev_seq;
} crispr_t;

void populate_cmap() {
    //create an array for all possible char values, with 4 as the default value
    fill( &cmap[0], &cmap[sizeof(cmap)/sizeof(cmap[0])], 4 );

    //fill in the 8 entries we're actually interested in with their two bit values
    cmap['a'] = cmap['A'] = 0; //00
    cmap['c'] = cmap['C'] = 1; //01
    cmap['g'] = cmap['G'] = 2; //10
    cmap['t'] = cmap['T'] = 3; //11
}

//count the number of bits set to 1
template<unsigned int k>
struct PopCount {
    static unsigned int popcnt(uint64_t x) {
        //as everything is two bit we must convert them all to one bit,
        //to do this we must turn off all MSBs, but before we can do that
        //we need to ensure that when an MSB is set to 1, the LSB is also set.
        //the 4 will change pam_right to 0 so it doesnt get counted.
        //5 is 0101, 4 is 0100
        //x = (x | (x >> 1)) & (0x5555545555555555ull);
        x = (x | (x >> 1)) & (0x5555555555555555ull);

        x = x-( (x>>1) & 0x5555555555555555ull);
        x = (x & 0x3333333333333333ull) + ((x >> 2) & 0x3333333333333333ull);
        x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0full;
        return  (0x0101010101010101ull*x >> 56);
    }
};

uint64_t revcom(uint64_t text, int size) {
    //for a size of 23 we & with 23 1s to undo the complement on the part of
    //the integer we aren't interested in, for consitency
    //i.e. set the unused bits back to 0.
    //assumes size is smaller than sizeof(text)...
    unsigned int num_bits = sizeof(text) * CHAR_BIT;
    //we have to -1 from here to account for the pam_right bit which we DO want flipped
    uint64_t mask = 0xFFFFFFFFFFFFFFFFull >> ( (num_bits - (size * 2)) - 1 );
    //bit complement here maps A -> T, G -> C etc.
    text = ~text & mask;

    uint64_t reversed = 0;
    int shift = 0; 

    //now reverse the sequence, 2 bits at a time
    //could just do shift = 0; shift < size*2; shift += 2
    for ( int i = 0; i < size; i++, shift += 2 ) {
        reversed <<= 2;
        reversed |= ( text >> shift ) & 0x3;
    }

    return reversed;
}


string bits_to_string(uint64_t text, int size) {
    //have to & 0x3 to turn off all bits but what we actually want.
    string s (size, ' '); //make an empty string of the right size
    int shift = 2 * ( size - 1 ); //there are twice as many bits as there are characters

    //extract each character from the text
    for ( int i = 0; i < size; i++, shift -= 2 ) {
        //put the character we're interested in at the very end
        //of the integer, and switch all remaining bits to 0 with & 0x3
        uint8_t character = (text >> shift) & 0x3;
        switch ( character ) {
            case 0: s[i] = 'A'; break;
            case 1: s[i] = 'C'; break;
            case 2: s[i] = 'G'; break;
            case 3: s[i] = 'T'; break;
            default: break;
        }
    }

    return s;
}

uint64_t string_to_bits(const string & seq, short base = 0) {
    //loop through each character
    //optionally allow the user to prepend the sequence with a short
    uint64_t bits = base;
    for ( string::size_type j = 0; j < seq.size(); j++ ) {
        uint8_t const c = seq[j]; //get char

        if ( cmap[c] == 4 ) {
            //cerr << "Skipping N: " << seq << "\n";
            bits = 0; //set to all 0s
            break;
        }
        else {
            bits <<= 2; //shift left to make room for new char
            bits |= cmap[c]; //add our new char to the end 
        }
    }

    return bits;
}

vector<string> split(const string & text) {
    stringstream stream (text);
    string segment;
    vector<string> items;

    while( getline(stream, segment, ',') ) {
        items.push_back(segment);
    }

    return items;
}

void print_binary(uint64_t data) {
    bitset<64> a (data);
    cerr << data << "\n";
}

//apparently stringstreams are slow
//it would be quicker to create a string with a reserved size,
//filling in as necessary.
template <typename T>
string array_to_string(T begin, T end, bool include_idx) {
    ostringstream oss;

    //if its empty return
    if ( begin == end )
        return oss.str();

    oss << "{" << (include_idx ? "0: " : "")  << *(begin++);

    //iterate over array, appending each value to the string
    for ( int i = 1; begin < end; ++i ) {
        oss << ",";
        if ( include_idx )
            oss << " " << i << ": ";
        oss << *(begin++);
    }

    oss << "}";

    return oss.str();
}

//view this with:
//xxd -c 8 -b -l 10000 /lustre/scratch110/sanger/ah19/chr1-10_crisprs.dat | less
void text_to_binary(const vector<string> & infiles, const string & outfile) {
    ofstream out ( outfile, ios::binary );

    //write out a 1 as the first 8 bits. we will read this back in, and make sure its a 1
    uint8_t endian_test = 1;
    out.write( (char *) &endian_test, sizeof(endian_test) );
    out.write( (char *) &version, sizeof(version) );

    size_t offset = sizeof(endian_test) + sizeof(version);

    cerr << "Writing metadata\n";

    //this should really be passed to this method, and we add the last fields
    metadata_t data;
    strcpy(data.assembly, "GRCh37");
    strcpy(data.species, "Human");
    data.species_id = 1;
    data.num_seqs = 0;
    data.seq_length = 20;

    //leave space at the beginning for the number of sequences and the sequence length
    out.seekp( offset + sizeof(data), ios::beg );

    uint64_t total_skipped = 0;

    cerr << "Version is " << version << "\n";

    for ( vector<string>::size_type i = 0; i < infiles.size(); i++ ) {
        cerr << "Processing " << infiles[i] << "\n";
        ifstream text ( infiles[i] );

        /*
            NOTE: apparently C style FILE* access with fread is quicker
        */

        if ( text.is_open() && out.is_open() ) {
            string line;
        
            while ( getline( text, line ) ) {
                if ( line.size() ) {
                    //split the line up so we can get the parts we want
                    vector<string> spl = split(line);

                    //treated as a boolean
                    short pam_right = stoi( spl[3] );
                    if ( ! (pam_right == 0 || pam_right == 1) )
                        throw runtime_error( "pam_right field must be 1 or 0" );

                    //if its pam right remove last 3 chars, if its pam left remove first 3
                    string seq = pam_right ? spl[2].substr(0, 20) : spl[2].substr(3);

                    if ( seq.length() != data.seq_length )
                        throw runtime_error( "Different seq lengths in file!" );

                    //get a binary representation of the string and write it to the file
                    uint64_t bits = string_to_bits( seq, pam_right );
                    out.write( (char *) &bits, sizeof(bits) );

                    //we use a string of all A as an error string
                    //which is fine for now as crisprs need CC/GG.
                    if ( bits == 0 )
                        total_skipped++;

                    //keep track of how many we have so we can write it to the file
                    if ( ++data.num_seqs % 50000000 == 0 ) { 
                        cerr << "Converted " << data.num_seqs << " sequences\n";
                    }
                }
            }
            text.close();
        }
        else {
            cerr << "Error opening files: " << infiles[i] << ", " << outfile << "\n";
            throw runtime_error("Failed to open files");
        }
    }

    //write the metadata after the endianness and version numbers
    out.seekp( offset, ios::beg );
    out.write( (char *) &data, sizeof(data) );

    out.close();

    cerr << "Sequence length is " << data.seq_length << endl;
    cerr << "Converted " << data.num_seqs << " sequences" << endl;
    cerr << "Skipped " << total_skipped << " sequences" << endl;
}

void load_binary(const string & filename) {
    ifstream text ( filename );

    if ( text.is_open() ) {
        //make sure that the binary file has the same endianness as our hardware
        uint8_t endian_test;
        text.read( (char *) &endian_test, sizeof(uint8_t) );

        if ( endian_test != 1 ) {
            throw runtime_error("Endianness of the file does not match your hardware!");
        }

        //make sure the versions matched or things will go very badly.
        uint32_t file_version;
        text.read( (char *) &file_version, sizeof(uint32_t) );
        
        if ( file_version != version )
            throw runtime_error("File is the wrong version! Please regenerate.");

        cerr << "Version is " << version << "\n";

        //load the metadata
        metadata_t data;
        text.read( (char *) &data, sizeof(data) );

        cerr << "Loading CRISPRs for " << data.assembly << " (" << data.species << ")" << "\n";
        cerr << "File has " << data.num_seqs << " sequences" << "\n";
        cerr << "Sequence length is " << data.seq_length << "\n";
        cerr << "Species id is " << int(data.species_id) << endl;

        uint64_t memory_required = data.num_seqs * sizeof(uint64_t);
        cerr << "Will require " << (memory_required/1024)/1024 << "MB of memory" << endl;

        //dont use more than 3gb of memory, safety thing for now.
        if ( memory_required > 3221225472 ) {
            throw runtime_error("More than 3gb of memory required, aborting.");
        }

        //allocate heap memory as its too big for stack
        uint64_t * sequences = new uint64_t [data.num_seqs];
        for ( uint64_t i = 0; i < data.num_seqs; i++ ) {
            if ( i % 50000000 == 0 ) {
                cerr << "Loaded " << i << " sequences" << endl;
            }
            //load each 'integer'
            text.read( (char *) &sequences[i], sizeof(uint64_t) );
        }

        cerr << "Loaded " << data.num_seqs << " sequences" << endl;

        text.close();

        //this will be in another method

        clock_t t;

        //we'll know the size of these, so should make them arrays
        vector<crispr_t> queries;

        uint64_t start = 263164224;
        for ( uint64_t i = 0; i < 2; i++ ) {
            crispr_t crispr;
            crispr.id = start+i;
            crispr.seq = sequences[crispr.id];
            crispr.rev_seq = revcom(sequences[crispr.id], data.seq_length);

            queries.push_back( crispr );
        }

        //queries.push_back( sequences[263164224] );
        //queries_rev.push_back( revcom(sequences[263164224], data.seq_length) );

        //cout << bits_to_string( queries[0], data.seq_length ) << "\n\n";

        //use inside the loop
        unsigned int mm;
        
        //this can be used to check if the pam_right bit is set.
        //the complement of it can be used to turn off the flag.
        uint64_t pam_on = 0x1ull << data.seq_length*2;
        uint64_t pam_off = ~pam_on;

        array<uint64_t, max_mismatches+1> summary;

        cerr << "Searching for off targets\n";

        t = clock();
        //calculate each query individually
        for ( vector<crispr_t>::size_type i = 0; i < queries.size(); i++ ) {
            //make vector with 5000 slots reserved
            vector<uint64_t> off_targets;
            off_targets.reserve( max_offs );
            
            //new crispr so reset the summary to all 0
            summary.fill(0);

            uint64_t total_matches = 0;

            //iterate over every crispr checking for of targets
            for ( uint64_t j = 0; j < data.num_seqs; j++ ) {
                //xor the two bit strings
                //try forward first
                uint64_t match = queries[i].seq ^ sequences[j];              
                

                //we only care about the orientation that matches the orientation
                //of the sequence. the other can safely be ignored.
                //one of these will always be true as its a single bit

                //if this is true its because the pams did NOT match.
                //when they don't match we need to use the reverse complemented sequence,
                //which will match.
                if ( match & pam_on ) {
                    //orientations didnt match so we need to try the reverse
                    uint64_t match_r = queries[i].rev_seq ^ sequences[j];
                    mm = PopCount<sizeof(unsigned long)>::popcnt(match_r & pam_off);
                }
                else {
                    //the PAMs matched so count the mismatches            
                    mm = PopCount<sizeof(unsigned long)>::popcnt(match & pam_off);
                }

                if ( mm <= max_mismatches ) {
                    summary[mm]++;

                    if ( ++total_matches < max_offs) {
                        off_targets.push_back( j + 1 );
                    }
                }
            }

            cerr << "Found " << total_matches << " off targets.\n";

            cerr << queries[i].id << "\n";

            //print off targets
            if ( total_matches < max_offs ) {
                cerr << array_to_string(off_targets.begin(), off_targets.end(), 0) << "\n";
            }

            //string test = array_to_string(summary.begin(), summary.end(), 1);

            cerr << array_to_string(summary.begin(), summary.end(), 1) << "\n";

            //cerr << "After\n";

            // print the off target summary

            //print opening brace and first item, so we can print the commas
            //cerr << "{" << summary[0];
            //for ( unsigned int i = 1; i <= max_mismatches; i++ ) {
            //    cerr << ", " << i << ": " << summary[i];
            //}
            //cerr << "}\n";

            //cerr << "Found " << total_matches << " matches" << endl;

            //cout << "\n";
        }

        t = clock() - t;

        fprintf( stderr, "Took %f seconds\n", ((float)t)/CLOCKS_PER_SEC );

        //cerr << "fwd:" << tot_fwd << "\n" << "rev:" << tot_rev << "\n";

        delete[] sequences;
    }
    else {
        cerr << "Error opening file" << filename << endl;
        throw runtime_error("Failed to open input file");
    }
}

int main(int argc, char * argv[]) {
    populate_cmap(); //used for base to bit conversion

    vector<string> infiles = {  
        "/lustre/scratch110/sanger/ah19/chr1-10_crisprs_new.csv",
        "/lustre/scratch110/sanger/ah19/chr11_onwards_crisprs_new.csv"
    };

    //string infile = "/lustre/scratch109/sanger/ah19/chr1-10_crisprs_new.csv";
    string outfile = "/lustre/scratch110/sanger/ah19/chr1-10_crisprs.dat";
    //text_to_binary( infiles, outfile );

    load_binary( outfile );

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
*/

/* 
    Written by Alex Hodgkins (ah19@sanger.ac.uk) in 2014 
    Some code taken from scanham written by German Tischler
*/
