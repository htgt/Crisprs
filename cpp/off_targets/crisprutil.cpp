#include <string>
#include <cstring>
#include <sstream>
#include <vector>
#include <array>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <ctime>
#include <climits>

#include "utils.h"
#include "crisprutil.h"

using namespace std;

CrisprUtil::CrisprUtil() {
    //initialise variables
    max_offs = 2000;

    //so we can see if data has been loaded
    crispr_data.num_seqs = 0;

    //populate the array we use to map char -> int
    _populate_cmap();    

    //cerr << bits_to_string( crisprs[263164224], 20 ) << "\n";
}

void CrisprUtil::_populate_cmap() {
    //create an array for all possible char values, with 4 as the default value
    std::fill( &cmap[0], &cmap[sizeof(cmap)/sizeof(cmap[0])], 4 );

    //fill in the 8 entries we're actually interested in with their two bit values
    cmap['a'] = cmap['A'] = 0; //00
    cmap['c'] = cmap['C'] = 1; //01
    cmap['g'] = cmap['G'] = 2; //10
    cmap['t'] = cmap['T'] = 3; //11
}

string CrisprUtil::get_crispr(uint64_t id) {
    //make sure load_binary has been called
    if ( crispr_data.num_seqs == 0 )
        throw runtime_error( "CRISPRs must be loaded before calling get_crispr" );

    if ( id >= crispr_data.num_seqs )
        throw runtime_error( "Index out of range" );

    return util::bits_to_string( crisprs[id], crispr_data.seq_length );
}

void CrisprUtil::load_binary(const string & filename) {
    cerr << "Loading binary data from:\n" << filename << "\n";
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
        
        if ( file_version != VERSION )
            throw runtime_error("File is the wrong version! Please regenerate.");

        cerr << "Version is " << VERSION << "\n";

        //load the metadata
        text.read( (char *) &crispr_data, sizeof(crispr_data) );

        cerr << "Loading CRISPRs for " << crispr_data.assembly 
             << " (" << crispr_data.species << ")" << "\n"
             << "File has " << crispr_data.num_seqs << " sequences" << "\n"
             << "Sequence length is " << crispr_data.seq_length << "\n"
             << "Species id is " << int(crispr_data.species_id) << endl;

        uint64_t memory_required = crispr_data.num_seqs * sizeof(uint64_t);
        cerr << "Will require " << (memory_required/1024)/1024 << "MB of memory" << endl;

        //dont use more than 3gb of memory, safety thing for now.
        if ( memory_required > 3221225472 ) {
            throw runtime_error("More than 3gb of memory required, aborting.");
        }

        clock_t t = clock();

        /*
            should consider creating a struct with bitfields around the uint64_t,
            as we can then name sections and stop passing the lengths around.
            so, something like:
            typdef struct {
                uint64_t seq:40, pam_right:1, blank:23;
            } crispr_seq_t;

            might be slower, however.
        */

        //allocate heap memory as its too big for stack
        //store pointer in the class variable
        crisprs = new uint64_t [crispr_data.num_seqs+1];
        crisprs[0] = 0; //we don't use 0 so that all ids match the db ids
        for ( uint64_t i = 1; i <= crispr_data.num_seqs; i++ ) {
            if ( i % 50000000 == 0 ) {
                cerr << "Loaded " << i << " sequences" << endl;
            }
            //load each 'integer'
            text.read( (char *) &crisprs[i], sizeof(uint64_t) );
        }

        t = clock() - t;

        cerr << "Loaded " << crispr_data.num_seqs << " sequences" << endl;

        text.close();

        fprintf( stderr, "Loading took %f seconds\n", ((float)t)/CLOCKS_PER_SEC );
    }
    else {
        cerr << "Error opening file" << filename << endl;
        throw runtime_error("Failed to open input file");
    }
}

//view this with:
//xxd -c 8 -b -l 10000 /lustre/scratch110/sanger/ah19/crisprs.bin | less
void CrisprUtil::text_to_binary(const vector<string> & infiles, const string & outfile, metadata_t * data) {
    ofstream out ( outfile, ios::binary );

    //write out a 1 as the first 8 bits. we will read this back in, and make sure its a 1
    uint8_t endian_test = 1;
    out.write( (char *) &endian_test, sizeof(endian_test) );
    out.write( (char *) &VERSION, sizeof(VERSION) );

    size_t offset = sizeof(endian_test) + sizeof(VERSION);

    cerr << "Writing metadata\n";

    //this should really be passed to this method, and we add the last fields
    // metadata_t data;
    // strcpy(data.assembly, "GRCh37");
    // strcpy(data.species, "Human");
    // data.species_id = 1;
    // data.num_seqs = 0;
    // data.seq_length = 20;

    //leave space at the beginning for the number of sequences and the sequence length
    out.seekp( offset + sizeof(*data), ios::beg );

    uint64_t total_skipped = 0;

    cerr << "Version is " << VERSION << "\n";

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
                    vector<string> spl = util::split(line);

                    //treated as a boolean
                    short pam_right = stoi( spl[3] );
                    if ( ! (pam_right == 0 || pam_right == 1) )
                        throw runtime_error( "pam_right field must be 1 or 0" );

                    //if its pam right remove last 3 chars, if its pam left remove first 3
                    string seq = pam_right ? spl[2].substr(0, 20) : spl[2].substr(3);

                    if ( seq.length() != data->seq_length )
                        throw runtime_error( "Different seq lengths in file!" );

                    //get a binary representation of the string and write it to the file
                    uint64_t bits = util::string_to_bits( cmap, seq, pam_right );
                    out.write( (char *) &bits, sizeof(bits) );

                    //we use a string of all A as an error string
                    //which is fine for now as crisprs need CC/GG.
                    if ( bits == 0 )
                        total_skipped++;

                    //keep track of how many we have so we can write it to the file
                    if ( ++data->num_seqs % 50000000 == 0 ) { 
                        cerr << "Converted " << data->num_seqs << " sequences\n";
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
    out.write( (char *) data, sizeof(*data) );

    out.close();

    cerr << "Sequence length is " << data->seq_length << endl;
    cerr << "Converted " << data->num_seqs << " sequences" << endl;
    cerr << "Skipped " << total_skipped << " sequences" << endl;
}

void CrisprUtil::find_off_targets(vector<uint64_t> ids) {
    vector<crispr_t> queries;

    for ( auto i = ids.begin(); i < ids.end(); i++ ) {
        crispr_t crispr;
        crispr.id = *i;
        crispr.seq = crisprs[crispr.id];
        crispr.rev_seq = util::revcom(crisprs[crispr.id], crispr_data.seq_length);

        queries.push_back( crispr );
    }

    _find_off_targets( queries );
}

void CrisprUtil::find_off_targets(uint64_t start, uint64_t amount) {
    vector<crispr_t> queries;

    //create a crispr_t for every requested id
    for ( uint64_t i = 0; i < amount; i++ ) {
        crispr_t crispr;
        crispr.id = start+i;
        crispr.seq = crisprs[crispr.id];
        crispr.rev_seq = util::revcom( crisprs[crispr.id], crispr_data.seq_length );

        queries.push_back( crispr );
    }

    _find_off_targets( queries );
}

void CrisprUtil::_find_off_targets(vector<crispr_t> queries) {
    unsigned int mm;
    
    //this can be used to check if the pam_right bit is set.
    //the complement of it can be used to turn off the flag.
    const uint64_t pam_on = 0x1ull << crispr_data.seq_length*2;
    const uint64_t pam_off = ~pam_on;

    array<uint64_t, max_mismatches+1> summary;

    cerr << "Searching for off targets\n";

    clock_t t = clock();
    //calculate each query individually
    for ( vector<crispr_t>::size_type i = 0; i < queries.size(); i++ ) {
        //make vector with 5000 slots reserved
        //should we do this outside the loop like summary?
        vector<uint64_t> off_targets;
        off_targets.reserve( max_offs );
        
        //new crispr so reset the summary to all 0
        summary.fill(0);

        uint64_t total_matches = 0;

        //iterate over every crispr checking for of targets
        for ( uint64_t j = 1; j <= crispr_data.num_seqs; j++ ) {
            //xor the two bit strings
            //try forward first
            uint64_t match = queries[i].seq ^ crisprs[j];         
            

            //we only care about the orientation that matches the orientation
            //of the sequence. the other can safely be ignored.
            //one of these will always be true as its a single bit

            //if this is true its because the pams did NOT match.
            //when they don't match we need to use the reverse complemented sequence,
            //which will match.
            if ( match & pam_on ) {
                //orientations didnt match so we need to try the reverse
                uint64_t match_r = queries[i].rev_seq ^ crisprs[j];
                mm = util::pop_count( match_r & pam_off );
            }
            else {
                //the PAMs matched so count the mismatches            
                mm = util::pop_count( match & pam_off );
            }

            if ( mm <= max_mismatches ) {
                summary[mm]++;

                if ( ++total_matches < max_offs ) {
                    off_targets.push_back( j );
                }
            }
        }

        cerr << "Found " << total_matches << " off targets.\n";

        cout << queries[i].id << "\n";

        //print off targets
        if ( total_matches < max_offs ) {
            cout << util::array_to_string(off_targets.begin(), off_targets.end(), 0) << "\n";
        }

        cout << util::array_to_string(summary.begin(), summary.end(), 1) << "\n\n";
    }

    t = clock() - t;

    fprintf( stderr, "Took %f seconds\n", ((float)t)/CLOCKS_PER_SEC );
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