/*
    get_all_crisprs
    Copyright (C) 2013 Genome Research Limited

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

/* Modified version of scanham originally written by German Tischler (gt1@sanger.ac.uk) in 2013 */

/*
    this is an awful hack. it should take an option to split the file, species, etc.
    it should use 2bit instead of a list of chars, etc. etc.
    its bad.
*/

#include <iostream>
#include <cstdlib>
#include <memory>
#include <vector>
#include <list>
#include <map>
#include <string>
#include <fstream>
#include <stdexcept>
#include <new>
#include <sstream>
#include <cstdint>
#include <cassert>

int species_id;

void println(std::list<char> & current, std::string & seqname, int64_t start, int pam_right) {
    std::cout << seqname << "," << start << ","; //chr
    //display seq
    for ( auto it=current.begin(); it != current.end(); ++it ) {
        std::cout << *it;
    }

    //the final 1 is the species id
    std::cout << "," << pam_right << "," << species_id << '\n';
}

int main(int argc, char * argv[]) {
    if ( argc < 3 ) {
        std::cerr << "[U] usage: " << argv[0] << " <species_id> <text>" << std::endl;
        return EXIT_FAILURE;
    }

    species_id = atoi( argv[1] );
    
    // current match end position
    std::string seqname;
    int64_t seqpos = 0;
    int64_t patlen = 23;

    // character mapping table
    uint8_t cmap[256];
    // fill character mapping table
    std::fill(&cmap[0],&cmap[sizeof(cmap)/sizeof(cmap[0])],4);
    cmap['a'] = cmap['A'] = 0;
    cmap['c'] = cmap['C'] = 1;
    cmap['g'] = cmap['G'] = 2;
    cmap['t'] = cmap['T'] = 3;

    // space table
    bool smap[256];
    // fill space table
    std::fill(&smap[0],&smap[sizeof(cmap)/sizeof(smap[0])],0);
    for ( unsigned int i = 0; i < 256; ++i )
        if ( isspace(i) )
            smap[i] = 1;
    
    

    // open the text file
    std::ifstream text (argv[2]);
    if ( ! text.is_open() ) {
        std::cerr << "[D] cannot open file " << argv[2] << std::endl;
        throw std::runtime_error("Cannot open text file");
    }

    uint64_t total = 0;

    //bool skip = true;
    //int64_t num_done = 0;

    std::list<char> current (23, 'N');
    std::string line;

    while ( std::getline(text, line) ) {
        //make sure the line isn't blank
        if ( line.size() ) {
            // start of new sequence?
            if ( line[0] == '>' ) {
                //extract seqname as all the textup to the first space,
                //not including the >
                seqname = line.substr(1, line.size()-1);
                std::string::size_type first = seqname.find(" ");

                if ( first != std::string::npos ) 
                    seqname = seqname.substr(0, first);

                //strip Chr if its the first 3 characters
                if ( seqname.find("Chr") == 0 )
                    seqname = seqname.substr(3, seqname.size()-1);

                seqpos = 0;

                /*
                    NOTE: we should reset current here or we might accidentally
                          look at sequence spanning two chromosomes
                */

                //change this to choose which half you want
                //set to > 10 to get everything before and including 10
                //set to <= 10 to get everything AFTER 10
                // skip = ( ++num_done > 10 );
                // if ( skip )
                //    std::cerr << "Skipping chromosome " << seqname << std::endl;
                // else
                //    std::cerr << "Processing chromosome " << seqname << std::endl;

                std::cerr << "Processing chromosome " << seqname << std::endl;
                continue;
            }

            //if ( skip )
            //    continue;

            // scan the line
            for ( std::string::size_type i = 0; i < line.size(); ++i ) {
                // next character
                uint8_t const c = line[i];

                //skip if its a whitespace character
                if ( smap[c] )
                    continue;
                
                //remove first element and add new char to the end
                //should use an integer instead it would be quicker
                current.pop_front();
                current.push_back(c);

                // if we have a full pattern length worth of text
                if ( ++seqpos >= patlen ) {
                    int64_t seq_start = seqpos - patlen;

                    //we do this as two separate checks because if a CRISPR
                    //is CC-GG we want to print it twice.

                    //check if this is a valid crispr
                    std::list<char>::iterator start = current.begin();
                    if ( *(start) == 'C' && *(++start) == 'C' ) {
                        total++;

                        //last field is pam_right
                        println(current, seqname, seq_start+1, 0);
                    }

                    std::list<char>::iterator end = current.end();
                    if ( *(--end) == 'G' && *(--end) == 'G' ) {
                        total++;
                        println(current, seqname, seq_start+1, 1);
                    }
                }
            }
        }
    }

    text.close();

    std::cerr << "Found a total of " << total << " crisprs" << std::endl;
}
