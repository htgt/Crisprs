/*
    scanham
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

/* Written by German Tischler (gt1@sanger.ac.uk) in 2013 */

#include <iostream>
#include <cstdlib>
#include <memory>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <stdexcept>
#include <new>
#include <sstream>
#include <cstdint>
#include <cassert>



/**
 * load patterns from a FastA style file
 **/
//std::vector<std::string> loadPatterns(std::string const & filename)
std::map<std::string, std::string> loadPatterns(std::string const & filename)
{
	std::map<std::string, std::string> M;
	std::ifstream istr(filename.c_str());
	if ( ! istr.is_open() )
	{
		std::cerr << "[D] Failed to open file " << filename << std::endl;
		throw std::runtime_error("Failed to open input file");
	}
	
	// buffer for building patterns from multi line input
	::std::shared_ptr<std::ostringstream> ostr(new std::ostringstream);
	std::string fname;
	while ( istr )
	{
		std::string line;
		std::getline(istr,line);
		if ( line.size() )
		{
			if ( line[0] == '@' || line[0] == '>' ) //fq or fa
			{
				if ( ostr->str().size() ) {
					M.insert( std::pair<std::string,std::string>(fname, ostr->str()) );
				}

				fname = line.substr(1, line.size()-1); //get fasta name
					
				ostr = ::std::shared_ptr<std::ostringstream>(new std::ostringstream);
			}
			else
			{
				(*ostr) << line;
			}
		}
	}
	
	//the very last entry wont have been added yet as there's no > line after it
	//fname will be correct though
	if ( ostr->str().size() ) {
		M.insert( std::pair<std::string,std::string>( fname, ostr->str()) );
	}
		
	std::map<std::string,std::string>::iterator it;
	//for ( std::vector<std::string>::size_type i = 0; i < V.size(); ++i )
	for ( it=M.begin(); it!=M.end(); ++it )
	{
		std::string pat = it->second; 
		std::ostringstream fstr;
		// keep non space symbols
		for ( std::string::size_type j = 0; j < pat.size(); ++j )
			if ( !isspace(pat[j]) )
				fstr << pat[j];
				
		pat = fstr.str();
		
		// unify
		for ( std::string::size_type j = 0; j < pat.size(); ++j )
			switch ( pat[j] )
			{
				case 'a': case 'A': pat[j] = 'A'; break;
				case 'c': case 'C': pat[j] = 'C'; break;
				case 'g': case 'G': pat[j] = 'G'; break;
				case 't': case 'T': pat[j] = 'T'; break;
				default: pat[j] = 'N'; break;
			}
	}
	
	return M;
}

static std::string masksToString(uint64_t const textmask, uint64_t const errmask, unsigned int const k)
{
	std::string s(k,' ');
	unsigned int shift = 2*(k-1);
	
	for ( unsigned int i = 0; i < k; ++i, shift -= 2 )
	{
		if ( (errmask >> shift) & 0x3 )
		{
			s[i] = 'N';
		}
		else
		{
			uint8_t const m = (textmask >> shift) & 0x3;
			switch ( m )
			{
				case 0: s[i] = 'A'; break;
				case 1: s[i] = 'C'; break;
				case 2: s[i] = 'G'; break;
				case 3: s[i] = 'T'; break;
				default: break;
			}
		}
	}
	
	return s;
}

//count the number of bits set to 1
template<unsigned int k>
struct PopCount
{
	static unsigned int popcnt(uint64_t x)
	{
		x = x-( (x>>1) & 0x5555555555555555ull);
		x = (x & 0x3333333333333333ull) + ((x >> 2) & 0x3333333333333333ull);
		x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0full;
		return  (0x0101010101010101ull*x >> 56);
	}
};

#if defined(__GNUC__)
template<>
struct PopCount<4>
{
	static unsigned int popcnt(uint64_t const v)
	{
		return
			__builtin_popcountl((v >> 32)&0xFFFFFFFFULL)
			+
			__builtin_popcountl((v >>  0)&0xFFFFFFFFULL);
	}
};

template<>
struct PopCount<8>
{
	static unsigned int popcnt(uint64_t const v)
	{
		return __builtin_popcountl(v);
	}
};
#endif

int main(int argc, char * argv[])
{
	try
	{
		if ( argc < 4 )
		{
			std::cerr << "[U] usage: " << argv[0] << " <patterns> <text> <mismatches>" << std::endl;
			return EXIT_FAILURE;
		}
	
		// parse number of mismatches	
		unsigned int h = 0;
		std::istringstream histr(argv[3]);
		histr >> h;
		
		if ( ! histr )
		{
			std::cerr << "[D] Cannot parse " << argv[3] << " as number." << std::endl;
			throw std::runtime_error("Cannot parse input parameter");
		}
		
		// load the pattern file
		std::map<std::string,std::string> patterns = loadPatterns(argv[1]);
		
		std::cerr << "[V] loaded " << patterns.size() << " patterns." << std::endl;

		//get first item to check sizes against
		std::map<std::string,std::string>::iterator first_pat = patterns.begin();

		// check that all patterns have the same length
		//for ( std::vector<std::string>::size_type i = 1; i < patterns.size(); ++i )
		std::map<std::string,std::string>::iterator it;
		for ( it=patterns.begin(); it!=patterns.end(); ++it )
			if ( it->second.size() != first_pat->second.size() )
			{
				std::cerr << "[D] Pattern length needs to be consistent." << std::endl;
				throw std::runtime_error("Inconsistent pattern size");
			}
		
		// we need at least one pattern
		if ( ! patterns.size() )
		{
			std::cerr << "[V] no patterns given." << std::endl;
			return EXIT_SUCCESS;
		}
		
		// check pattern length
		unsigned int const patlen = first_pat->second.size();
		if ( patlen > 32 )
		{
			std::cerr << "[D] Cannot handle pattern length of more than 32" << std::endl;
			throw std::runtime_error("Invalid pattern length");
		}
		if ( ! patlen )
		{
			std::cerr << "[D] Cannot handle empty patterns" << std::endl;
			throw std::runtime_error("Invalid pattern length");
		}
		else 
		{
			std::cerr << "[D] Pattern length is " << patlen << std::endl;
		}
		
		// pattern masks and error masks
		std::vector<uint64_t> patmasks;
		std::vector<uint64_t> paterr;
		std::vector<std::string> names;
		unsigned int const ishift = 2*(patlen-1);
		
		// build bit masks
		//for ( std::vector<std::string>::size_type i = 0; i < patterns.size(); ++i )
		for ( it=patterns.begin(); it!=patterns.end(); ++it )
		{
			uint64_t mask = 0, err = 0;
			uint64_t imask = 0, ierr = 0;
			std::string const & pat = it->second;
			
			for ( std::string::size_type j = 0; j < pat.size(); ++j )
			{
				mask <<= 2;
				err <<= 2;
				imask >>= 2; //push the other way so its the reverse complement
				ierr >>= 2;
				
				switch ( pat[j] )
				{
					case 'A': 
						mask |= 0;
						imask |= (3ull << ishift);
						break;
					case 'C': 
						mask |= 1; 
						imask |= (2ull << ishift);
						break;
					case 'G': 
						mask |= 2; 
						imask |= (1ull << ishift);
						break;
					case 'T':
						mask |= 3;
						imask |= (0ull << ishift);
						break;
					default: 
						err |= 1; 
						ierr |= (1ull << ishift);
						break;
				}
			}

			/*
				TOMORROW: take value from names array and put it in bed file output. finish bash script
			*/
			
			patmasks.push_back(mask);
			paterr.push_back(err);
			names.push_back(it->first); //wahtever
			patmasks.push_back(imask); //to check the reverse complement
			paterr.push_back(ierr);
			names.push_back(it->first);
		}
		
		//note: as far as i can see this doesn't do anything. 
		//we end up doing & against all 1s.
		
		// mask with 2*patlen lower bits set
		uint64_t andmask = 0;
		for ( unsigned int i = 0; i < patlen; ++i )
		{
			andmask <<= 2;
			andmask |= 0x3;
		}
		
		// current match end position
		int64_t seqid = -1;
		std::string seqname;
		int64_t seqpos = 0;
		// current bit mask of window on text
		uint64_t textmask = 0;
		uint64_t texterr = 0;
		// character mapping table
		uint8_t cmap[256];
		// space table
		bool smap[256];
		// fill character mapping table
		std::fill(&cmap[0],&cmap[sizeof(cmap)/sizeof(cmap[0])],4);
		cmap['a'] = cmap['A'] = 0;
		cmap['c'] = cmap['C'] = 1;
		cmap['g'] = cmap['G'] = 2;
		cmap['t'] = cmap['T'] = 3;
		// fill space table
		std::fill(&smap[0],&smap[sizeof(cmap)/sizeof(smap[0])],0);
		for ( unsigned int i = 0; i < 256; ++i )
			if ( isspace(i) )
				smap[i] = 1;
		
		// open the text file
		std::ifstream textistr(argv[2]);
		if ( ! textistr.is_open() )
		{
			std::cerr << "[D] cannot open file " << argv[2] << std::endl;
			throw std::runtime_error("Cannot open text file");
		}
		// read text file
		while ( textistr )
		{
			// get line
			std::string line;
			std::getline(textistr,line);
			// if line is not empty
			if ( line.size() )
			{
				// start of new sequence?
				if ( line[0] == '>' )
				{
					seqid++;
					seqname = line.substr(1, line.size()-1);
					std::string::size_type first = seqname.find(" ");

					if ( first != std::string::npos ) 
						seqname = seqname.substr(0, first);

					seqpos = 0;
					textmask = 0;
					texterr = 0;
				}
				else
				{
					// scan the line
					for ( std::string::size_type i = 0; i < line.size(); ++i )
					{
						// next character
						uint8_t const c = line[i];
						
						// if character is not white space
						if ( ! smap[c] )
						{
							// insert new symbol in text and error bit masks
							textmask <<= 2;
							texterr <<= 2;
							
							if ( cmap[c] == 4 )
							{
								texterr  |= 1;
							}
							else
							{
								assert ( cmap[c] < 4 );
								textmask |= cmap[c];
							}

							textmask &= andmask;
							texterr &= andmask;

							// if we have a full pattern length worth of text
							if ( ++seqpos >= patlen )
							{
								// iterate over patterns
								for ( std::vector<uint64_t>::size_type j = 0; j < patmasks.size(); ++j )
								{
									uint64_t const m = 
										// "regular" mismatches between text and pattern
										(patmasks[j] ^ textmask) 
										// indeterminate bases in text
										| texterr 
										// indeterminate bases in pattern
										| paterr[j];
								

									//look up SWAR algorithm to explain this kind of thing properly.
									//2 bits represents a single match, but if both bits are 'on' (if its a 3)
									//we only want to count it once.
									//what this does:
									//1. make sure the least significant bit is set if the most significant bit is
									//		(so that 10 becomes 01)
									//2. set all most significant bits to 0. then the number of 1s is the pop count
									uint64_t const m1 = (m | (m >> 1)) & (0x5555555555555555ull);
									
									// count number of mismatches + indeterminate bases
									unsigned int const mismatches = PopCount<sizeof(unsigned long)>::popcnt(m1);
										
									// if number is low enough, then print the match
									if ( mismatches <= h )
									{
										//std::cout << "match for pattern " << j/2 << " in chr" << seqname << " pos " << (seqpos-patlen) << " with " << mismatches << " mismatches, strand " << ((j&1)?'-':'+') << std::endl;
										//std::cout << "[M] " << masksToString(textmask,texterr,patlen) << std::endl;
										//could add the actual mismatch alignment here
										//std::cout << "[M] " << masksToString(patmasks[j],paterr[j],patlen) << std::endl;

										//need to change whatever to the seq name from the pattern fasta
										//store the mismatches in the score field- naughty
										std::cout << seqname << "\t" << (seqpos-patlen) << "\t" << seqpos << "\t" << names[j] << "-" << masksToString(textmask,texterr,patlen) << "\t" << mismatches << "\t" << ((j&1)?'-':'+') << std::endl;
									}
								}
							}
						}
						
					}
				}
			}
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << "[D] " << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
