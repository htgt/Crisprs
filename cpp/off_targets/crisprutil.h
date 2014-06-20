#ifndef __CRISPR_UTIL__
#define __CRISPR_UTIL__

#define MAX_CHAR_SIZE 30

const uint32_t VERSION = 3;

#include <string>
#include <vector>

//struct that is stored in the binary index
typedef struct {
    uint64_t num_seqs;
    uint64_t seq_length;
    uint64_t offset; //so we can give real mouse db ids
    uint8_t species_id;
    char species[MAX_CHAR_SIZE]; //use fixed char arrays so we don't have to store size
    char assembly[MAX_CHAR_SIZE];
} metadata_t;

//crisprs we want to compare are inserted into this
typedef struct {
    uint64_t id;
    uint64_t seq;
    uint64_t rev_seq;
} crispr_t;

class CrisprUtil {
    private:
        uint64_t * crisprs;
        metadata_t crispr_data;
        unsigned int max_offs;
        static const unsigned int max_mismatches = 4; //make this non static
        void _populate_cmap();
        void _find_off_targets(std::vector<crispr_t> queries);
    public:
        CrisprUtil();
        ~CrisprUtil() { if ( crispr_data.num_seqs > 0 ) delete[] crisprs; }
        uint8_t cmap[256];
        void load_binary(const std::string & outfile);
        std::string get_crispr(uint64_t id);
        uint64_t get_crispr_int(uint64_t id);
        uint64_t num_seqs();
        uint64_t seq_length();
        void find_off_targets(std::vector<uint64_t> ids);
        void find_off_targets(uint64_t start, uint64_t amount);
        void text_to_binary(const std::vector<std::string> & infiles, const std::string & outfile, metadata_t * data);
        void analyse_binary(const std::string & filename);
        /* void load_metadata(std::ifstream & text); */
        /* void search_by_seq(const std::string & filename, std::string seq, short pam_right); */
};

#endif