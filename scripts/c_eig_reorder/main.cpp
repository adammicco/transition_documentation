#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <filesystem>
#include <sys/stat.h>
using namespace std;

const std::string WHITESPACE = " \n\r\t\f\v";

std::string ltrim(const std::string &s)
{
    size_t start = s.find_first_not_of(WHITESPACE);
    return (start == std::string::npos) ? "" : s.substr(start);
}

std::vector<std::vector<char> > read_geno_2_transpose(std::string geno_path) {
    std::ifstream infile;
    infile.open(geno_path);
    std::string line;
    std::vector<std::vector<char> > geno_vect;

    if (!infile) {
        std::cerr << "Unable to open " << geno_path << "!\n";
        exit(1);
    }
    // Use first row of input geno file to initialize the inner vectors of the 2D geno_vec and populate the first column of the transpose.
    getline(infile, line);
    for (char &ch : line) {
        vector<char> v;
        v.push_back(ch);
        geno_vect.push_back(v);
    }
    std::cout << "Transpose array contstructed at first position.\n";
    // Transpose each following row in the input geno file.
    int snp_pos = 1;
    while (getline(infile, line)) {
        int idx = 0;
        for (char &ch : line) {
            geno_vect[idx].push_back(ch);
            idx++;
        }
        if (snp_pos % 10000 == 0) {
            std::cout << snp_pos << "\n";
        }
        snp_pos++;
    }
    infile.close();
    std::cout << "Transposed geno loaded.\n";
    return geno_vect;
}

void write_untransposed_geno(std::vector<std::vector<char> >& transposed_geno, std::string output_path) {
    std::ofstream outfile;
    outfile.open(output_path + ".geno");
    int snp_count = transposed_geno[0].size();
    int ind_count = transposed_geno.size();
    
    for (int snp_idx = 0; snp_idx < snp_count; snp_idx++) {
        for (int ind_idx = 0; ind_idx < ind_count; ind_idx++) {
            outfile << transposed_geno[ind_idx][snp_idx];
        }
        outfile << "\n";
    }
    outfile.close();
}

std::vector<std::string> read_ind_file(std::string ind_path) {
    std::ifstream infile;
    infile.open(ind_path);
    std::string line;
    std::vector<std::string> ind_vec;

    if (!infile) {
        std::cerr << "Unable to open " << ind_path << "!\n";
        exit(1);
    }
    // Read ind file and generate a vector like: {"SAMPLE_ID\tSEX\tPOP_ID", ...} where indices should match row indices in the transposed geno data matrix
    while (getline(infile, line)) {
        if (line.size() > 0) {
            ind_vec.push_back(ltrim(line));
        }
    }
    infile.close();
    return ind_vec;
}

void write_ind_file(std::vector<std::string> ind_data, std::string output_path) {
    std::ofstream outfile;
    outfile.open(output_path + ".ind");
    
    for (string& line : ind_data) {
        outfile << line << "\n";
    }
    outfile.close();
}

std::vector<std::string> ind_2_list(std::vector<std::string> ind_vec) {
    std::vector<std::string> list_vec;

    for (string& line : ind_vec) {
        list_vec.push_back(line.substr(0, line.find_first_of(" \t")));
    }
    return list_vec;
}

template< typename order_iterator, typename value_iterator_vec, typename value_iterator_str >
void reorder_destructive(order_iterator order_begin, order_iterator order_end, value_iterator_vec &v, value_iterator_str &l)  {
    // typedef typename std::iterator_traits< value_iterator >::value_type value_t;
    typedef typename std::iterator_traits< order_iterator >::value_type index_t;
    typedef typename std::iterator_traits< order_iterator >::difference_type diff_t;
    
    diff_t remaining = order_end - 1 - order_begin;
    for ( index_t s = index_t(); remaining > 0; ++ s ) {
        index_t d = order_begin[s];
        if ( d == (diff_t) -1 ) continue;
        -- remaining;
        std::vector<char> tempv = v[s];
        std::string templ = l[s];
        for ( index_t d2; d != s; d = d2 ) {
            swap( tempv, v[d] );
            swap( templ, l[d] );
            swap( order_begin[d], d2 = (diff_t) -1 );
            -- remaining;
        }
        v[s] = tempv;
        l[s] = templ;
    }
}

int main (int argc, char* argv[])
{
     if (argc < 4) {
         std::cerr << "Usage: " << argv[0] << " INPUT_STEM DIRVER_FILEPATH OUTPUT_STEM" << std::endl;
         return 1;
     }

    // TODO Resolve file inputs before calling fuctions below
//    std::string geno_in_path = "/Users/adm515/eig_reorder/v50.geno";
//    std::string ind_in_path = "/Users/adm515/eig_reorder/v50.ind";
//    std::string snp_in_path = "/Users/adm515/eig_reorder/v50.snp";
//    std::string driver_path = "/Users/adm515/eig_reorder/driver";
//    std::string output_path = "/Users/adm515/eig_reorder/v50_sorted";
    filesystem::path geno_in_rel = argv[1] + std::string(".geno");
    filesystem::path ind_in_rel = argv[1] + std::string(".ind");
    filesystem::path snp_in_rel = argv[1] + std::string(".snp");
    filesystem::path driver_rel = argv[2];
    std::string geno_in_path = filesystem::absolute(geno_in_rel);
    std::string ind_in_path = filesystem::absolute(ind_in_rel);
    std::string snp_in_path = filesystem::absolute(snp_in_rel);
    std::string driver_path = filesystem::absolute(driver_rel);
    std::string output_path = argv[3];

    std::vector<std::vector<char> > transposed_geno = read_geno_2_transpose(geno_in_path);
    std::vector<std::string> ind_in_vec = read_ind_file(ind_in_path);
    std::vector<std::string> ind_in_list_vec = ind_2_list(ind_in_vec);
    std::vector<std::string> driver_vec = read_ind_file(driver_path);
    std::cout << "All file dependancies loaded.\n";

    std::vector<int> indIdx_vec;
    for (string& indID : ind_in_list_vec) {
        // get index of the indID in the driver
        auto pos = std::find(driver_vec.begin(), driver_vec.end(), indID);

        if(pos != driver_vec.end()) {
            int indIdx = pos - driver_vec.begin();
            std::cout << indIdx << "\n";
            indIdx_vec.push_back(indIdx);
        }
        else {
            exit(1);
        }
        // Now that we know the index of the indID in the driver file and the row of the transposed geno file, we can take those and push them to the back of their respective vectors. Once we've interated through the whole driver_vec, we'll have put these in the proper order.
    }
    reorder_destructive(indIdx_vec.begin(), indIdx_vec.end(), transposed_geno, ind_in_vec);
    write_untransposed_geno(transposed_geno, output_path);
    write_ind_file(ind_in_vec, output_path);
    filesystem::copy_file(snp_in_path, output_path + ".snp");
    
    return 0;
}
