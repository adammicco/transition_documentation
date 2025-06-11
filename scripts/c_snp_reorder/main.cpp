//
//  c_snp_reorder.cpp
//  c_eig_reorder
//
//  Created by Adam Micco on 11/9/21.
//
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <map>
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

std::vector<std::string> read_geno(std::string geno_path) {
    std::ifstream infile;
    infile.open(geno_path);
    std::string line;
    std::vector<std::string> geno_vect;
    
    if (!infile) {
        std::cerr << "Unable to open " << geno_path << "!\n";
        exit(1);
    }
    
    int snp_pos = 0;
    while (getline(infile, line)) {
        geno_vect.push_back(line);
        if (snp_pos % 10000 == 0) {
            std::cout << snp_pos << "\n";
        }
        snp_pos++;
    }
    infile.close();
    std::cout << "Geno file loaded.\n";
    return geno_vect;
}

void write_geno(std::vector<std::string>& geno, std::string output_path) {
    std::ofstream outfile;
    outfile.open(output_path + ".geno");
    int snp_count = geno.size();
    
    for (int snp_idx = 0; snp_idx < snp_count; snp_idx++) {
        outfile << geno[snp_idx];
        outfile << "\n";
    }
    outfile.close();
}

std::vector<std::string> read_snp_file(std::string snp_path) {
    std::ifstream infile;
    infile.open(snp_path);
    std::string line;
    std::vector<std::string> snp_vec;

    if (!infile) {
        std::cerr << "Unable to open " << snp_path << "!\n";
        exit(1);
    }
    
    while (getline(infile, line)) {
        if (line.size() > 0) {
            snp_vec.push_back(ltrim(line));
        }
    }
    infile.close();
    return snp_vec;
}

std::map<std::string, int> read_driver(std::string driver_path) {
    std::ifstream infile;
    infile.open(driver_path);
    std::string line;
    std::map<std::string, int> driver_map;
    
    int idx = 0;
    while (getline(infile, line)) {
        if (line.size() > 0) {
            driver_map.insert({line, idx});
            idx++;
        }
    }
    return driver_map;
}

void write_snp_file(std::vector<std::string> snp_data, std::string output_path) {
    std::ofstream outfile;
    outfile.open(output_path + ".snp");
    
    for (string& line : snp_data) {
        outfile << line << "\n";
    }
    outfile.close();
}

std::vector<std::string> snp_2_list(std::vector<std::string> snp_vec) {
    std::vector<std::string> list_vec;

    for (string& line : snp_vec) {
        list_vec.push_back(line.substr(0, line.find_first_of(" \t")));
    }
    return list_vec;
}

template< typename order_iterator, typename value_iterator_str >
void reorder_destructive(order_iterator order_begin, order_iterator order_end, value_iterator_str &v, value_iterator_str &l)  {
    // typedef typename std::iterator_traits< value_iterator >::value_type value_t;
    typedef typename std::iterator_traits< order_iterator >::value_type index_t;
    typedef typename std::iterator_traits< order_iterator >::difference_type diff_t;
    
    diff_t remaining = order_end - 1 - order_begin;
    for ( index_t s = index_t(); remaining > 0; ++ s ) {
        index_t d = order_begin[s];
        if ( d == (diff_t) -1 ) continue;
        -- remaining;
        std::string tempv = v[s];
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

int main (int argc, char* argv []) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " INPUT_STEM SNP_DIRVER_FILEPATH OUTPUT_STEM" << std::endl;
        return 1;
    }
//
//        std::string geno_in_path = "/Users/adm515/eig_reorder/v50.geno";
//        std::string ind_in_path = "/Users/adm515/eig_reorder/v50.ind";
//        std::string snp_in_path = "/Users/adm515/eig_reorder/v50.snp";
//        std::string driver_path = "/Users/adm515/eig_reorder/v50_snp_driver";
//        std::string output_path = "/Users/adm515/eig_reorder/v50_snp_sorted";
    filesystem::path geno_in_rel = argv[1] + std::string(".geno");
    filesystem::path ind_in_rel = argv[1] + std::string(".ind");
    filesystem::path snp_in_rel = argv[1] + std::string(".snp");
    filesystem::path driver_rel = argv[2];
    std::string geno_in_path = filesystem::absolute(geno_in_rel);
    std::string ind_in_path = filesystem::absolute(ind_in_rel);
    std::string snp_in_path = filesystem::absolute(snp_in_rel);
    std::string driver_path = filesystem::absolute(driver_rel);
    std::string output_path = argv[3];
    
    std::vector<std::string> geno = read_geno(geno_in_path);
    std::vector<std::string> snp_in_vec = read_snp_file(snp_in_path);
    std::vector<std::string> snp_in_list_vec = snp_2_list(snp_in_vec);
    std::map<std::string, int> driver_map = read_driver(driver_path);
    std::cout << "All file dependancies loaded.\n";
    
    std::vector<int> snpIdx_vec;
    for (string& snpID : snp_in_list_vec) {
        // get index of the snpID in the driver
        snpIdx_vec.push_back(driver_map.at(snpID));
    }
        
//        auto pos = std::find(driver_vec.begin(), driver_vec.end(), snpID);
//
//        if(pos != driver_vec.end()) {
//            int snpIdx = pos - driver_vec.begin();
////          std::cout << snpIdx << "\n";
//            snpIdx_vec[idx] = snpIdx;
//            idx++;
//            if (idx % 1000 == 0) {
//                std::cout << idx << "\n";
//            }
//        }
//        else {
//            exit(1);
//        }
//    }
    
    reorder_destructive(snpIdx_vec.begin(), snpIdx_vec.end(), geno, snp_in_vec);
    write_geno(geno, output_path);
    write_snp_file(snp_in_vec, output_path);
    filesystem::copy_file(ind_in_path, output_path + ".ind");
    return 0;
}
