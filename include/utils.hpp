#pragma once

#include <string>
#include <vector>

struct PathParts {
    std::string directory;  // parent directory
    std::string basename;   // filename with extension
    std::string stem;       // filename without extension
    std::string extension;  // extension (with leading dot, if any)
};

void print_help();
void print_version(const std::string &version, const std::string &date);
std::vector<std::string> split_by_spaces(const std::string &);
std::vector<std::string> split_by_commas(const std::string &);
std::vector<int> string_to_int(const std::vector<std::string> &);
std::vector<int> list_positions(const std::vector<std::string> &, const std::vector<std::string> &);
PathParts splitPath(const std::string& );
std::vector<std::string> parse_extent(const std::string &);
bool directory_exists(const std::string&);
bool stats(const std::vector<int16_t> &, int16_t &, int16_t &, int16_t &, int16_t &, int16_t &, bool = false);
bool more_stats(const std::vector<int16_t> &, int16_t &, int16_t &, int16_t &, int16_t &, int16_t &, int16_t &, double &, bool = false);
bool ends_with(const std::string &, const std::string &);
std::vector<std::string> list_files(const std::string &, const std::string & = ".hdf");
std::vector<std::string> split_list(const std::string& str);
std::vector<int16_t> split_list_int(const std::string& str);
std::string extract_date(const std::string &);
std::string doy2date(const int&, const int&);
void extract_doy(const std::string&, int&, int&);
int get_position(const std::vector<std::string>&, const std::string&);
bool find_pattern(const std::string&, const std::string&);

bool files2ts(const std::vector<std::string> &list_files_1,
              const std::vector<std::string> &list_files_2,
              const std::vector<std::string> &list_input_dirs,
              const std::vector<std::string> & datasets,
              const std::string &output_dir,
              const std::string &output_file_name,
              const double &vmin,
              const double &vmax,
              const bool &slack,
              const int16_t &no_data
);
std::string match_file_doy(const int &year, const int &doy,
                           const std::vector<std::string> &list_files
);
int find_int_position(const std::vector<int>& vec, int value);
void extract_year(const std::string& filename, int& year);