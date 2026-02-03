#include <cstdlib>
#include <iostream>
#include <cstdint>
#include <limits>
#include <cmath>
#include <sys/stat.h>
#include <dirent.h>
#include <cstring>
#include <sstream>
#include <regex>
#include <ctime>
#include <iomanip>
#include <algorithm>
#include <fstream>

#include "utils.hpp"

#include "trans_coords.hpp"
#include "hdf4_utils.hpp"

void print_help() {
    std::cout << "Usage: tss [OPTIONS]\n"
              << "Extraction of time series from HDF files in multiple dirs.\n"
              << "Options:\n"
              << "  -f, --format FORMAT  Specify HDF format: hdf4 or hdf5 (default: hdf4)\n"
              << "  -i, --idir DIR       Input directories to scan HDF files (Landsat,MODIS)\n"
              << "  -o, --odir DIR       Output directory to write CSV time series files\n"
              << "  -d, --datasets LIST  List of datasets to read (e.g. \"NDVI,EVI (calc)\" or 0,1,2)\n"
              << "  -p, --pattern STR    Pattern to filter files to process (e.g. '.A2008')\n"
              << "  -n, --name STR       Basename for the CSV time series files\n"
              << "  -b, --bounds LIST    Row,col to start and row,col to finish extraction (e.g. 0,0,100,100)\n"
              << "  -s, --slack          Matches dates to -1 and +1 days of DOY\n"
              << "  -l, --lower VALUE    Exclude values less than VALUE from extraction (default: -10000)\n"
              << "  -u, --upper VALUE    Exclude values greater than VALUE from extraction (default: 10000)\n"
              << "  -v, --version        Show version information\n"
              << "  -h, --help           Show this help message\n"
              << std::endl;
}

void print_version(const std::string &version, const std::string &date) {
    std::cout << "tss - extraction of time series from HDF files in multiple dirs v"
              << version << " release " << date << "\n"
              << "\nCopyright (C) 2025 Eduardo Jimenez Hernandez <eduardojh@arizona.edu>.\n"
              << "License GPLv3+: GNU GPL version 3 or later <https://gnu.org/licenses/gpl.html>.\n"
              << "This is free software: you are free to change and redistribute it.\n"
              << "There is NO WARRANTY, to the extent permitted by law."
              << std::endl;
}

std::vector<std::string> split_by_spaces(const std::string &input) {
    std::vector<std::string> result;
    std::string word;

    for (size_t i = 0; i < input.size(); ++i) {
        char c = input[i];
        if (c == ',') {
            if (!word.empty()) {
                result.push_back(word);
                word.clear();
            }
        } else {
            word += c;
        }
    }
    if (!word.empty()) {
        result.push_back(word);
    }

    std::cout << "String split into " << result.size() << " elements." << std::endl;

    return result;
}

std::vector<std::string> parse_extent(const std::string &input) {
    std::string s = input;

    // Check and remove surrounding brackets
    if (!s.empty() && s.front() == '[' && s.back() == ']') {
        s = s.substr(1, s.size() - 2);
    }

    return split_by_commas(s);
}

std::vector<std::string> split_by_commas(const std::string &input) {
    std::vector<std::string> result;
    std::string token;

    for (size_t i = 0; i < input.size(); ++i) {
        char c = input[i];
        if (c == ',') {
            if (!token.empty()) {
                result.push_back(token);
                // std::cout << "word:" << token << "\n";
                token.clear();
            }
        } else {
            token += c;
        }
    }
    if (!token.empty()) {
        result.push_back(token);
    }

    // std::cout << "String split into " << result.size() << " elements." << std::endl;

    return result;
}

std::vector<int> string_to_int(const std::vector<std::string> &strings) {
    std::vector<int> numbers;
    numbers.reserve(strings.size());

    for (const auto &s : strings) {
        try {
            size_t pos;
            int value = std::stoi(s, &pos);
            if (pos != s.size()) {
                return {}; // invalid conversion
            }
            numbers.push_back(value);
        } catch (...) {
            return {}; // conversion failed
        }
    }
    
    return numbers;
}

std::vector<int> list_positions(const std::vector<std::string> &list1,
                                const std::vector<std::string> &list2) {
    // list1 should be larger than list2
    std::vector<int> list_pos = string_to_int(list2);

    if (list_pos.empty()) {

        // Conversion to integer failed, thus convert names to positions
        std::cout << "Requested dataset names\n";
        for (size_t i = 0; i < list2.size(); ++i) {
            bool found = false;
            for (size_t j = 0; j < list1.size(); ++j) {
                if (list2[i] == list1[j]) {
                    found = true;
                    list_pos.push_back(j);
                    break;
                }
            }
            std::cout << "  " << list2[i] << " found: " << (found ? "true" : "false") << "\n";
        }
        
        std::cout << "Positions found: " << list_pos.size() << std::endl;

    } else {

        // Check for values out of range
        std::vector<int> dataset_pos_temp(list_pos);
        list_pos.clear();

        std::cout << "Checking for values out or range, len1=" << dataset_pos_temp.size() << " len2=" << list_pos.size() << "\n";

        for (size_t i = 0; i < dataset_pos_temp.size(); ++i) {
            size_t p = dataset_pos_temp[i];
            // std::cout << p << "\n";
            if (p >= list1.size()) {
                std::cerr << "  ERROR: index " << p << " out of range (" << list1.size() << "). Discarding." << std::endl;
            } else {
                std::cout << "  Index " << p << " within range, saving." << std::endl;
                list_pos.push_back(p);
            }
        }
    }

    return list_pos;
}

PathParts splitPath(const std::string& path) {
    PathParts parts;

    // Find last path separator (works for Unix '/' and Windows '\\')
    size_t slashPos = path.find_last_of("/\\");
    if (slashPos == std::string::npos) {
        parts.basename = path;  // no directory
    } else {
        parts.directory = path.substr(0, slashPos);
        parts.basename = path.substr(slashPos + 1);
    }

    // Find extension in basename
    size_t dotPos = parts.basename.find_last_of('.');
    if (dotPos == std::string::npos) {
        parts.stem = parts.basename;  // no extension
    } else {
        parts.stem = parts.basename.substr(0, dotPos);
        parts.extension = parts.basename.substr(dotPos);
    }

    return parts;
}

bool directory_exists(const std::string& path) {
    struct stat info;
    if (stat(path.c_str(), &info) != 0) {
        return false; // cannot access path
    }
    return (info.st_mode & S_IFDIR) != 0;
}

bool stats(const std::vector<int16_t> &data,
           int16_t &minVal,
           int16_t &maxVal,
           int16_t &avg,
           int16_t &stddev,
           int16_t &count,
           bool useWelford)
{
    // Calculates min, max, average, standard deviation and count (number of pixels with data)
    if (data.empty()) {
        return false; // No data, cannot compute stats
    }

    count = static_cast<int16_t>(data.size());

    // Initialize min/max
    minVal = std::numeric_limits<int16_t>::max();
    maxVal = std::numeric_limits<int16_t>::min();

    if (!useWelford) {
        // --- Simple two-pass method ---
        int64_t sum = 0;
        for (size_t i = 0; i < data.size(); ++i) {
            if (data[i] < minVal) minVal = data[i];
            if (data[i] > maxVal) maxVal = data[i];
            sum += data[i];
        }

        double mean = static_cast<double>(sum) / data.size();

        double variance = 0.0;
        for (size_t i = 0; i < data.size(); ++i) {
            double diff = data[i] - mean;
            variance += diff * diff;
        }

        // // biased
        // variance /= data.size();

        // unbiased: divide by (N-1) if more than 1 element
        if (data.size() > 1) {
            variance /= (data.size() - 1);
        } else {
            variance = 0.0;
        }

        avg = static_cast<int16_t>(std::lround(mean));
        stddev = static_cast<int16_t>(std::lround(std::sqrt(variance)));
    } else {
        // --- Welford’s online algorithm ---
        double mean = 0.0;
        double M2 = 0.0;

        for (size_t i = 0; i < data.size(); ++i) {
            int16_t x = data[i];
            if (x < minVal) minVal = x;
            if (x > maxVal) maxVal = x;

            double delta = x - mean;
            mean += delta / (i + 1);
            M2 += delta * (x - mean);
        }

        double variance = M2 / data.size();

        avg = static_cast<int16_t>(std::lround(mean));
        stddev = static_cast<int16_t>(std::lround(std::sqrt(variance)));
    }

    return true;
}

bool more_stats(const std::vector<int16_t> &data,
           int16_t &minVal,
           int16_t &maxVal,
           int16_t &avg,
           int16_t &stddev,
           int16_t &count,
           int16_t &median,
           double &zscore,
           bool useWelford)
{
    // Calculates the same as stats plus z-score
    if (data.empty()) {
        return false; // No data, cannot compute stats
    }

    count = static_cast<int16_t>(data.size());

    // Initialize min/max
    minVal = std::numeric_limits<int16_t>::max();
    maxVal = std::numeric_limits<int16_t>::min();

    double mean;

    if (!useWelford) {
        // --- Simple two-pass method ---
        int64_t sum = 0;
        for (size_t i = 0; i < data.size(); ++i) {
            if (data[i] < minVal) minVal = data[i];
            if (data[i] > maxVal) maxVal = data[i];
            sum += data[i];
        }

        mean = static_cast<double>(sum) / data.size();

        double variance = 0.0;
        for (size_t i = 0; i < data.size(); ++i) {
            double diff = data[i] - mean;
            variance += diff * diff;
        }

        // // biased
        // variance /= data.size();

        // unbiased: divide by (N-1) if more than 1 element
        if (data.size() > 1) {
            variance /= (data.size() - 1);
        } else {
            variance = 0.0;
        }

        avg = static_cast<int16_t>(std::lround(mean));
        stddev = static_cast<int16_t>(std::lround(std::sqrt(variance)));
    } else {
        // --- Welford’s online algorithm ---
        mean = 0.0;
        double M2 = 0.0;

        for (size_t i = 0; i < data.size(); ++i) {
            int16_t x = data[i];
            if (x < minVal) minVal = x;
            if (x > maxVal) maxVal = x;

            double delta = x - mean;
            mean += delta / (i + 1);
            M2 += delta * (x - mean);
        }

        double variance = M2 / data.size();

        avg = static_cast<int16_t>(std::lround(mean));
        stddev = static_cast<int16_t>(std::lround(std::sqrt(variance)));
    }

    // Estimate the Z-Score
    double cum_z = 0.0;
    for (size_t i = 0; i < data.size(); ++i) {
        cum_z += (static_cast<double>(data[i]) - mean) / static_cast<double>(stddev);
    }
    // Average the Z-Score of all pixels
    // zscore = static_cast<int16_t>(std::lround(cum_z / data.size()));
    zscore = cum_z / data.size();

    // Calculate median, original data will be sorted
    size_t mid = count / 2;
    std::vector<int16_t> values = data;
    std::sort(values.begin(), values.end());
    if (count % 2 == 0) {
        // Even: average two middle values
        median = (values[mid - 1] + values[mid]) / 2.0;
    } else {
        median = values[mid];
    }

    return true;
}

bool ends_with(const std::string &str, const std::string &suffix) {
    if (suffix.size() > str.size()) return false;
    return std::equal(suffix.rbegin(), suffix.rend(), str.rbegin());
}

std::vector<std::string> list_files(const std::string &directory,
                                   const std::string &extension)
{
    // Returns a list of files in the directory matching the extension
    std::vector<std::string> files;

    DIR *dir = opendir(directory.c_str());
    if (!dir) {
        return files; // directory not found or not accessible
    }

    struct dirent *entry;
    while ((entry = readdir(dir)) != nullptr) {
        std::string name(entry->d_name);

        // Skip . and ..
        if (name == "." || name == "..") continue;

        std::string fullPath = directory + "/" + name;

        // Check if it's a file (not a directory)
        struct stat s;
        if (stat(fullPath.c_str(), &s) == 0 && S_ISREG(s.st_mode)) {
            if (extension.empty() || ends_with(name, extension)) {
                files.push_back(name); // store filename only
            }
        }
    }

    closedir(dir);
    return files;
}

std::vector<std::string> split_list(const std::string& str) {
    // Splits a list of comma separated strings into a vector
    std::vector<std::string> result;
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, ',')) {
        if (!token.empty()) result.push_back(token);
    }
    return result;
}

std::vector<int16_t> split_list_int(const std::string& str) {
    // Splits a list of comma separated strings into a vector of int16
    std::vector<int16_t> result;
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, ',')) {
        if (!token.empty()) result.push_back(std::stoi(token));
    }
    return result;
}

std::string extract_date(const std::string& filename) {
    // Regex: look for ".A" followed by 7 digits (YYYYDDD)
    std::regex pattern("\\.A(\\d{4})(\\d{3})\\.");
    std::smatch match;

    if (!std::regex_search(filename, match, pattern)) {
        // If no match found, return empty string or throw
        return "";  
    }

    int year = std::stoi(match[1].str());
    int dayOfYear = std::stoi(match[2].str());

    // Convert year + day-of-year to YYYY-MM-DD
    std::tm time_in = {};
    time_in.tm_year = year - 1900;
    time_in.tm_mday = dayOfYear;

    std::mktime(&time_in);  // normalize to calendar date

    std::ostringstream oss;
    oss << std::put_time(&time_in, "%Y-%m-%d");
    return oss.str();
}

std::string doy2date(const int& year, const int& doy) {
    // Convert year + day-of-year to YYYY-MM-DD
    std::tm time_in = {};
    time_in.tm_year = year - 1900;
    time_in.tm_mday = doy;

    std::mktime(&time_in);  // normalize to calendar date

    std::ostringstream oss;
    oss << std::put_time(&time_in, "%Y-%m-%d");
    return oss.str();
}

void extract_doy(const std::string& filename, int& year, int& dayOfYear) {
    // Regex: look for ".A" followed by 7 digits (YYYYDDD)
    std::regex pattern("\\.A(\\d{4})(\\d{3})\\.");
    std::smatch match;

    if (!std::regex_search(filename, match, pattern)) {
        // If no match found, return empty string or throw
        return;  
    }

    year = std::stoi(match[1].str());
    dayOfYear = std::stoi(match[2].str());
}

int get_position(const std::vector<std::string>& list, const std::string& target) {
    for (size_t i = 0; i < list.size(); ++i) {
        if (list[i] == target) {
            return static_cast<int>(i);
        }
    }
    return -1; // Not found
}

bool find_pattern(const std::string& filename, const std::string& pattern) {
    // Returns true if 'pattern' is found anywhere inside 'filename'.
    return filename.find(pattern) != std::string::npos;
}

bool files2ts(const std::vector<std::string> &list_files_1,
              const std::vector<std::string> &list_files_2,
              const std::vector<std::string> &list_input_dirs,
              const std::vector<std::string> &datasets,
              const std::string &output_dir,
              const std::string &output_file_name,
              const double &vmin,
              const double &vmax,
              const bool &slack,
              const int16_t &no_data
) {
    // Create time series from list of files
    // Landsat=list_files_1, MODIS=list_files_2
    
    const std::vector<std::string> stats = {"avg", "min", "max", "stdev", "count", "percent", "sza"};
    const std::vector<int> slack_days = {0, -1, 1};
    const std::string EXTENSION = ".csv";

    // Initialize CSV file to save time series
    std::ofstream csvFile;
    std::string stats_file1 = output_dir + output_file_name + "_" + datasets[0] + EXTENSION;
    if (output_dir.back() != '/') {
        stats_file1 = output_dir + "/" + output_file_name + "_" + datasets[0] + EXTENSION;
    }
    std::cout << "Creating: " << stats_file1 << "\n";
    csvFile.open(stats_file1);
    if (!csvFile.is_open()) {
        std::cerr << "ERROR: Failed to open stats CSV file" << std::endl;
        return EXIT_FAILURE;
    }
    // CSV headers
    csvFile << "Date,Year,DOY,Pixel,Row,Col,L_Min,L_Max,L_Avg,L_Stdev,L_Count,L_Percent";
    if (list_input_dirs.size() > 1 ) {
        csvFile << ",MODIS_DOY,M_Min,M_Max,M_Avg,M_Stdev,M_Count,M_Percent,M_SZA";
    }
    csvFile << std::endl;

    // Create a second CSV if second dataset provided
    std::ofstream csvFile2;
    std::string stats_file2;
    if (datasets.size() > 1) {
        stats_file2 = output_dir + output_file_name + "_" + datasets[1] + EXTENSION;
        if (output_dir.back() != '/') {
            stats_file2 = output_dir + "/" + output_file_name + "_" + datasets[1] + EXTENSION;
        }
        std::cout << "Creating: " << stats_file2 << "\n";
        csvFile2.open(stats_file2);
        if (!csvFile2.is_open()) {
            std::cerr << "ERROR: Failed to open stats CSV file" << std::endl;
            return EXIT_FAILURE;
        }
        csvFile2 << "Date,Year,DOY,Pixel,Row,Col,L_Min,L_Max,L_Avg,L_Stdev,L_Count,L_Percent";
        if (list_input_dirs.size() > 1 ) {
            csvFile2 << ",MODIS_DOY,M_Min,M_Max,M_Avg,M_Stdev,M_Count,M_Percent,M_SZA";
        }
        csvFile2 << std::endl;
    }

    // Create a third CSV if third dataset provided
    std::ofstream csvFile3;
    std::string stats_file3;
    if (datasets.size() > 2) {
        stats_file3 = output_dir + output_file_name + "_" + datasets[2] + EXTENSION;
        if (output_dir.back() != '/') {
            stats_file3 = output_dir + "/" + output_file_name + "_" + datasets[2] + EXTENSION;
        }
        std::cout << "Creating: " << stats_file3 << "\n";
        csvFile3.open(stats_file3);
            if (!csvFile3.is_open()) {
            std::cerr << "ERRIR: Failed to open stats CSV file" << std::endl;
            return EXIT_FAILURE;
        }
        csvFile3 << "Date,Year,DOY,Pixel,Row,Col,L_Min,L_Max,L_Avg,L_Stdev,L_Count,L_Percent";
        if (list_input_dirs.size() > 1 ) {
            csvFile3 << ",MODIS_DOY,M_Min,M_Max,M_Avg,M_Stdev,M_Count,M_Percent,M_SZA";
        }
        csvFile3 << std::endl;
    }

    // Iterate over Landsat files in list_files_1
    for (size_t i=0; i < list_files_1.size(); ++i) {
        int year;
        int doy;
        std::string date;
        std::string file_name;
        std::string abs_file_name;
        std::vector<std::string> dataset_names;
        std::vector<std::string> modis_dataset_names;
        std::vector<std::string> match_file_names;
        std::vector<int> match_day;

        std::vector<int16_t> landsat_dataset_avg;
        std::vector<int16_t> landsat_dataset_min;
        std::vector<int16_t> landsat_dataset_max;
        std::vector<int16_t> landsat_dataset_std;
        std::vector<int16_t> landsat_dataset_cnt;
        std::vector<double> landsat_dataset_per;

        std::vector<int16_t> modis_dataset_avg;
        std::vector<int16_t> modis_dataset_min;
        std::vector<int16_t> modis_dataset_max;
        std::vector<int16_t> modis_dataset_std;
        std::vector<int16_t> modis_dataset_cnt;
        std::vector<double> modis_dataset_per;
        std::vector<int16_t> modis_dataset_sza;

        // In case there is a second DOY (-1: the day before)
        std::vector<int16_t> modis_dataset_avg_2;
        std::vector<int16_t> modis_dataset_min_2;
        std::vector<int16_t> modis_dataset_max_2;
        std::vector<int16_t> modis_dataset_std_2;
        std::vector<int16_t> modis_dataset_cnt_2;
        std::vector<double> modis_dataset_per_2;
        std::vector<int16_t> modis_dataset_sza_2;

        // In case there is a third DOY (+1: the day after)
        std::vector<int16_t> modis_dataset_avg_3;
        std::vector<int16_t> modis_dataset_min_3;
        std::vector<int16_t> modis_dataset_max_3;
        std::vector<int16_t> modis_dataset_std_3;
        std::vector<int16_t> modis_dataset_cnt_3;
        std::vector<double> modis_dataset_per_3;
        std::vector<int16_t> modis_dataset_sza_3;
        
        // Get the current Landsat file name
        file_name = list_files_1[i];
        
        // Prepare file name and extract string date, year, and DOY
        date = extract_date(file_name);
        extract_doy(file_name, year, doy);

        std::cout << std::endl;
        std::cout << "============================================================================================\n";
        std::cout << "=== Landsat File (" << i+1 << "): " << file_name << " " << date << " ===\n";
        std::cout << "============================================================================================\n";

        // Split each file name
        PathParts file_parts = splitPath(file_name);
        std::cout << "  Directory: " << file_parts.directory << "\n";
        std::cout << "  Basename:  " << file_parts.basename << "\n";
        std::cout << "  Stem:      " << file_parts.stem << "\n";
        std::cout << "  Extension: " << file_parts.extension << "\n";

        std::cout << "Reading datasets from HDF file...\n";
        // Create file name with absolute path
        abs_file_name = list_input_dirs[0] + file_name;
        if (list_input_dirs[0].back() != '/') {
            abs_file_name =list_input_dirs[0] + "/"  + file_name;
        }
        dataset_names = read_hdf4_info(abs_file_name); // Names in file
        // std::cout << "Datasets in file: \n"; 
        // for (size_t j=0; j < dataset_names.size(); ++j) {
        //     std::cout << " " << j << " " << dataset_names[j] << "\n";
        // }
        
        if (slack) {
            std::cout << "Slack=true, matching dates with +/- 1 day.\n";
        }

        // Look for matches one day before and after, if slack is true
        // 'slack_days' sets priority as 0=same day, -1=day before, 1=day after
        for (int day : slack_days) {
            int slack_doy = doy + day;

            std::cout << "Looking for matches with " << year << "-" << slack_doy << " ";
            // Skip invalid DOYs
            if (slack_doy < 1 || slack_doy > 366) {
                std::cout << "invalid DOY. Skipping.";
                continue;
            }
            std::cout << "\n";

            // Also skip when slack is false
            if ((slack && day == -1) || (slack && day == 1)) {
                continue;
            }

            // Match Landsat date to MODIS file name
            std::string match_file_name = match_file_doy(year, slack_doy, list_files_2);

            if (match_file_name.empty()) {
                std::cout << "No MODIS file matched the date: " << year << "-" << slack_doy << ", skipping...\n";
                continue;
            }

            // Valid match, save the found file and the day (-1,0,+1)
            match_file_names.push_back(match_file_name);
            match_day.push_back(day);

            std::cout << "============================================================================================\n";
            std::cout << "=== MODIS File: " << match_file_name << " " << year << "-" << slack_doy << " ===\n";
            std::cout << "============================================================================================\n";

        }

        // // Check the lengths match
        // if (match_file_names.size() != match_day.size()) {
        //     std::cout << "ERROR: Something went wrong match_file_names=" << match_file_names.size() << ", match_day=" << match_day.size() <<"\n";
        //     break;
        // }

        // If no matched files skip
        if (match_file_names.empty()) {
            std::cout << "No matching MODIS files for this Landsat file.\n";
            continue;
        }

        // Open the datasets
        int ds = 0;
        for (std::string dataset : datasets) {
            
            int dataset_rows;
            int dataset_cols;
            int used_slack_doy = 367; // -1,0,+1 MODIS day match

            for (std::string stat : stats) {
                // Create the dataset name with the stat
                std::string single_dataset = dataset + "_" + stat;

                // Find the dataset's position
                int pos = get_position(dataset_names, single_dataset);
                std::cout << "\n*************************************************************************************\n";
                if (stat == "sza") {
                    std::cout << "--IGNORING Landsat dataset: " << single_dataset << " pos=(" << pos << ")" << "\n";
                } else {
                    std::cout << "--Reading Landsat dataset: " << single_dataset << " pos=(" << pos << ")" << "\n";
                }

                // For the special case of "EVI2 (calc)"
                if (dataset == "EVI2" && pos == -1) {
                    // Try again
                    // dataset = "EVI2 (calc)";
                    single_dataset = "EVI2 (calc)_" + stat;
                    pos = get_position(dataset_names, single_dataset);
                    std::cout << "'EVI2' not found, trying '" << single_dataset << "'... \n";
                }

                // if (pos == -1) {
                //     std::cout << "WARNING! Dataset " << single_dataset << " pos=(" << pos << ") not found! Skipping.\n";
                //     continue;
                // }

                // *************** Read all Landsat datasets ***************

                if (stat == "avg") {
                    if (!read_hdf4_by_index(abs_file_name, pos, landsat_dataset_avg, dataset_rows, dataset_cols)) {
                        std::cerr << "ERROR: Cannot read HDF4 dataset " << pos << std::endl;
                        return EXIT_FAILURE;
                    }
                    std::cout << "--Dataset: " << single_dataset << " [" << dataset_rows << "x"
                        << dataset_cols << "] (" << landsat_dataset_avg.size() << ")\n";
                }

                if (stat == "min") {
                    if (!read_hdf4_by_index(abs_file_name, pos, landsat_dataset_min, dataset_rows, dataset_cols)) {
                        std::cerr << "ERROR: Cannot read HDF4 dataset " << pos << std::endl;
                        return EXIT_FAILURE;
                    }
                    std::cout << "--Dataset: " << single_dataset << " [" << dataset_rows << "x"
                        << dataset_cols << "] (" << landsat_dataset_min.size() << ")\n";
                }

                if (stat == "max") {
                    if (!read_hdf4_by_index(abs_file_name, pos, landsat_dataset_max, dataset_rows, dataset_cols)) {
                        std::cerr << "ERROR: Cannot read HDF4 dataset " << pos << std::endl;
                        return EXIT_FAILURE;
                    }
                    std::cout << "--Dataset: " << single_dataset << " [" << dataset_rows << "x"
                        << dataset_cols << "] (" << landsat_dataset_max.size() << ")\n";
                }

                if (stat == "stdev") {
                    if (!read_hdf4_by_index(abs_file_name, pos, landsat_dataset_std, dataset_rows, dataset_cols)) {
                        std::cerr << "ERROR: Cannot read HDF4 dataset " << pos << std::endl;
                        return EXIT_FAILURE;
                    }
                    std::cout << "--Dataset: " << single_dataset << " [" << dataset_rows << "x"
                        << dataset_cols << "] (" << landsat_dataset_std.size() << ")\n";
                }

                if (stat == "count") {
                    if (!read_hdf4_by_index(abs_file_name, pos, landsat_dataset_cnt, dataset_rows, dataset_cols)) {
                        std::cerr << "ERROR: Cannot read HDF4 dataset " << pos << std::endl;
                        return EXIT_FAILURE;
                    }
                    std::cout << "--Dataset: " << single_dataset << " [" << dataset_rows << "x"
                        << dataset_cols << "] (" << landsat_dataset_cnt.size() << ")\n";
                }

                if (stat == "percent") {
                    if (!read_hdf4_by_index(abs_file_name, pos, landsat_dataset_per, dataset_rows, dataset_cols)) {
                        std::cerr << "ERROR: Cannot read HDF4 dataset " << pos << std::endl;
                        return EXIT_FAILURE;
                    }
                    std::cout << "--Dataset: " << single_dataset << " [" << dataset_rows << "x"
                        << dataset_cols << "] (" << landsat_dataset_per.size() << ")\n";
                }

                // *************** Read all MODIS datasets for FIRST DOY ***************
                
                // first position is highest priority
                // used_slack_doy = match_day[0];
                // Create file name with absolute path
                std::string abs_modis_file_name = list_input_dirs[1] + match_file_names[0];
                if (list_input_dirs[1].back() != '/') {
                    abs_modis_file_name =list_input_dirs[1] + "/"  + match_file_names[0];
                }
                modis_dataset_names = read_hdf4_info(abs_modis_file_name);

                // Reset for MODIS in case of EVI2 (calc)
                if (single_dataset == "EVI2 (calc)_" + stat) {
                    single_dataset = dataset + "_" + stat;
                }

                // Special case for the Sensor Zenith Angle
                if (stat == "sza") {
                    single_dataset = "SensorZenith_avg";
                }

                // Find the dataset's position
                int modis_pos = get_position(modis_dataset_names, single_dataset);
                std::cout << "--Reading MODIS dataset: " << single_dataset << " pos=(" << modis_pos << ")" << "\n";

                if (modis_pos == -1) {
                    std::cout << "WARNING! Dataset " << single_dataset << " pos=(" << modis_pos << ") not found! Skipping.\n";
                    continue;
                }
                
                if (stat == "avg") {
                    if (!read_hdf4_by_index(abs_modis_file_name, modis_pos, modis_dataset_avg, dataset_rows, dataset_cols)) {
                        std::cerr << "ERROR: Cannot read HDF4 dataset " << modis_pos << std::endl;
                        return EXIT_FAILURE;
                    }
                    std::cout << "--Dataset: " << single_dataset << " [" << dataset_rows << "x"
                        << dataset_cols << "] (" << modis_dataset_avg.size() << ")\n";
                }

                if (stat == "min") {
                    if (!read_hdf4_by_index(abs_modis_file_name, modis_pos, modis_dataset_min, dataset_rows, dataset_cols)) {
                        std::cerr << "ERROR: Cannot read HDF4 dataset " << modis_pos << std::endl;
                        return EXIT_FAILURE;
                    }
                    std::cout << "--Dataset: " << single_dataset << " [" << dataset_rows << "x"
                        << dataset_cols << "] (" << modis_dataset_min.size() << ")\n";
                }

                if (stat == "max") {
                    if (!read_hdf4_by_index(abs_modis_file_name, modis_pos, modis_dataset_max, dataset_rows, dataset_cols)) {
                        std::cerr << "ERROR: Cannot read HDF4 dataset " << modis_pos << std::endl;
                        return EXIT_FAILURE;
                    }
                    std::cout << "--Dataset: " << single_dataset << " [" << dataset_rows << "x"
                        << dataset_cols << "] (" << modis_dataset_max.size() << ")\n";
                }

                if (stat == "stdev") {
                    if (!read_hdf4_by_index(abs_modis_file_name, modis_pos, modis_dataset_std, dataset_rows, dataset_cols)) {
                        std::cerr << "ERROR: Cannot read HDF4 dataset " << modis_pos << std::endl;
                        return EXIT_FAILURE;
                    }
                    std::cout << "--Dataset: " << single_dataset << " [" << dataset_rows << "x"
                        << dataset_cols << "] (" << modis_dataset_std.size() << ")\n";
                }

                if (stat == "count") {
                    if (!read_hdf4_by_index(abs_modis_file_name, modis_pos, modis_dataset_cnt, dataset_rows, dataset_cols)) {
                        std::cerr << "ERROR: Cannot read HDF4 dataset " << modis_pos << std::endl;
                        return EXIT_FAILURE;
                    }
                    std::cout << "--Dataset: " << single_dataset << " [" << dataset_rows << "x"
                        << dataset_cols << "] (" << modis_dataset_cnt.size() << ")\n";
                }

                if (stat == "percent") {
                    if (!read_hdf4_by_index(abs_modis_file_name, modis_pos, modis_dataset_per, dataset_rows, dataset_cols)) {
                        std::cerr << "ERROR: Cannot read HDF4 dataset " << modis_pos << std::endl;
                        return EXIT_FAILURE;
                    }
                    std::cout << "--Dataset: " << single_dataset << " [" << dataset_rows << "x"
                        << dataset_cols << "] (" << modis_dataset_per.size() << ")\n";
                }

                if (stat == "sza") {
                    if (!read_hdf4_by_index(abs_modis_file_name, modis_pos, modis_dataset_sza, dataset_rows, dataset_cols)) {
                        std::cerr << "ERROR: Cannot read HDF4 dataset " << modis_pos << std::endl;
                        return EXIT_FAILURE;
                    }
                    std::cout << "--Dataset: " << single_dataset << " [" << dataset_rows << "x"
                        << dataset_cols << "] (" << modis_dataset_sza.size() << ")\n";
                }

                // *************** Read all MODIS datasets for the SECOND DOY ***************

                if (match_day.size() > 1) {
                    // first position is highest priority
                    // used_slack_doy = match_day[1];
                    // Create file name with absolute path
                    std::string abs_modis_file_name = list_input_dirs[1] + match_file_names[1];
                    if (list_input_dirs[1].back() != '/') {
                        abs_modis_file_name =list_input_dirs[1] + "/"  + match_file_names[1];
                    }
                    modis_dataset_names = read_hdf4_info(abs_modis_file_name);

                    // Reset for MODIS in case of EVI2 (calc)
                    if (single_dataset == "EVI2 (calc)_" + stat) {
                        single_dataset = dataset + "_" + stat;
                    }

                    // Special case for the Sensor Zenith Angle
                    if (stat == "sza") {
                        single_dataset = "SensorZenith_avg";
                    }

                    // Find the dataset's position
                    int modis_pos = get_position(modis_dataset_names, single_dataset);
                    std::cout << "--Reading MODIS dataset: " << single_dataset << " pos=(" << modis_pos << ")" << "\n";

                    if (modis_pos == -1) {
                        std::cout << "WARNING! Dataset " << single_dataset << " pos=(" << modis_pos << ") not found! Skipping.\n";
                        continue;
                    }
                    
                    if (stat == "avg") {
                        if (!read_hdf4_by_index(abs_modis_file_name, modis_pos, modis_dataset_avg_2, dataset_rows, dataset_cols)) {
                            std::cerr << "ERROR: Cannot read HDF4 dataset " << modis_pos << std::endl;
                            return EXIT_FAILURE;
                        }
                        std::cout << "--Dataset: " << single_dataset << " [" << dataset_rows << "x"
                            << dataset_cols << "] (" << modis_dataset_avg_2.size() << ")\n";
                    }

                    if (stat == "min") {
                        if (!read_hdf4_by_index(abs_modis_file_name, modis_pos, modis_dataset_min_2, dataset_rows, dataset_cols)) {
                            std::cerr << "ERROR: Cannot read HDF4 dataset " << modis_pos << std::endl;
                            return EXIT_FAILURE;
                        }
                        std::cout << "--Dataset: " << single_dataset << " [" << dataset_rows << "x"
                            << dataset_cols << "] (" << modis_dataset_min_2.size() << ")\n";
                    }

                    if (stat == "max") {
                        if (!read_hdf4_by_index(abs_modis_file_name, modis_pos, modis_dataset_max_2, dataset_rows, dataset_cols)) {
                            std::cerr << "ERROR: Cannot read HDF4 dataset " << modis_pos << std::endl;
                            return EXIT_FAILURE;
                        }
                        std::cout << "--Dataset: " << single_dataset << " [" << dataset_rows << "x"
                            << dataset_cols << "] (" << modis_dataset_max_2.size() << ")\n";
                    }

                    if (stat == "stdev") {
                        if (!read_hdf4_by_index(abs_modis_file_name, modis_pos, modis_dataset_std_2, dataset_rows, dataset_cols)) {
                            std::cerr << "ERROR: Cannot read HDF4 dataset " << modis_pos << std::endl;
                            return EXIT_FAILURE;
                        }
                        std::cout << "--Dataset: " << single_dataset << " [" << dataset_rows << "x"
                            << dataset_cols << "] (" << modis_dataset_std_2.size() << ")\n";
                    }

                    if (stat == "count") {
                        if (!read_hdf4_by_index(abs_modis_file_name, modis_pos, modis_dataset_cnt_2, dataset_rows, dataset_cols)) {
                            std::cerr << "ERROR: Cannot read HDF4 dataset " << modis_pos << std::endl;
                            return EXIT_FAILURE;
                        }
                        std::cout << "--Dataset: " << single_dataset << " [" << dataset_rows << "x"
                            << dataset_cols << "] (" << modis_dataset_cnt_2.size() << ")\n";
                    }

                    if (stat == "percent") {
                        if (!read_hdf4_by_index(abs_modis_file_name, modis_pos, modis_dataset_per_2, dataset_rows, dataset_cols)) {
                            std::cerr << "ERROR: Cannot read HDF4 dataset " << modis_pos << std::endl;
                            return EXIT_FAILURE;
                        }
                        std::cout << "--Dataset: " << single_dataset << " [" << dataset_rows << "x"
                            << dataset_cols << "] (" << modis_dataset_per_2.size() << ")\n";
                    }

                    if (stat == "sza") {
                        if (!read_hdf4_by_index(abs_modis_file_name, modis_pos, modis_dataset_sza_2, dataset_rows, dataset_cols)) {
                            std::cerr << "ERROR: Cannot read HDF4 dataset " << modis_pos << std::endl;
                            return EXIT_FAILURE;
                        }
                        std::cout << "--Dataset: " << single_dataset << " [" << dataset_rows << "x"
                            << dataset_cols << "] (" << modis_dataset_sza_2.size() << ")\n";
                    }
                }

                // *************** Read all MODIS datasets for the THIRD DOY ***************

                if (match_day.size() > 2) {
                    // first position is highest priority
                    // used_slack_doy = match_day[2];
                    // Create file name with absolute path
                    std::string abs_modis_file_name = list_input_dirs[1] + match_file_names[2];
                    if (list_input_dirs[1].back() != '/') {
                        abs_modis_file_name =list_input_dirs[1] + "/"  + match_file_names[2];
                    }
                    modis_dataset_names = read_hdf4_info(abs_modis_file_name);

                    // Reset for MODIS in case of EVI2 (calc)
                    if (single_dataset == "EVI2 (calc)_" + stat) {
                        single_dataset = dataset + "_" + stat;
                    }

                    // Special case for the Sensor Zenith Angle
                    if (stat == "sza") {
                        single_dataset = "SensorZenith_avg";
                    }

                    // Find the dataset's position
                    int modis_pos = get_position(modis_dataset_names, single_dataset);
                    std::cout << "--Reading MODIS dataset: " << single_dataset << " pos=(" << modis_pos << ")" << "\n";

                    if (modis_pos == -1) {
                        std::cout << "WARNING! Dataset " << single_dataset << " pos=(" << modis_pos << ") not found! Skipping.\n";
                        continue;
                    }
                    
                    if (stat == "avg") {
                        if (!read_hdf4_by_index(abs_modis_file_name, modis_pos, modis_dataset_avg_3, dataset_rows, dataset_cols)) {
                            std::cerr << "ERROR: Cannot read HDF4 dataset " << modis_pos << std::endl;
                            return EXIT_FAILURE;
                        }
                        std::cout << "--Dataset: " << single_dataset << " [" << dataset_rows << "x"
                            << dataset_cols << "] (" << modis_dataset_avg_3.size() << ")\n";
                    }

                    if (stat == "min") {
                        if (!read_hdf4_by_index(abs_modis_file_name, modis_pos, modis_dataset_min_3, dataset_rows, dataset_cols)) {
                            std::cerr << "ERROR: Cannot read HDF4 dataset " << modis_pos << std::endl;
                            return EXIT_FAILURE;
                        }
                        std::cout << "--Dataset: " << single_dataset << " [" << dataset_rows << "x"
                            << dataset_cols << "] (" << modis_dataset_min_3.size() << ")\n";
                    }

                    if (stat == "max") {
                        if (!read_hdf4_by_index(abs_modis_file_name, modis_pos, modis_dataset_max_3, dataset_rows, dataset_cols)) {
                            std::cerr << "ERROR: Cannot read HDF4 dataset " << modis_pos << std::endl;
                            return EXIT_FAILURE;
                        }
                        std::cout << "--Dataset: " << single_dataset << " [" << dataset_rows << "x"
                            << dataset_cols << "] (" << modis_dataset_max_3.size() << ")\n";
                    }

                    if (stat == "stdev") {
                        if (!read_hdf4_by_index(abs_modis_file_name, modis_pos, modis_dataset_std_3, dataset_rows, dataset_cols)) {
                            std::cerr << "ERROR: Cannot read HDF4 dataset " << modis_pos << std::endl;
                            return EXIT_FAILURE;
                        }
                        std::cout << "--Dataset: " << single_dataset << " [" << dataset_rows << "x"
                            << dataset_cols << "] (" << modis_dataset_std_3.size() << ")\n";
                    }

                    if (stat == "count") {
                        if (!read_hdf4_by_index(abs_modis_file_name, modis_pos, modis_dataset_cnt_3, dataset_rows, dataset_cols)) {
                            std::cerr << "ERROR: Cannot read HDF4 dataset " << modis_pos << std::endl;
                            return EXIT_FAILURE;
                        }
                        std::cout << "--Dataset: " << single_dataset << " [" << dataset_rows << "x"
                            << dataset_cols << "] (" << modis_dataset_cnt_3.size() << ")\n";
                    }

                    if (stat == "percent") {
                        if (!read_hdf4_by_index(abs_modis_file_name, modis_pos, modis_dataset_per_3, dataset_rows, dataset_cols)) {
                            std::cerr << "ERROR: Cannot read HDF4 dataset " << modis_pos << std::endl;
                            return EXIT_FAILURE;
                        }
                        std::cout << "--Dataset: " << single_dataset << " [" << dataset_rows << "x"
                            << dataset_cols << "] (" << modis_dataset_per_3.size() << ")\n";
                    }

                    if (stat == "sza") {
                        if (!read_hdf4_by_index(abs_modis_file_name, modis_pos, modis_dataset_sza_3, dataset_rows, dataset_cols)) {
                            std::cerr << "ERROR: Cannot read HDF4 dataset " << modis_pos << std::endl;
                            return EXIT_FAILURE;
                        }
                        std::cout << "--Dataset: " << single_dataset << " [" << dataset_rows << "x"
                            << dataset_cols << "] (" << modis_dataset_sza_3.size() << ")\n";
                    }
                }

                // *************** Done reading all stats for dataset ***************
            }

            // Assuming the rows and columns are the same for all the datasets
            for (int row = 0; row < dataset_rows; ++row) {
                for (int col = 0; col < dataset_cols; ++col) {

                    // Change from row & columns to pixel number
                    int pixel =  row * dataset_cols + col;
                    // std::cout << row << "," << col << "\n";

                    int16_t landsat_avg = landsat_dataset_avg[pixel];
                    int16_t landsat_min = landsat_dataset_min[pixel];
                    int16_t landsat_max = landsat_dataset_max[pixel];
                    int16_t landsat_std = landsat_dataset_std[pixel];
                    int16_t landsat_cnt = landsat_dataset_cnt[pixel];
                    double landsat_per = landsat_dataset_per[pixel];

                    used_slack_doy = match_day[0]; // Use the first
                    int16_t modis_avg = modis_dataset_avg[pixel];
                    int16_t modis_min = modis_dataset_min[pixel];
                    int16_t modis_max = modis_dataset_max[pixel];
                    int16_t modis_std = modis_dataset_std[pixel];
                    int16_t modis_cnt = modis_dataset_cnt[pixel];
                    double modis_per = modis_dataset_per[pixel];
                    int16_t modis_sza = modis_dataset_sza[pixel];

                    // Filter Landsat by range and no_data
                    if ((landsat_avg == no_data) || (landsat_avg < vmin) || (landsat_avg > vmax)) {
                        // Ignore this pixel if no_data
                        continue;
                    }

                    // Filter MODIS by range and no_data
                    if ((modis_avg == no_data) || (modis_avg < vmin) || (modis_avg > vmax)) {
                        // Fallback to the second DOY MODIS dataset
                        if (match_day.size() > 1) {
                            used_slack_doy = match_day[1];
                            // std::cout << "MODIS out of range or NA, fallback to DOY=" << used_slack_doy << "\n";
                            
                            modis_avg = modis_dataset_avg_2[pixel];
                            modis_min = modis_dataset_min_2[pixel];
                            modis_max = modis_dataset_max_2[pixel];
                            modis_std = modis_dataset_std_2[pixel];
                            modis_cnt = modis_dataset_cnt_2[pixel];
                            modis_per = modis_dataset_per_2[pixel];
                            modis_sza = modis_dataset_sza_2[pixel];
                        }

                    }

                    // Filter MODIS again by range and no_data
                    if ((modis_avg == no_data) || (modis_avg < vmin) || (modis_avg > vmax)) {
                        // Fallback to the third DOY MODIS dataset
                        if (match_day.size() > 2) {
                            used_slack_doy = match_day[2];
                            // std::cout << "MODIS out of range or NA, fallback to DOY=" << used_slack_doy << "\n";

                            modis_avg = modis_dataset_avg_3[pixel];
                            modis_min = modis_dataset_min_3[pixel];
                            modis_max = modis_dataset_max_3[pixel];
                            modis_std = modis_dataset_std_3[pixel];
                            modis_cnt = modis_dataset_cnt_3[pixel];
                            modis_per = modis_dataset_per_3[pixel];
                            modis_sza = modis_dataset_sza_3[pixel];
                        }
                    }

                    // Filter MODIS for the last time by range and no_data
                    if ((modis_avg == no_data) || (modis_avg < vmin) || (modis_avg > vmax)) {
                        // MODIS still invalid. Ignore this pixel if no_data
                        // std::cout << "MODIS still out of range or NA, skipping this pixel.\n";
                        continue;
                    }
                    
                    // Concatenate stats and write line to CSV
                    // use dataset counter to decide what file to write into
                    if (ds == 0) {
                        // Write Landsat
                        csvFile << date << "," << year << "," << doy << ","
                                << pixel << "," << row << "," << col
                                << "," << landsat_min << "," << landsat_max << ","
                                << landsat_avg << "," << landsat_std << "," <<
                                landsat_cnt << "," << landsat_per;
                        // Write MODIS
                        if (list_input_dirs.size() > 1) {
                            csvFile << "," << doy +  used_slack_doy << "," << modis_min << ","
                                    << modis_max << "," << modis_avg << "," << modis_std << ","
                                    << modis_cnt << "," << modis_per << "," << modis_sza;
                        }
                        csvFile << "\n";
                    }
                    if (ds == 1) {
                        // Write Landsat
                        csvFile2 << date << "," << year << "," << doy << ","
                                 << pixel << "," << row << "," << col
                                 << "," << landsat_min << "," << landsat_max << ","
                                 << landsat_avg << "," << landsat_std << ","
                                 << landsat_cnt << "," << landsat_per;
                        // Write MODIS
                        if (list_input_dirs.size() > 1) {
                            csvFile2 << "," << doy + used_slack_doy << "," << modis_min << ","
                                     << modis_max << "," << modis_avg << "," << modis_std << ","
                                     << modis_cnt << "," << modis_per << "," << modis_sza;
                        }
                        csvFile2 << "\n";
                    }
                    if (ds == 2) {
                        // Write Landsat
                        csvFile3 << date << "," << year << "," << doy << ","
                                 << pixel << "," << row << "," << col
                                 << "," << landsat_min << "," << landsat_max << ","
                                 << landsat_avg << "," << landsat_std << ","
                                 << landsat_cnt << "," << landsat_per;
                        // Write MODIS
                        if (list_input_dirs.size() > 1) {
                            csvFile3 << "," << doy +  used_slack_doy << "," << modis_min << ","
                                     << modis_max << "," << modis_avg << "," << modis_std << ","
                                     << modis_cnt << "," << modis_per << "," << modis_sza;
                        }
                        csvFile3 << "\n";
                    }

                }
            }
            // Done processing all pixels

            // Dataset counter
            ds++;
        }
        std::cout << "*************************************************************************************\n";
    }

    // Close the CSV files
    csvFile.close();
    std::cout << "***CSV file: " << stats_file1 << " created successfully.\n";

    if (datasets.size() > 1) {
        csvFile2.close();
        std::cout << "***CSV file: " << stats_file2 << " created successfully.\n";
    }

    if (datasets.size() > 2) {
        csvFile3.close();
        std::cout << "***CSV file: " << stats_file3 << " created successfully.\n";
    }

    return true;
}

std::string match_file_doy(const int &year, const int &doy,
                           const std::vector<std::string> &list_files
                          ) {
    std::string year_str;
    std::string doy_str;
    std::string date_pattern;
    std::string matched_file;

    // Look for the file that contains the year and doy in file name
    year_str = std::to_string(year);
    doy_str  = std::to_string(doy);
    while (doy_str.length() < 3) {  // Pad day of year to 3 digits
        doy_str = "0" + doy_str;
    }

    // Combine into final string to make date filter: AYYYYDOY
    date_pattern = ".A" + year_str + doy_str;

    // std::cout << "--Matching for date pattern: " << date_pattern << "\n";

    for (size_t j = 0; j < list_files.size(); ++j) {
        // Find the file with the date filter
        if (find_pattern(list_files[j], date_pattern)) {
            // Date found in file name

            matched_file = list_files[j];
            break;

        }
    }

    return matched_file;
}

int find_int_position(const std::vector<int>& vec, int value) {
    auto it = std::find(vec.begin(), vec.end(), value);
    if (it != vec.end())
        return static_cast<int>(std::distance(vec.begin(), it));
    return -1;
}

void extract_year(const std::string& filename, int& year) {
    // Regex: look for ".A" followed by 4 digits (YYYY)
    std::regex pattern("\\.A(\\d{4})");
    std::smatch match;

    if (!std::regex_search(filename, match, pattern)) {
        // If no match found, return empty string or throw
        year = 0;
        return;  
    }

    year = std::stoi(match[1].str());
}
