/************************************************************************************************
 * Extraction of time series from HDF files in multiple dirs                                    *
 * Eduardo Jimenez Hernandez <eduardojh@arizona.edu> - October 2025                             *
 * Changelog:                                                                                   *
 * - 2025-10-26: Initial code                                                                   *
 * - 2025-10-27: First version two  simultaneous (Landsat & MODIS) time series                  *
 * - 2025-10-30: Fixed DOY not present in EVI & EVI2                                            *
 * - 2025-11-14: Fallback to retrieve values from +1/-1 days when indalid MODIS data, pixelwise *
 * - 2025-11-19: Fixed formula to calculate pixel position                                      * 
 ************************************************************************************************/

#include <cstdint>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

#include "gdal_priv.h"

#include "utils.hpp"

const size_t MAX_DATASETS = 3;
const size_t MAX_BOUNDS = 4;

int main(int argc, char* argv[]) {

    std::string format = "hdf4";  // default to hdf4
    std::string input_dirs;
    std::string output_dir;
    std::string datasets;
    std::string bounds;
    std::string pattern;
    std::string name;
    std::string lbound;
    std::string ubound;
    // bool debug = false;
    bool slack = false;

    std::vector<std::string> list_input_dirs;
    std::vector<std::string> list_bounds;
    std::vector<std::string> list_datasets;
    std::vector<std::string> list_patterns;
    
    double lower_bound = -10000;
    double upper_bound = 10000;
    int16_t no_data = -13000;

    // IMPORTANT: when addding a new option, don't forget to update the
    // short option, followed by colon if the argument is required.
    const char* const short_opts = "hvf:i:o:d:p:n:b:sg";

    const option long_opts[] = {
         {"help",       no_argument,       nullptr, 'h'},
         {"version",    no_argument,       nullptr, 'v'},
         {"format",     required_argument, nullptr, 'f'},
         {"idir",       required_argument, nullptr, 'i'},
         {"odir",       required_argument, nullptr, 'o'},
         {"datasets",   required_argument, nullptr, 'd'},
         {"pattern",    required_argument, nullptr, 'p'},
         {"name",       required_argument, nullptr, 'n'},
         {"bounds",     required_argument, nullptr, 'b'},
         {"slack",      no_argument,       nullptr, 's'},
        {"lower",      no_argument,       nullptr, 'l'},
        {"upper",      no_argument,       nullptr, 'u'},
        // {"debug",      no_argument,       nullptr, 'g'},
        {nullptr,      0,                 nullptr,  0 }
    };

    while (true) {
        const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

        if (opt == -1)
            break;

        switch (opt) {
            case 'h':
                print_help();
                return EXIT_SUCCESS;
            case 'v':
                print_version(APP_VERSION, APP_DATETIME);
                return EXIT_SUCCESS;
            case 'f':
                format = optarg;
                break;
            case 'i':
                input_dirs = optarg;
                break;
            case 'o':
                output_dir = optarg;
                break;
            case 'd':
                datasets = optarg;
                break;
            case 'p':
                pattern = optarg;
                break;
            case 'n':
                name = optarg;
                break;
            case 'b':
                bounds = optarg;
                break;
            case 's':
                slack = true;
                break;
            case 'l':
                lbound = optarg;
                break;
            case 'u':
                ubound = optarg;
                break;
            // case 'g':
            //     debug = true;
            //     break;
            case '?': // Unrecognized option
            default:
                print_help();
                return EXIT_FAILURE;
        }
    }

    std::cout << "\nExtraction of time series from HDF files in multiple dirs.\n";

    // Check required arguments

    if (format != "hdf4" && format != "hdf5") {
        std::cerr << "Invalid format: " << format << "\n";
        print_help();
        return EXIT_FAILURE;
    }

    if (input_dirs.empty()) {
        std::cerr << "ERROR: Input directory paths is required.\n\n";
        print_help();
        return EXIT_FAILURE;
    }

    if (output_dir.empty()) {
        std::cerr << "ERROR: Output directory path is required.\n\n";
        print_help();
        return EXIT_FAILURE;
    }

    if (datasets.empty()) {
        std::cout << "ERROR: No datasets to read.\n";
        print_help();
        return EXIT_FAILURE;
    }

    // if (bounds.empty()) {
    //     std::cerr << "ERROR: Bounds for (initial and final) rows and columns are required.\n\n";
    //     print_help();
    //     return EXIT_FAILURE;
    // }

    if (name.empty()) {
        std::cout << "ERROR: Outpit file name base is required.\n";
        print_help();
        return EXIT_FAILURE;
    }

    if (!lbound.empty()) {
        lower_bound = std::stod(lbound);
    } else {
        lower_bound = -10000;
    }

    if (!ubound.empty()) {
        upper_bound = std::stod(ubound);
    } else {
        upper_bound = 10000;
    }

    // Check directories
    list_input_dirs = split_by_commas(input_dirs);
    std::cout << "Processing " << input_dirs.size() << " input directories.\n";
    std::cout << "(assuming 1=LANDSAT, 2=MODIS, 3=VIIRS. 2 & 3 are optional)\n";

    for (std::string idir : list_input_dirs) {
        std::cout << "Input directory: " << idir << "... okay!\n";
        if (!directory_exists(idir)) {
            std::cerr << "ERROR: Input directory not found: " << idir << std::endl;
            return EXIT_FAILURE;
        }
    }

    std::cout << "Output directory: " << output_dir << "\n";
    if (!directory_exists(output_dir)) {
        std::cerr << "ERROR: Output directory not found: " << output_dir << std::endl;
        return EXIT_FAILURE;
    }

    // Check the number of datasets
    list_datasets = split_by_commas(datasets); // Dataset names from user
    std::cout << "Datasets: ";
    for (std::string ds : list_datasets) {
        std::cout << ds << " ";
    }
    std::cout << "(" << list_datasets.size() << " total datasets)." << std::endl;
    if (list_datasets.size() > MAX_DATASETS) {
        std::cerr << "ERROR: Sorry only " << MAX_DATASETS << " datasets allowed." << std::endl;
        return EXIT_FAILURE;
    }
    if (list_datasets.empty()) {
        std::cout << "No datasets to extract." << std::endl;
        return EXIT_SUCCESS;
    }

    // Check the bounds: initial/final rows and columns to extract time series from
    if (!bounds.empty()) {
        list_bounds = split_by_commas(bounds);
        if (list_bounds.size() != MAX_BOUNDS) {
            std::cerr << "ERROR: Mismatch in bounds: " << list_bounds.size()
                    << " provided, but" << MAX_BOUNDS << " expected." << std::endl;
            return EXIT_FAILURE;
        }
        std::cout << "Bounds: ";
        std::cout << " -Initial row:" << list_bounds[0] << "\n"
                  << " -Initial col: " << list_bounds[1] << "\n"
                  << " -Final row: " << list_bounds[2] << "\n"
                  << " -Final col: " << list_bounds[3] << "\n"
                  << std::endl;
    } else {
        std::cout << "Bounds not provided using all rows and columns." << std::endl;
    }

    // Pattern is optional
    if (!pattern.empty()) {
        list_patterns = split_by_commas(pattern);
        std::cout << "Patterns: ";
        for (std::string pat : list_patterns) {
            std::cout << pat << " ";
            std::cout << "(" << list_patterns.size() << " total patterns)." << std::endl;
        }
    }

    // Show initial arguments
    std::cout << "File format: " << format << "\n";
    std::cout << "Output file name base: " << name << "\n";

    // Get all HDF files in Landsat directory
    std::cout << "Scanning dir: " << list_input_dirs[0] << ". ";
    std::vector<std::string> file_list_1 = list_files(list_input_dirs[0], ".hdf");
    // Sort in alphabetical order
    std::sort(file_list_1.begin(), file_list_1.end());
    std::cout << "Found " << file_list_1.size() << " files.\n";

    // Get all HDF files in MODIS directory
    std::vector<std::string> file_list_2;
    if (list_input_dirs.size() > 1) {
        std::cout << "Scanning dir: " << list_input_dirs[1] << ". ";
        file_list_2 = list_files(list_input_dirs[1], ".hdf");
        std::sort(file_list_2.begin(), file_list_2.end());
        std::cout << "Found " << file_list_2.size() << " files.\n";
    }

    if (!list_patterns.empty()) {
        // Process each pattern
        for (std::string single_pattern : list_patterns) {
            // Filter the files by each pattern and process separately
            
            std::vector<std::string> filtered_list_1;
            std::vector<std::string> filtered_list_2;
            
            // First list of files
            std::cout << "Filtering files by: " << single_pattern << "...";
            for (size_t i=0; i < file_list_1.size(); ++i) {
                // If pattern provided and found in file name, keep it
                if (!single_pattern.empty() && find_pattern(file_list_1[i], single_pattern)) {
                    filtered_list_1.push_back(file_list_1[i]);
                }
            }
            std::cout << " found " << filtered_list_1.size() << " files.\n";

            if (filtered_list_1.empty()) {
                std::cout << "Nothing to do.\n";
                continue;
            }

            // Second list of files
            std::cout << "Filtering files by: " << single_pattern << "...";
            for (size_t i=0; i < file_list_2.size(); ++i) {
                // If pattern provided and found in file name, keep it
                if (!single_pattern.empty() && find_pattern(file_list_2[i], single_pattern)) {
                    filtered_list_2.push_back(file_list_2[i]);
                }
            }
            std::cout << " found " << filtered_list_2.size() << " files.\n";

            if (filtered_list_2.empty()) {
                std::cout << "Nothing to do.\n";
                continue;
            }

            // Create a file name per year, try using year
            std::string name_pattern = name + "_" + single_pattern;

            int year = 0;
            extract_year(single_pattern, year);
            if (year != 0) {
                name_pattern = name + "_" + std::to_string(year);
            }

            // Create a time series with the filtered files
            files2ts(filtered_list_1,
                     filtered_list_2,
                     list_input_dirs,
                     list_datasets,
                     output_dir,
                     name_pattern,
                     lower_bound,
                     upper_bound,
                     slack,
                     no_data);
        }
    } else {
        // Process all the files in dir
        files2ts(file_list_1,
                 file_list_2,
                 list_input_dirs,
                 list_datasets,
                 output_dir,
                 name,
                 lower_bound,
                 upper_bound,
                 slack,
                 no_data);
    }
    std::cout << "Ice never dies!" << std::endl;
    return EXIT_SUCCESS;
}
