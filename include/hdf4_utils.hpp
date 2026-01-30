#pragma once

#include <string>
#include <vector>

#include "mfhdf.h"  // HDF4 Scientific Dataset API

std::string get_hdf4_type_name(int32_t data_type);
std::vector<std::string> read_hdf4_info(const std::string& file_path);
bool read_hdf4_dataset_by_index(const std::string &, const int &, std::vector<int16_t> &, int &, int &);

template <typename T>
bool read_hdf4_by_index(const std::string &filename, const int &index, std::vector<T> &output_data, int &rows, int &cols);