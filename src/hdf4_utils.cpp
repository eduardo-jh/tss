#include <iostream>
#include <iomanip>
#include "hdf4_utils.hpp"

std::string get_hdf4_type_name(int32_t data_type) {
    // Warning: doublecheck
    switch (data_type) {
        case DFNT_CHAR:        return "DFNT_CHAR";
        case DFNT_UCHAR:       return "DFNT_UCHAR";
        case DFNT_INT8:        return "DFNT_INT8";
        case DFNT_UINT8:       return "DFNT_UINT8";
        case DFNT_INT16:       return "DFNT_INT16";
        case DFNT_UINT16:      return "DFNT_UINT16";
        case DFNT_INT32:       return "DFNT_INT32";
        case DFNT_UINT32:      return "DFNT_UINT32";
        case DFNT_FLOAT32:     return "DFNT_FLOAT32";
        case DFNT_FLOAT64:     return "DFNT_FLOAT64";
        case DFNT_NINT8:       return "DFNT_NINT8";
        case DFNT_NUINT8:      return "DFNT_NUINT8";
        case DFNT_NINT16:      return "DFNT_NINT16";
        case DFNT_NUINT16:     return "DFNT_NUINT16";
        case DFNT_NINT32:      return "DFNT_NINT32";
        case DFNT_NUINT32:     return "DFNT_NUINT32";
        case DFNT_NFLOAT32:    return "DFNT_NFLOAT32";
        case DFNT_NFLOAT64:    return "DFNT_NFLOAT64";
        default:               return "UNKNOWN_TYPE";
    }
}

std::vector<std::string> read_hdf4_info(const std::string& file_path) {
    std::vector<std::string> dataset_names;

    int32_t sd_id = SDstart(file_path.c_str(), DFACC_READ);
    if (sd_id < 0) {
        std::cerr << "ERROR: Failed to open HDF4 file: " << file_path << std::endl;
        return dataset_names;
    }

    int32_t num_datasets = 0, num_file_attrs = 0;
    if (SDfileinfo(sd_id, &num_datasets, &num_file_attrs) != 0) {
        std::cerr << "ERROR: Failed to get HDF4 file info: " << file_path << std::endl;
        SDend(sd_id);
        return dataset_names;
    }

    std::cout << "  HDF4 file contains: "
              << num_datasets << " datasets, "
              << num_file_attrs << " file attributes.\n";

    for (int32_t i = 0; i < num_datasets; ++i) {
        int32_t sds_id = SDselect(sd_id, i);
        if (sds_id == -1) continue;

        char name[H4_MAX_NC_NAME];
        int32_t rank, data_type, n_attrs;
        int32_t dimsizes[H4_MAX_VAR_DIMS];

        if (SDgetinfo(sds_id, name, &rank, dimsizes, &data_type, &n_attrs) == 0) {
            dataset_names.emplace_back(name);

            std::cout << "  " << i << " → " << name << " [";
            for (int j = 0; j < rank; ++j) {
                std::cout << dimsizes[j];
                if (j < rank - 1) std::cout << " x ";
            }
            std::cout << "] TYPE=" << get_hdf4_type_name(data_type) << "(" << data_type << ")" << std::endl;
        }

        SDendaccess(sds_id);
    }

    SDend(sd_id);
    return dataset_names;
}

bool read_hdf4_dataset_by_index(const std::string &filename, const int &index, std::vector<int16_t> &output_data, int &rows, int &cols) {
    int32_t sd_id = SDstart(filename.c_str(), DFACC_READ);
    if (sd_id == -1) {
        std::cerr << "ERROR: Failed to open HDF file: " << filename << std::endl;
        return false;
    }

    int32_t n_datasets = 0, n_attrs = 0;
    if (SDfileinfo(sd_id, &n_datasets, &n_attrs) != 0) {
        std::cerr << "ERROR: Failed to get HDF file info." << std::endl;
        SDend(sd_id);
        return false;
    }

    if (index < 0 || index >= n_datasets) {
        std::cerr << "ERROR: Index " << index << " out of range. Total datasets: " << n_datasets << std::endl;
        SDend(sd_id);
        return false;
    }

    // Select dataset by index
    int32_t sds_id = SDselect(sd_id, index);
    if (sds_id == -1) {
        std::cerr << "ERROR: Failed to select dataset at index " << index << std::endl;
        SDend(sd_id);
        return false;
    }

    // Get dataset info
    char name[H4_MAX_NC_NAME];
    int32_t rank, dimsizes[H4_MAX_VAR_DIMS], data_type, n_attrs_ds;
    if (SDgetinfo(sds_id, name, &rank, dimsizes, &data_type, &n_attrs_ds) != 0) {
        std::cerr << "ERROR: Failed to get dataset info." << std::endl;
        SDendaccess(sds_id);
        SDend(sd_id);
        return false;
    }

    std::cout << "==>Reading dataset: [" << name << "]\n";

    // Only rank 2 is supported
    if (rank != 2) {
        std::cerr << "ERROR: Rank is " << rank << ", only 2D are supported. Exiting." << std::endl;
        SDendaccess(sds_id);
        SDend(sd_id);
        return false;
    }

    std::cout << "  Rank: " << rank << "\n  Dimensions: ";
    // int rows = dimsizes[0];
    // int cols = dimsizes[1];
    rows = dimsizes[0];
    cols = dimsizes[1];
    std::cout << "(" << rows << "x" << cols << ")" << "\n";

    // Values for display
    int rmax = std::min(rows, 5);
    int cmax = std::min(cols, 5);

    // Read all the rows and columns
    int32_t start[2] = {0, 0};
    int32_t edges[2] = {rows, cols};

    // Allocate buffer depending on data type
    switch (data_type) {
        case DFNT_INT8: {
            std::cout << "  Data type: INT8\n";
            std::vector<int8_t> buffer(rows * cols);
            if (SDreaddata(sds_id, start, nullptr, edges, buffer.data()) == 0) {
                std::cout << "  Preview: [" << rmax << "x" << cmax << "]\n";
                for (int r = 0; r < rmax; ++r) {
                    for (int c = 0; c < cmax; ++c) {
                        std::cout << std::setw(6) << static_cast<int>(buffer[r * cols + c]);
                    }
                    std::cout << std::endl;
                }
            } else {
                std::cerr << "\n  ERROR: Failed to read INT8 data." << std::endl;
            }
            break;
        }

        case DFNT_INT16: {
            std::cout << "  Data type: INT16\n";
            std::vector<int16_t> buffer(rows * cols);
            if (SDreaddata(sds_id, start, nullptr, edges, buffer.data()) == 0) {
                std::cout << "  Preview: [" << rmax << "x" << cmax << "]\n";
                for (int r = 0; r < rmax; ++r) {
                    for (int c = 0; c < cmax; ++c) {
                        std::cout << std::setw(10) << std::fixed << std::setprecision(4)
                                  << buffer[r * cols + c];
                    }
                    std::cout << std::endl;
                }
                output_data = buffer;
                std::cout << "Returning data...ok.\n";
            } else {
                std::cerr << "\n  ERROR: Failed to read INT16 data.\n";
            }
            break;
        }

        case DFNT_FLOAT32: {
            std::cout << "  Data type: FLOAT32\n";
            std::vector<float> buffer(rows * cols);
            if (SDreaddata(sds_id, start, nullptr, edges, buffer.data()) == 0) {
                std::cout << "  Preview: [" << rmax << "x" << cmax << "]\n";
                for (int r = 0; r < rmax; ++r) {
                    for (int c = 0; c < cmax; ++c) {
                        std::cout << std::setw(8) << buffer[r * cols + c];
                    }
                    std::cout << std::endl;
                }
            } else {
                std::cerr << "\n  ERROR: Failed to read INT16 data." << std::endl;
            }
            break;
        }

        default: {
            std::cerr << "  ERROR: Dataset type " << get_hdf4_type_name(data_type) << " (" << data_type << ") unsupported." << std::endl;
            // SDendaccess(sds_id);
            // SDend(sd_id);
            // return false;
        }

    }

    SDendaccess(sds_id);
    SDend(sd_id);
    return true;
}

template <typename T>
bool read_hdf4_by_index(const std::string &filename, const int &index, std::vector<T> &output_data, int &rows, int &cols) {
    int32_t sd_id = SDstart(filename.c_str(), DFACC_READ);
    if (sd_id == -1) {
        std::cerr << "ERROR: Failed to open HDF file: " << filename << std::endl;
        return false;
    }

    int32_t n_datasets = 0, n_attrs = 0;
    if (SDfileinfo(sd_id, &n_datasets, &n_attrs) != 0) {
        std::cerr << "ERROR: Failed to get HDF file info." << std::endl;
        SDend(sd_id);
        return false;
    }

    if (index < 0 || index >= n_datasets) {
        std::cerr << "ERROR: Index " << index << " out of range. Total datasets: " << n_datasets << std::endl;
        SDend(sd_id);
        return false;
    }

    // Select dataset by index
    int32_t sds_id = SDselect(sd_id, index);
    if (sds_id == -1) {
        std::cerr << "ERROR: Failed to select dataset at index " << index << std::endl;
        SDend(sd_id);
        return false;
    }

    // Get dataset info
    char name[H4_MAX_NC_NAME];
    int32_t rank, dimsizes[H4_MAX_VAR_DIMS], data_type, n_attrs_ds;
    if (SDgetinfo(sds_id, name, &rank, dimsizes, &data_type, &n_attrs_ds) != 0) {
        std::cerr << "ERROR: Failed to get dataset info." << std::endl;
        SDendaccess(sds_id);
        SDend(sd_id);
        return false;
    }

    // std::cout << "==>Reading dataset: [" << name << "]\n";

    // Only rank 2 is supported
    if (rank != 2) {
        std::cerr << "ERROR: Rank is " << rank << ", only 2D are supported. Exiting." << std::endl;
        SDendaccess(sds_id);
        SDend(sd_id);
        return false;
    }

    std::cout << "  Rank: " << rank << "\n  Dimensions: ";
    // int rows = dimsizes[0];
    // int cols = dimsizes[1];
    rows = dimsizes[0];
    cols = dimsizes[1];
    std::cout << "(" << rows << "x" << cols << ")" << "\n";

    // Values for display
    int rmax = std::min(rows, 5);
    int cmax = std::min(cols, 5);

    // Read all the rows and columns
    int32_t start[2] = {0, 0};
    int32_t edges[2] = {rows, cols};

    // Allocate buffer depending on data type, then read
    std::vector<T> buffer(rows * cols);
    if (SDreaddata(sds_id, start, nullptr, edges, buffer.data()) == 0) {

        // Show a preview depending on type
        switch (data_type) {
            case DFNT_INT8: {
                std::cout << "  Data type: INT8\n";
                std::cout << "  Preview: [" << rmax << "x" << cmax << "]\n";
                for (int r = 0; r < rmax; ++r) {
                    for (int c = 0; c < cmax; ++c) {
                        std::cout << std::setw(6) << static_cast<int>(buffer[r * cols + c]);
                    }
                    std::cout << std::endl;
                }
                break;
            }

            case DFNT_INT16: {
                std::cout << "  Data type: INT16\n";
                std::cout << "  Preview: [" << rmax << "x" << cmax << "]\n";
                for (int r = 0; r < rmax; ++r) {
                    for (int c = 0; c < cmax; ++c) {
                        std::cout << std::setw(10) << std::fixed << std::setprecision(4)
                                  << buffer[r * cols + c];
                    }
                    std::cout << std::endl;
                }
                break;
            }

            case DFNT_FLOAT32: {
                std::cout << "  Data type: FLOAT32\n";
                std::cout << "  Preview: [" << rmax << "x" << cmax << "]\n";
                for (int r = 0; r < rmax; ++r) {
                    for (int c = 0; c < cmax; ++c) {
                        std::cout << std::setw(8) << buffer[r * cols + c];
                    }
                    std::cout << std::endl;
                }
                break;
            }

            default: {
                std::cout << "  Data type: " << data_type << ", no preview available.\n";
            }
        }

        // Return the data
        output_data = buffer;
        std::cout << "Returning data...ok.\n";

    } else {
        std::cerr << "\n  ERROR: Failed to read " << data_type << " dataset." << std::endl;
    }
    

    // switch (data_type) {
    //     case DFNT_INT8: {
    //         std::cout << "  Data type: INT8\n";
    //         // std::vector<int8_t> buffer(rows * cols);
    //         std::vector<T> buffer(rows * cols);
    //         if (SDreaddata(sds_id, start, nullptr, edges, buffer.data()) == 0) {
    //             std::cout << "  Preview: [" << rmax << "x" << cmax << "]\n";
    //             for (int r = 0; r < rmax; ++r) {
    //                 for (int c = 0; c < cmax; ++c) {
    //                     std::cout << std::setw(6) << static_cast<int>(buffer[r * cols + c]);
    //                 }
    //                 std::cout << std::endl;
    //             }
    //             output_data = buffer;
    //             std::cout << "Returning data...ok.\n";
    //         } else {
    //             std::cerr << "\n  ERROR: Failed to read INT8 data." << std::endl;
    //         }
    //         break;
    //     }

    //     case DFNT_INT16: {
    //         std::cout << "  Data type: INT16\n";
    //         // std::vector<int16_t> buffer(rows * cols);
    //         std::vector<T> buffer(rows * cols);
    //         if (SDreaddata(sds_id, start, nullptr, edges, buffer.data()) == 0) {
    //             std::cout << "  Preview: [" << rmax << "x" << cmax << "]\n";
    //             for (int r = 0; r < rmax; ++r) {
    //                 for (int c = 0; c < cmax; ++c) {
    //                     std::cout << std::setw(10) << std::fixed << std::setprecision(4)
    //                               << buffer[r * cols + c];
    //                 }
    //                 std::cout << std::endl;
    //             }
    //             output_data = buffer;
    //             std::cout << "Returning data...ok.\n";
    //         } else {
    //             std::cerr << "\n  ERROR: Failed to read INT16 data.\n";
    //         }
    //         break;
    //     }

    //     case DFNT_FLOAT32: {
    //         std::cout << "  Data type: FLOAT32\n";
    //         // std::vector<float> buffer(rows * cols);
    //         std::vector<T> buffer(rows * cols);
    //         if (SDreaddata(sds_id, start, nullptr, edges, buffer.data()) == 0) {
    //             std::cout << "  Preview: [" << rmax << "x" << cmax << "]\n";
    //             for (int r = 0; r < rmax; ++r) {
    //                 for (int c = 0; c < cmax; ++c) {
    //                     std::cout << std::setw(8) << buffer[r * cols + c];
    //                 }
    //                 std::cout << std::endl;
    //             }
    //             output_data = buffer;
    //             std::cout << "Returning data...ok.\n";
    //         } else {
    //             std::cerr << "\n  ERROR: Failed to read FLOAT32 data." << std::endl;
    //         }
    //         break;
    //     }

    //     default: {
    //         std::cerr << "  ERROR: Dataset type " << get_hdf4_type_name(data_type) << " (" << data_type << ") unsupported." << std::endl;
    //         // SDendaccess(sds_id);
    //         // SDend(sd_id);
    //         // return false;
    //     }

    // }

    SDendaccess(sds_id);
    SDend(sd_id);
    return true;
}

// Explicit template instantiations (so linker finds them)
template bool read_hdf4_by_index(const std::string &, const int &, std::vector<int8_t> &, int &, int &);
template bool read_hdf4_by_index(const std::string &, const int &, std::vector<uint8_t> &, int &, int &);
template bool read_hdf4_by_index(const std::string &, const int &, std::vector<int16_t> &, int &, int &);
template bool read_hdf4_by_index(const std::string &, const int &, std::vector<uint16_t> &, int &, int &);
template bool read_hdf4_by_index(const std::string &, const int &, std::vector<int32_t> &, int &, int &);
template bool read_hdf4_by_index(const std::string &, const int &, std::vector<uint32_t> &, int &, int &);
template bool read_hdf4_by_index(const std::string &, const int &, std::vector<float> &, int &, int &);
template bool read_hdf4_by_index(const std::string &, const int &, std::vector<double> &, int &, int &);
