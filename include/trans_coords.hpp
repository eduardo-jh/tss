#pragma once

#include <string>
#include <vector>
#include <gdal_priv.h>
#include <cpl_conv.h>

struct GDALTranslateOpts {
    int in_rows;
    int in_cols;
    int out_rows;
    int out_cols;
    std::string resampling = "average"; // default
    std::string a_srs = ""; // empty means no reprojection
    float no_data = -13000.0f;
};

// struct GDALWarpOpts {
//     int inRows;
//     int inCols;
//     int outRows;
//     int outCols;
//     std::string s_srs;
//     std::string t_srs;
//     std::vector<double> extent;       // [minX, minY, maxX, maxY]
//     std::string resampling = "average";
// };
struct GDALWarpOpts {
    int in_rows;   // input rows
    int in_cols;   // input cols
    int out_rows;  // desired output rows
    int out_cols;  // desired output cols
    std::string s_srs;        // source CRS
    std::string t_srs;        // target CRS
    std::vector<double> extent;
    std::string resampling = "average"; // resampling method
    double nodata = -13000.0; // optional no-data values
};


struct UTMInfo {
    bool isUTM;
    int zone;
    std::string hemisphere; // "N" or "S"
};

void utm_to_latlon(double, double, int, char);
void latlon_to_utm(double, double, int, char);
std::vector<double> utm2geo(double, double, int, char);
std::vector<double> geo2utm(double, double);
void process_coordinates(const std::string &);
void parse_coordinates(const std::string &);
std::vector<double> utm2geo_extent(const std::string &, const std::string &);
double angle_to_distance(const double &, const double &);
double degrees2meters(const double &);
double distance_to_angle(const double &, const double &);
double meters2degrees(const double &);
bool saveVectorAsGeoTIFFLatLon(const std::string &, const std::vector<int> &, const size_t &, const size_t &, const double[6]);
std::vector<double> adjustExtentMODIS(const std::vector<double> &, const std::vector<double> &, const size_t &, const size_t &, const bool & = false);
std::vector<int> matchExtentRowCol(std::vector<double> &, const std::vector<double> &, const size_t &, const size_t &, const bool & = false);
// std::vector<double> extent_from_tile(const std::string &);
// bool reproject_tif(const std::string &, const std::string &, const std::string &, const std::string &);
bool reproject_tif(const std::string &, const std::string &, const std::string &, const std::string &, const std::string & = "nearest");
bool reproject_to_memory(const std::string &, std::vector<float> &, double[6], const GDALWarpOpts &);
bool reproject_tif_ds(const std::string &, const std::string &, std::vector<int16_t> &, int &, int &, const std::string &, const std::string &, const std::string & = "nearest");
bool reproject_tif2ds(const std::string &, const std::string &, std::vector<int16_t> &, double[6], GDALWarpOpts &);

// ---- Type mapping (C++ type -> GDALDataType) ----
template <typename T> struct GDALTypeMap;

template <> struct GDALTypeMap<int8_t>   { static constexpr GDALDataType type = GDT_Byte;   };
template <> struct GDALTypeMap<uint8_t>  { static constexpr GDALDataType type = GDT_Byte;   };
template <> struct GDALTypeMap<int16_t>  { static constexpr GDALDataType type = GDT_Int16;  };
template <> struct GDALTypeMap<uint16_t> { static constexpr GDALDataType type = GDT_UInt16; };
template <> struct GDALTypeMap<int32_t>  { static constexpr GDALDataType type = GDT_Int32;  };
template <> struct GDALTypeMap<uint32_t> { static constexpr GDALDataType type = GDT_UInt32; };
template <> struct GDALTypeMap<float>    { static constexpr GDALDataType type = GDT_Float32;};
template <> struct GDALTypeMap<double>   { static constexpr GDALDataType type = GDT_Float64;};

// ---- Generic resample function ----
template <typename T>
bool resample_dataset(const std::vector<T> &, std::vector<T> &, const GDALTranslateOpts &);

template <typename T>
bool saveToTiff(const std::string &, const std::vector<T> &, int, int, const double[6], const std::string &, double);

template <typename T>
bool reproject_dataset(const std::vector<T> &inputData,
                    std::vector<T> &outputData,
                    double geoTransform[6],
                    const GDALWarpOpts &opts);
UTMInfo interpretSRS(const std::string &srsStr);