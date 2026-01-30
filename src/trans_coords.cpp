#include <cstdint>
#include <string>
#include <cmath>
#include <cstddef>
#include <ios>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <fstream>

#include <gdal_utils.h>
#include <gdalwarper.h>
#include <ogr_spatialref.h>

#include "trans_coords.hpp"
#include "utils.hpp"

void utm_to_latlon(double x, double y, int zone, char hemisphere) {
    OGRSpatialReference utmSRS;
    utmSRS.SetUTM(zone, (hemisphere == 'N' || hemisphere == 'n'));
    // utmSRS.SetUTM(zone, hemisphere);
    utmSRS.SetWellKnownGeogCS("WGS84");

    OGRSpatialReference latlonSRS;
    latlonSRS.SetWellKnownGeogCS("WGS84");

    OGRCoordinateTransformation* coordTrans =
        OGRCreateCoordinateTransformation(&utmSRS, &latlonSRS);

    if (!coordTrans) {
        std::cerr << "Failed to create coordinate transformation.\n";
        return;
    }

    if (!coordTrans->Transform(1, &y, &x)) {
        std::cerr << "Coordinate transformation failed.\n";
        OCTDestroyCoordinateTransformation(coordTrans);
        return;
    }

    std::cout << "Latitude: " << y << ", Longitude: " << x << std::endl;
    OCTDestroyCoordinateTransformation(coordTrans);
}

void latlon_to_utm(double lat, double lon, int zone, char hemisphere) {
    // std::cout << "Received lon=" << lon << ", lat=" << lat << "\n";

    OGRSpatialReference latlonSRS;
    latlonSRS.SetWellKnownGeogCS("WGS84");

    OGRSpatialReference utmSRS;
    utmSRS.SetUTM(zone, (hemisphere == 'N' || hemisphere == 'n'));
    // utmSRS.SetUTM(zone, hemisphere);
    utmSRS.SetWellKnownGeogCS("WGS84");

    OGRCoordinateTransformation* coordTrans =
        OGRCreateCoordinateTransformation(&latlonSRS, &utmSRS);

    if (!coordTrans) {
        std::cerr << "Failed to create coordinate transformation.\n";
        return;
    }

    // Note: GDAL expects lat, lon as input
    // but returns Easting and Northing, WTF?
    // https://gdal.org/en/stable/tutorials/osr_api_tut.html
    if (!coordTrans->Transform(1, &lat, &lon)) {
        std::cerr << "Coordinate transformation failed.\n";
        OCTDestroyCoordinateTransformation(coordTrans);
        return;
    }

    char hemisphere_char = hemisphere ? 'N' : 'S';
    std::cout << std::fixed << std::setprecision(2)
            //   << "Easting: " << lon << ", Northing: " << lat;
              << "Northing: " << lon << ", Easting: " << lat;
    std::cout << ", Zone: " << zone << hemisphere_char << std::endl;
    OCTDestroyCoordinateTransformation(coordTrans);
}

std::vector<double> utm2geo(double x, double y, int zone, char hemisphere) {

    std::vector<double> geo_coords = {0, 0};

    OGRSpatialReference utmSRS;
    utmSRS.SetUTM(zone, (hemisphere == 'N' || hemisphere == 'n'));
    // utmSRS.SetUTM(zone, hemisphere);
    utmSRS.SetWellKnownGeogCS("WGS84");

    OGRSpatialReference latlonSRS;
    latlonSRS.SetWellKnownGeogCS("WGS84");

    OGRCoordinateTransformation* coordTrans =
        OGRCreateCoordinateTransformation(&utmSRS, &latlonSRS);

    if (!coordTrans) {
        std::cerr << "Failed to create coordinate transformation.\n";
        return geo_coords;
    }

    if (!coordTrans->Transform(1, &y, &x)) {
        std::cerr << "Coordinate transformation failed.\n";
        OCTDestroyCoordinateTransformation(coordTrans);
        return geo_coords;
    }

    OCTDestroyCoordinateTransformation(coordTrans);

    // std::cout << "Latitude: " << y << ", Longitude: " << x << "\n";

    geo_coords[0] = x;
    geo_coords[1] = y;

    return geo_coords;
}

std::vector<double> geo2utm(double lat, double lon) {
    
    std::vector<double> utm_coords = {0, 0, 0, 0};
    int zone = static_cast<int>((lon / 6.0) + 31);  // works for +/- lon
    bool hemisphere = lat >= 0; // N=true, S=false

    OGRSpatialReference latlonSRS;
    latlonSRS.SetWellKnownGeogCS("WGS84");

    OGRSpatialReference utmSRS;
    // utmSRS.SetUTM(zone, hemisphere == 'N' || hemisphere == 'n'));
    utmSRS.SetUTM(zone, hemisphere);
    utmSRS.SetWellKnownGeogCS("WGS84");

    OGRCoordinateTransformation* coordTrans =
        OGRCreateCoordinateTransformation(&latlonSRS, &utmSRS);

    if (!coordTrans) {
        std::cerr << "ERROR: Failed to create coordinate transformation.\n";
        return utm_coords;
    }

    // Note: GDAL expects lat, lon as input
    // but returns Easting and Northing, WTF?
    // https://gdal.org/en/stable/tutorials/osr_api_tut.html
    if (!coordTrans->Transform(1, &lat, &lon)) {
        std::cerr << "Coordinate transformation failed.\n";
        OCTDestroyCoordinateTransformation(coordTrans);
        return utm_coords;
    }

    OCTDestroyCoordinateTransformation(coordTrans);

    // std::cout << std::fixed << std::setprecision(2)
    //           << "Northing: " << lon << ", Easting: " << lat;
    // std::cout << ", Zone: " << zone <<  hemisphere ? 'N' : 'S'; << std::endl;

    utm_coords[0] = lon; // northing returned by GDAL = y
    utm_coords[1] = lat; // easting returned by GDAL = x
    utm_coords[2] = zone;
    utm_coords[3] = hemisphere; // N=true, S=false

    return utm_coords;
}

void process_coordinates(const std::string &data) {
    std::istringstream iss(data);  // input stream string
    std::string record;

    // Split by whitespace (space or newline)
    while (iss >> record) {
        std::vector<std::string> parts;
        std::stringstream ss(record);
        std::string item;

        // Split by comma
        while (std::getline(ss, item, ',')) {
            parts.push_back(item);
        }

        if (parts.size() == 2) {
            // Lat/Lon
            double lon = std::stod(parts[0]);
            double lat = std::stod(parts[1]);
            // int zone = static_cast<int>(std::floor((lon + 180.0) / 6.0) + 1);
            int zone = static_cast<int>((lon / 6.0) + 31);  // works for +/- lon
            // int hemisphere = (lat >= 0) ? 1 : 0;  // 1 for N, 0 for S
            // char hemisphere_char = hemisphere ? 'N' : 'S';
            char hemisphere = (lat >= 0) ? 'N' : 'S';
            std::cout << "Lon=" << lon << ", Lat=" << lat << ", (zone="
                    //   << zone << hemisphere_char << " " << hemisphere << ") -> ";
                    << zone << hemisphere << ") -> ";
            latlon_to_utm(lat, lon, zone, hemisphere);
        }
        else if (parts.size() == 4) {
            // UTM
            double x = std::stod(parts[0]);
            double y = std::stod(parts[1]);
            int zone = std::stoi(parts[2]);
            // int hemisphere = std::stoi(parts[3]); // 1 or 0
            char hemisphere = parts[3].empty() ? 'N' : parts[3][0];
            utm_to_latlon(x, y, zone, hemisphere);
        }
        else {
            std::cerr << "Invalid coordinate format: " << record << std::endl;
        }
    }
}

void parse_coordinates(const std::string &coordinates) {
    // Check if the coordinates is file
    bool is_file;
    std::ifstream f(coordinates);
    is_file = f.good();

    if (is_file) {
        // If path is an existing file open it
        std::cout << "Openning file: " << coordinates << "\n";

        std::ifstream file(coordinates);
        if (!file.is_open()) {
            std::cerr << "ERROR: file not found or error opening file.\n";
            return;
        }
        
        std::string line;
        while (std::getline(file, line)) {
            if (!line.empty()) {
                // Process line by line
                process_coordinates(line);
            }
        }
    } else {
        // Process the inline coordinates
        process_coordinates(coordinates);
    }
}

std::vector<double> utm2geo_extent(const std::string &coord_list, const std::string &utmzone) {
    // Expected: x1,y1,x2,y2; for Upper-Left=(x1,y1) and Lower-Right=(x2,y2) corners
    // UTM Zone as: 12N
    std::vector<double> geo_coords = {0, 0, 0, 0};
    char northing = 'X';
    int zone = 0;

    // Get the UTM Zone
    switch (utmzone.length()) {
        case 2: {
            // std::cout << "UTM 2 chars, ";
            zone = std::stoi(utmzone.substr(0, 1));
            northing = utmzone.substr(1, 1)[0];
            break;
        }
        case 3: {
            // std::cout << "UTM 3 chars, ";
            zone = std::stoi(utmzone.substr(0, 2));
            northing = utmzone.substr(2, 1)[0];
            break;
        }
        default: {
            std::cerr << "\nERROR: Invalid UTM Zone of " << utmzone.length() << " chars.\n";
            return geo_coords;
        }
    }

    if (northing != 'N' && northing != 'n' && northing != 'S' && northing != 's') {
        std::cerr << "ERROR: Invalid northing, expected N or S, got " << northing << " instead.\n";
        return geo_coords;
    }
    if (zone <= 0 || zone > 60) {
        std::cerr << "ERROR: Invalid UTM Zone " << zone << "\n";
        return geo_coords;
    }
    // std::cout << "UTM Zone: " << zone << " " << northing << "\n";

    // Split coordinates
    std::vector<std::string> in_coords = split_by_commas(coord_list);
    // std::cout << " (Coordinates length: " << in_coords.size() << ")\n";
    double x1 = std::stod(in_coords[0]);
    double y1 = std::stod(in_coords[1]);
    double x2 = std::stod(in_coords[2]);
    double y2 = std::stod(in_coords[3]);
    // std::cout << std::fixed << std::setprecision(2) << "Extent UTM: ulx=" << x1
    //           << ", uly=" << y1 << ", lrx=" << x2 << ", lry=" << y2 
    //           <<  " (x-diff=" << abs(x2 - x1) << ", y-diff=" << abs(y2 - y1) << ")" << std::endl;

    std::vector<double> ul;
    std::vector<double> lr;
    // ul = utm2geo(x1, y1, zone, northing);
    // lr = utm2geo(x2, y2, zone, northing);
    // GDAL is crazy, Lat/Lon expected
    ul = utm2geo(y1, x1, zone, northing);
    lr = utm2geo(y2, x2, zone, northing);
    
    // std::cout << "Extent upper-left: lon=" << ul[0] << ", lat=" << ul[1] << "\n";
    // std::cout << "Extent lower-right: lon=" << lr[0] << ", lat=" << lr[1] << "\n";

    // geo_coords[0] = static_cast<float>(ul[0]);
    // geo_coords[1] = static_cast<float>(ul[1]);
    // geo_coords[2] = static_cast<float>(lr[0]);
    // geo_coords[3] = static_cast<float>(lr[1]);
    geo_coords[0] = ul[0];
    geo_coords[1] = ul[1];
    geo_coords[2] = lr[0];
    geo_coords[3] = lr[1];

    // std::cout << "Extent upper-left: lon=" << geo_coords[0] << ", lat=" << geo_coords[1] << "\n";
    // std::cout << "Extent lower-right: lon=" << geo_coords[2] << ", lat=" << geo_coords[3] << "\n";

    return geo_coords;
}

double angle_to_distance(const double &radius_m, const double &angle_d) {
    // Returns distance in meters, inputs radius (m) and angle (degrees)
    // std::cout << "PI=" << M_PI << "\n";
    return (radius_m * M_PI * angle_d) / 180.0;
}

double degrees2meters(const double &angle_d) {
    // Returns distance in meters using equatorial radius of Earth, angle (degrees)
    // double earth_radius = 12756000.0/2.0; // 12,756 kilometers
    double earth_radius = 6378137.0; // in meters
    return angle_to_distance(earth_radius, angle_d);
}

double distance_to_angle(const double &radius_m, const double &distance) {
    // Returns angle (degrees), inputs radius (m) and distance (m)
    // std::cout << "PI=" << M_PI << "\n";
    return (distance * 180.0) / (radius_m * M_PI);
}

double meters2degrees(const double &distance) {
    // Returns angle (degrees), inputs radius (m) and distance (m)
    double earth_radius = 6378137.0; // in meters
    return distance_to_angle(earth_radius, distance);
}

bool saveVectorAsGeoTIFFLatLon(const std::string &filename,
                         const std::vector<int> &data,
                         const size_t &width,
                         const size_t &height,
                         const double geoTransform[6]) {

    std::cout << "Saving GeoTIFF file, Lat/Lon projection.\n";
    std::cout << "  nrows=" << height << "\n";
    std::cout << "  ncols=" << width << "\n";
    std::cout << "  size=" << height*width << "\n";
    std::cout << "  gt=[" << geoTransform[0] << "," << geoTransform[1] << "," << geoTransform[2]
              << "," <<  geoTransform[3] << "," <<  geoTransform[4] << "," <<  geoTransform[5] << "]\n";

    if (data.size() != width * height) {
        std::cerr << "ERROR: Data vector size must be " << width * height << " elements.\n";
        return false;
    }

    // GDALAllRegister();

    const char *pszFormat = "GTiff";
    GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
    if (!poDriver) {
        std::cerr << "ERROR: GTiff driver not available.\n";
        return false;
    }

    char **papszOptions = NULL;
    GDALDataset *poDstDS = poDriver->Create(filename.c_str(), width, height, 1, GDT_Int32, papszOptions);
    if (!poDstDS) {
        std::cerr << "ERROR: Could not create output file.\n";
        return false;
    }
    poDstDS->SetGeoTransform(const_cast<double*>(geoTransform));

    // Set projection to WGS84
    OGRSpatialReference srs;
    srs.SetWellKnownGeogCS("WGS84");
    char *projectionWKT = NULL;
    srs.exportToWkt(&projectionWKT);
    poDstDS->SetProjection(projectionWKT);
    CPLFree(projectionWKT);

    GDALRasterBand *poBand = poDstDS->GetRasterBand(1);
    if (!poBand) {
        std::cerr << "ERROR: Could not get raster band.\n";
        GDALClose(poDstDS);
        return false;
    }

    // Write the entire dataset in one call
    CPLErr err = poBand->RasterIO(GF_Write, 0, 0, width, height,
                                  (void*)data.data(), width, height, GDT_Int32,
                                  0, 0);
    if (err != CE_None) {
        std::cerr << "ERROR: RasterIO write failed.\n";
        GDALClose(poDstDS);
        return false;
    }

    // Set NoData value (optional)
    poBand->SetNoDataValue(-9999);

    // Flush and close
    GDALClose(poDstDS);
    std::cout << "GeoTIFF saved: " << filename << std::endl;
    return true;
}

// std::vector<double> extent_from_tile(const std::string &modis_tile) {
//     // Returns the extent in Lat/Lon coordinates for the MODIS tile
//     // Adjust to existing MODIS grid
//     // h09v05: 40,-110,30,-100
//     // h11v08: 10,-70,0,-60
//     // h11v11: -20,-74,-30,-64
//     // h17v07: 20,-10,10,0

//     std::vector<double> extent = {0,0,0,0}; // {ulx, uly, lrx, lry}

//     if (modis_tile == "h09v05") {
//         extent[0] = -110;
//         extent[1] = 40;
//         extent[2] = -100;
//         extent[3] = 30;
//     } else if (modis_tile == "h11v08") {
//         extent[0] = -70;
//         extent[1] = 10;
//         extent[2] = -60;
//         extent[3] = 0;
//     } else if (modis_tile == "h11v11") {
//         extent[0] = -74;
//         extent[1] = -20;
//         extent[2] = -64;
//         extent[3] = -30;
//     } else if (modis_tile == "h17v07") {
//         extent[0] = -10;
//         extent[1] = 20;
//         extent[2] = 0;
//         extent[3] = 10;
//     } else {
//         std::cout << "ERROR: MODIS scene not recognized." << std::endl;
//         return extent;
//     }

//     return extent;
// }

std::vector<double> adjustExtentMODIS(const std::vector<double> &extent, const std::vector<double> &base_extent, const size_t &m_rows, const size_t &m_cols, const bool &at) {
    // at = all touch criteria, includes pixels touched by the edges,
    //      if false only pixels inside will be considered

    std::vector<double> adj_extent = {0,0,0,0};

    double ulx = base_extent[0];
    double uly = base_extent[1];
    double lrx = base_extent[2];
    double lry = base_extent[3];

    double x_res = abs(ulx - lrx) / m_cols;
    double y_res = abs(uly - lry) / m_rows;

    // Find pixel size resolution
    std::cout << "Adjusting extent (at=" << std::boolalpha << at << ")\n";
    std::cout << "  x_res=" << x_res << ", y_res=" << y_res << "\n";

    for (size_t col=0; col < m_cols-1; ++col) {
        // Increase coordinates by pixel size resolution
        double x_ini = ulx + (x_res * col);
        double x_end = ulx + (x_res * (col + 1));

        // Adjusted x-axis coordinate for upper-left corner
        if ((x_ini <= extent[0]) && (x_end >= extent[0])) {
            if (at) {
                adj_extent[0] = x_ini; // all-touch
                // std::cout << "  ulx_adj(at)=" << adj_extent[0] << " [" << x_ini << " < " << extent[0] << " < " << x_end << "]\n";
            } else {
                adj_extent[0] = x_end;  // inside pixels only
                // std::cout << "  ulx_adj=" << adj_extent[0] << " [" << x_ini << " < " << extent[0] << " < " << x_end << "]\n";
            }
        }

        // Adjusted x-axis coordinate for lower-right corner
        if ((x_ini <= extent[2]) && (x_end >= extent[2])) {
            if (at) {
                adj_extent[2] = x_end; // all-touch
                // std::cout << "  lrx_adj(at)=" << adj_extent[2] << " [" << x_ini << " < " << extent[2] << " < " << x_end << "]\n";
            } else {
                adj_extent[2] = x_ini;  // inside pixels only
                // std::cout << "  lrx_adj=" << adj_extent[2] << " [" << x_ini << " < " << extent[2] << " < " << x_end << "]\n";
            }
        }

    }

    for (size_t row=0; row < m_rows-1; ++row) {
        // Increase coordinates by pixel size resolution
        double y_ini = lry + (y_res * row);
        double y_end = lry + (y_res * (row + 1));

        // Adjusted y-axis coordinate for lower-right corner
        if ((y_ini <= extent[3]) && (y_end >= extent[3])) {
            if (at) {
                adj_extent[3] = y_ini; // all-touch
                // std::cout << "  lry_adj(at)=" << adj_extent[3] << " [" << y_ini << " < " << extent[3] << " < " << y_end << "]\n";
            } else {
                adj_extent[3] = y_end;  // inside pixels only
                // std::cout << "  lry_adj=" << adj_extent[3] << " [" << y_ini << " < " << extent[3] << " < " << y_end << "]\n";
            }
        }

        // Adjusted y-axis coordinate for upper-left corner
        if ((y_ini <= extent[1]) && (y_end >= extent[1])) {
            if (at) {
                adj_extent[1] = y_end; // all-touch
                // std::cout << "  uly_adj(at)=" << adj_extent[1] << " [" << y_ini << " < " << extent[1] << " < " << y_end << "]\n";
            } else {
                adj_extent[1] = y_ini;  // inside pixels only
                // std::cout << "  uly_adj=" << adj_extent[1] << " [" << y_ini << " < " << extent[1] << " < " << y_end << "]" << std::endl;
            }
        }
    }

    return adj_extent;
}

std::vector<int> matchExtentRowCol(std::vector<double> &extent, const std::vector<double> &base_extent, const size_t &m_rows, const size_t &m_cols, const bool &at) {
    // Return the row and columns that match the extent
    std::vector<int> indices = {0,0,0,0};
    std::vector<double> adj_extent = {0,0,0,0};

    double ulx = base_extent[0];
    double uly = base_extent[1];
    double lrx = base_extent[2];
    double lry = base_extent[3];

    // std::cout << "matching cols=" << std::fixed << std::setprecision(6) << (ulx - lrx) / static_cast<double>(m_cols) << "\n";
    // std::cout << "matching rows=" << std::fixed << std::setprecision(6) << (uly - lry) / static_cast<double>(m_rows) << "\n";

    // double x_res = abs(ulx - lrx) / m_cols;
    // double y_res = abs(uly - lry) / m_rows;
    double x_res = (ulx - lrx) / static_cast<double>(m_cols);
    double y_res = (uly - lry) / static_cast<double>(m_rows);

    // Find pixel size resolution
    std::cout << "Match extent indices (at=" << std::boolalpha << at << ")\n";
    std::cout << std::fixed << std::setprecision(6) << "  x_res=" << x_res << ", y_res=" << y_res << "\n";

    for (size_t col=0; col < m_cols-1; ++col) {
        // Increase coordinates by pixel size resolution
        double x_ini = ulx - (x_res * col);
        double x_end = ulx - (x_res * (col + 1));

        // Adjusted x-axis coordinate for upper-left corner
        if ((x_ini <= extent[0]) && (x_end >= extent[0])) {
            size_t col_index = col;
            if (at) {
                adj_extent[0] = x_ini; // all-touch
                indices[0] = col_index;
                std::cout << std::fixed << std::setprecision(6) 
                          << "  ulx_adj(at)=" << adj_extent[0]
                          << " [" << x_ini << " < " << extent[0]
                          << " < " << x_end << "] (Found col=" << col_index << ", Final col=" << indices[0] << ")\n";
            } else {
                adj_extent[0] = x_end;  // inside pixels only
                indices[0] = col_index+1;
                std::cout << std::fixed << std::setprecision(6)
                << "  ulx_adj=" << adj_extent[0]
                << " [" << x_ini << " < " << extent[0]
                << " < " << x_end << "] (Found col=" << col_index << ", Final col=" << indices[0] << ")\n";
            }
        }

        // Adjusted x-axis coordinate for lower-right corner
        if ((x_ini <= extent[2]) && (x_end >= extent[2])) {
            size_t col_index = col;
            if (at) {
                adj_extent[2] = x_end; // all-touch
                indices[2] = col_index+1;
                std::cout << std::fixed << std::setprecision(6)
                          << "  lrx_adj(at)=" << adj_extent[2]
                          << " [" << x_ini << " < " << extent[2]
                          << " < " << x_end << "] (Found col=" << col_index << ", Final col=" << indices[2] << ")\n";
            } else {
                adj_extent[2] = x_ini;  // inside pixels only
                indices[2] = col_index;
                std::cout << std::fixed << std::setprecision(6)
                          << "  lrx_adj=" << adj_extent[2]
                          << " [" << x_ini << " < " << extent[2]
                          << " < " << x_end << "] (Found col=" << col_index << ", Final col=" << indices[2] << ")\n";
            }
        }
    }

    for (size_t row=0; row < m_rows-1; ++row) {
        // Increase coordinates by pixel size resolution
        double y_ini = uly - (y_res * row);
        double y_end = uly - (y_res * (row + 1));

        // Adjusted y-axis coordinate for lower-right corner
        if ((y_ini >= extent[3]) && (y_end <= extent[3])) {
            size_t row_index = row;
            // std::cout << "  row=" << row_index << "\n";
            if (at) {
                adj_extent[3] = y_ini; // all-touch
                indices[3] = row_index+1;
                std::cout << std::fixed << std::setprecision(6) 
                          << "  lry_adj(at)=" << adj_extent[3]
                          << " [" << y_ini << " > " << extent[3]
                          << " > " << y_end << "] (Found row=" << row_index << ", Final row=" << indices[3] << ")\n";
            } else {
                adj_extent[3] = y_end;  // inside pixels only
                indices[3] = row_index;
                std::cout << std::fixed << std::setprecision(6)
                          << "  lry_adj=" << adj_extent[3]
                          << " [" << y_ini << " > " << extent[3]
                          << " > " << y_end << "] (Found row=" << row_index << ", Final row=" << indices[3] << ")\n";
            }
        }

        // Adjusted y-axis coordinate for upper-left corner
        if ((y_ini >= extent[1]) && (y_end <= extent[1])) {
            size_t row_index = row;
            // std::cout << "  row=" << row_index << "\n";
            if (at) {
                adj_extent[1] = y_end; // all-touch
                indices[1] = row_index;
                std::cout << std::fixed << std::setprecision(6)
                          << "  uly_adj(at)=" << adj_extent[1]
                          << " [" << y_ini << " < " << extent[1]
                          << " > " << y_end << "] (Found row=" << row_index << ", Final row=" << indices[1] << ")\n";
            } else {
                adj_extent[1] = y_ini;  // inside pixels only
                indices[1] = row_index+1;
                std::cout << std::fixed << std::setprecision(6)
                          << "  uly_adj=" << adj_extent[1]
                          << " [" << y_ini << " > " << extent[1]
                          << " > " << y_end << "] (Found row=" << row_index << ", Final row=" << indices[1] << ")" << std::endl;
            }
        }
    }

    // Update the extent
    extent = adj_extent;

    return indices;
}

template <typename T>
bool resample_dataset(const std::vector<T> &input,
                   std::vector<T> &output,
                   const GDALTranslateOpts &opts)
{
    if (input.size() != static_cast<size_t>(opts.in_rows * opts.in_cols)) {
        std::cerr << "ERROR: Input size does not match rows*cols" << std::endl;
        return false;
    }

    std::cout << "**Running resampling with GDALTranslate\n";

    GDALAllRegister();
    GDALDriver *memDriver = GetGDALDriverManager()->GetDriverByName("MEM");
    if (!memDriver) {
        std::cerr << "ERROR: MEM driver not available" << std::endl;
        return false;
    }

    // Create input dataset
    GDALDataset *srcDS = memDriver->Create("", opts.in_cols, opts.in_rows, 1,
                                           GDALTypeMap<T>::type, nullptr);
    if (!srcDS) return false;

    GDALRasterBand *srcBand = srcDS->GetRasterBand(1);
    srcBand->SetNoDataValue(opts.no_data);

    if (srcBand->RasterIO(GF_Write, 0, 0, opts.in_cols, opts.in_rows,
                          const_cast<T*>(input.data()), opts.in_cols, opts.in_rows,
                          GDALTypeMap<T>::type, 0, 0) != CE_None)
    {
        std::cerr << "ERROR: Failed writing input data" << std::endl;
        GDALClose(srcDS);
        return false;
    }

    int outCols = static_cast<int>(opts.out_cols);
    int outRows = static_cast<int>(opts.out_rows);

    // Build args for gdal_translate
    std::string of = "GTiff";
    char **papszArgv = nullptr;
    // papszArgv = CSLSetNameValue(papszArgv, "of", of.c_str());
    // std::cout << "  Setting GDALTranslate options: of=" << of.c_str() << "\n";
    papszArgv = CSLSetNameValue(papszArgv, "r", opts.resampling.c_str());
    std::cout << "  Setting GDALTranslate options: r=" << opts.resampling.c_str() << "\n";
    papszArgv = CSLAddString(papszArgv, "-outsize");
    std::cout << "  Setting GDALTranslate options: -outsize=";
    papszArgv = CSLAddString(papszArgv, std::to_string(outCols).c_str());
    std::cout << " " << std::to_string(outCols).c_str() << " ";
    papszArgv = CSLAddString(papszArgv, std::to_string(outRows).c_str());
    std::cout << std::to_string(outRows).c_str() << "\n";
    if (!opts.a_srs.empty()) {
        papszArgv = CSLSetNameValue(papszArgv, "a_srs", opts.a_srs.c_str());
        std::cout << "  Setting GDALTranslate options: a_srs=" << opts.a_srs.c_str() << "\n";
    }

    GDALTranslateOptions *translateOpts = GDALTranslateOptionsNew(papszArgv, nullptr);
    CSLDestroy(papszArgv);

    // ChatGPT suggestion, creates a file named "MEM:::"
    // GDALDataset *translated = (GDALDataset*) GDALTranslate("MEM:::", srcDS, translateOpts, nullptr);
    GDALDataset *translated = (GDALDataset*) GDALTranslate("/vsimem/resampled.tif", srcDS, translateOpts, nullptr);
    GDALTranslateOptionsFree(translateOpts);
    GDALClose(srcDS);

    if (!translated) {
        std::cerr << "ERROR: GDALTranslate failed" << std::endl;
        return false;
    }

    output.resize(outCols * outRows);
    GDALRasterBand *translatedBand = translated->GetRasterBand(1);

    if (translatedBand->RasterIO(GF_Read, 0, 0, outCols, outRows,
                                 output.data(), outCols, outRows,
                                 GDALTypeMap<T>::type, 0, 0) != CE_None)
    {
        std::cerr << "Failed reading translated data" << std::endl;
        GDALClose(translated);
        return false;
    }

    GDALClose(translated);
    return true;
}

template <typename T>
bool saveToTiff(const std::string &filename,
                const std::vector<T> &data,
                int rows, int cols,
                const double geoTransform[6],
                const std::string &a_srs,
                double no_data)
{
    // GDALAllRegister();
    GDALDriver *gtiffDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    if (!gtiffDriver) {
        std::cerr << "ERROR: GTiff driver not available" << std::endl;
        return false;
    }

    GDALDataset *dstDS = gtiffDriver->Create(filename.c_str(), cols, rows, 1,
                                             GDALTypeMap<T>::type, nullptr);
    if (!dstDS) {
        std::cerr << "ERROR: Failed to create output file: " << filename << std::endl;
        return false;
    }

    if (geoTransform) {
        if (dstDS->SetGeoTransform(const_cast<double*>(geoTransform)) != CE_None) {
            std::cerr << "Warning: Failed to set GeoTransform" << std::endl;
        }
    }

    if (!a_srs.empty()) {
        // dstDS->SetProjection(a_srs.c_str());
        OGRSpatialReference srs;
        if (srs.SetFromUserInput(a_srs.c_str()) == OGRERR_NONE) {
            char *wkt = nullptr;
            srs.exportToWkt(&wkt);
            dstDS->SetProjection(wkt);
            CPLFree(wkt);
        } else {
            std::cerr << "Warning: Invalid SRS string: " << a_srs << std::endl;
        }
    }

    GDALRasterBand *band = dstDS->GetRasterBand(1);
    band->SetNoDataValue(no_data);

    if (band->RasterIO(GF_Write, 0, 0, cols, rows,
                       const_cast<T*>(data.data()), cols, rows,
                       GDALTypeMap<T>::type, 0, 0) != CE_None)
    {
        std::cerr << "ERROR: Failed writing GeoTIFF" << std::endl;
        GDALClose(dstDS);
        return false;
    }

    GDALClose(dstDS);
    return true;
}

template <typename T>
bool reproject_dataset(const std::vector<T> &inputData,
                    std::vector<T> &outputData,
                    double geoTransform[6],
                    const GDALWarpOpts &opts)
{
    GDALAllRegister();

    std::cout << "**Running reprojection with GDALWarp\n";
    std::cout << "--nrows=" << opts.in_rows
              << "\n--ncols=" << opts.in_cols
              << "\n--out_rows=" << opts.out_rows
              << "\n--out_cols=" << opts.out_cols
              << "\n--s_srs=" << opts.s_srs
              << "\n--t_srs=" << opts.t_srs
              << "\n--resampling=" << opts.resampling
              << "\n--extent=[";
    // Extent should have 4 coordinates
    for (int i=0; i<4; ++i) {
        std::cout << opts.extent[i];
        if (i < 3) std::cout << " ";
    }
    std::cout << "]\n";
             

    // 1. Create the source in-memory dataset
    GDALDriver *memDriver = GetGDALDriverManager()->GetDriverByName("MEM");
    if (!memDriver) {
        std::cerr << "ERROR: MEM driver not available" << std::endl;
        return false;
    }

    GDALDataset *srcDS = memDriver->Create("", opts.in_cols, opts.in_rows, 1, GDALTypeMap<T>::type, nullptr);
    if (!srcDS) {
        std::cerr << "ERROR: Failed to create input MEM dataset" << std::endl;
        return false;
    }

    // // GeoTransform
    // double srcGT[6];
    // srcGT[0] = opts.extent[0];                           // top-left X (minX)
    // srcGT[1] = abs(opts.extent[2] - opts.extent[0]) / opts.inCols;    // pixel width
    // srcGT[2] = 0;                                   // rotation (0 for north-up)
    // srcGT[3] = opts.extent[3];                           // top-left Y (maxY)
    // srcGT[4] = 0;                                   // rotation
    // srcGT[5] = -abs(opts.extent[3] - opts.extent[1]) / opts.inRows;   // pixel height (negative -> north-up)

    // Set geotransform to source dataset
    std::cout << "  gt=[" << geoTransform[0]
              << " " << geoTransform[1]
              << " " << geoTransform[2]
              << " " << geoTransform[3]
              << " " << geoTransform[4]
              << " " << geoTransform[5] << "]\n";

    srcDS->SetGeoTransform(geoTransform);

    // Set projection to source dataset
    if (!opts.s_srs.empty()) {
        OGRSpatialReference srs;
        if (srs.SetFromUserInput(opts.s_srs.c_str()) == OGRERR_NONE) {
            char *wkt = nullptr;
            srs.exportToWkt(&wkt);
            srcDS->SetProjection(wkt);
            // std::cout << "  Setting source projection: " << wkt << "\n";
            CPLFree(wkt);
        } else {
            std::cerr << "ERROR: Invalid source SRS: " << opts.s_srs << std::endl;
            return false;
        }
    }

    // if (!opts.t_srs.empty()) {
    //     OGRSpatialReference srs;
    //     if (srs.SetFromUserInput(opts.t_srs.c_str()) == OGRERR_NONE) {
    //         char *wkt = nullptr;
    //         srs.exportToWkt(&wkt);
    //         std::cout << "  Target projection: " << wkt << "\n";
    //         CPLFree(wkt);
    //     } else {
    //         std::cerr << "ERROR: Invalid target SRS: " << opts.t_srs << std::endl;
    //         return false;
    //     }
    // }


    // Write input data into source dataset
    GDALRasterBand *srcBand = srcDS->GetRasterBand(1);
    if (srcBand->RasterIO(GF_Write, 0, 0, opts.in_cols, opts.in_rows,
                          const_cast<T*>(inputData.data()), opts.in_cols, opts.in_rows,
                          GDALTypeMap<T>::type, 0, 0) != CE_None) {
        std::cerr << "ERROR: Failed to write input raster" << std::endl;
        GDALClose(srcDS);
        return false;
    }

    // 2. Prepare Warp options
    char **warpArgs = nullptr;
    // warpArgs = CSLSetNameValue(warpArgs, "r", opts.resampling.c_str());
    // std::cout << "  Setting GDALWarp options: r=" << opts.resampling.c_str() << "\n";
    // if (opts.extent.size() == 4) {
    //     std::string teStr = std::to_string(opts.extent[0]) + " " +
    //                         std::to_string(opts.extent[1]) + " " +
    //                         std::to_string(opts.extent[2]) + " " +
    //                         std::to_string(opts.extent[3]);
    //     warpArgs = CSLSetNameValue(warpArgs, "te", teStr.c_str());
    //     std::cout << "  Setting GDALWarp options: te=" << teStr.c_str() << "\n";
    // }

    std::string tsStr = std::to_string(opts.out_cols) + " " + std::to_string(opts.out_rows);
    // std::string tsStr = std::to_string(opts.inCols) + " " + std::to_string(opts.inRows);
    warpArgs = CSLSetNameValue(warpArgs, "ts", tsStr.c_str());
    std::cout << "  Setting GDALWarp options: ts=" << tsStr.c_str() << "\n";

    // if (!opts.s_srs.empty()) {
    //     warpArgs = CSLSetNameValue(warpArgs, "s_srs", opts.s_srs.c_str());
    //     std::cout << "  Setting GDALWarp options: s_srs=" << opts.s_srs.c_str() << "\n";
    // }
    if (!opts.t_srs.empty()) {
        warpArgs = CSLSetNameValue(warpArgs, "t_srs", opts.t_srs.c_str());
        std::cout << "  Setting GDALWarp options: t_srs=" << opts.t_srs.c_str() << "\n";
    }
    
    // // -to SRC_METHOD=NO_GEOTRANSFORM
    // std::string toStr = "SRC_METHOD=NO_GEOTRANSFORM";
    // warpArgs = CSLSetNameValue(warpArgs, "to", toStr.c_str());
    // std::cout << "  Setting GDALWarp options: to=" << toStr.c_str() << "\n";
    // std::string toStr2 = "DST_METHOD=NO_GEOTRANSFORM";
    // warpArgs = CSLSetNameValue(warpArgs, "to", toStr2.c_str());
    // std::cout << "  Setting GDALWarp options: to=" << toStr2.c_str() << "\n";

    // // Trying something...
    // warpArgs = CSLSetNameValue(warpArgs, "src_method", "NO_GEOTRANSFORM");
    // std::cout << "  Setting GDALWarp options: src_method=NO_GEOTRANSFORM\n";

    GDALWarpAppOptions *warpOpts = GDALWarpAppOptionsNew(warpArgs, nullptr);

    // 3. Create output in-memory dataset
    // GDALDataset *dstDS = memDriver->Create("/vsimem/output_warp.tif", opts.outCols, opts.outRows, 1, GDALTypeMap<T>::type, nullptr);
    GDALDataset *dstDS = memDriver->Create("", opts.out_cols, opts.out_rows, 1, GDALTypeMap<T>::type, nullptr);
    if (!dstDS) {
        std::cerr << "ERROR: Failed to create output MEM dataset" << std::endl;
        GDALWarpAppOptionsFree(warpOpts);
        GDALClose(srcDS);
        return false;
    }

    // Perform warp
    // GDALDataset *warpDS = (GDALDataset*) GDALWarp("/vsimem/warped.tif", dstDS, 1, (GDALDatasetH*)&srcDS, warpOpts, nullptr);
    GDALDataset *warpDS = (GDALDataset*) GDALWarp("", dstDS, 1, (GDALDatasetH*)&srcDS, warpOpts, nullptr);
    if (!warpDS) {
        std::cerr << "ERROR: Warp failed" << std::endl;
        GDALWarpAppOptionsFree(warpOpts);
        GDALClose(srcDS);
        GDALClose(dstDS);
        return false;
    }

    // 4. Extract results back into output dataset
    outputData.resize(opts.out_rows * opts.out_cols);
    GDALRasterBand *outBand = warpDS->GetRasterBand(1);
    if (outBand->RasterIO(GF_Read, 0, 0, opts.out_cols, opts.out_rows,
                          outputData.data(), opts.out_cols, opts.out_rows,
                          GDALTypeMap<T>::type, 0, 0) != CE_None) {
        std::cerr << "ERROR: Failed to read output raster" << std::endl;
        GDALWarpAppOptionsFree(warpOpts);
        GDALClose(warpDS);
        return false;
    }

    // Retrieve geotransform from the warp data
    if (warpDS->GetGeoTransform(geoTransform) != CE_None) {
        std::cerr << "Warning: GeoTransform not available, filling default" << std::endl;
        for (int i = 0; i < 6; ++i) geoTransform[i] = 0.0;
    }

    // Cleanup
    GDALWarpAppOptionsFree(warpOpts);
    GDALClose(warpDS);
    GDALClose(srcDS);

    return true;
}

UTMInfo interpretSRS(const std::string &srsStr)
{
    UTMInfo info = {false, 0, ""};

    OGRSpatialReference srs;
    if (srs.SetFromUserInput(srsStr.c_str()) != OGRERR_NONE) {
        std::cerr << "ERROR: Failed to parse SRS string: " << srsStr << std::endl;
        return info;
    }

    int zone = 0;
    int north = 0;
    // Is it UTM?
    if (srs.GetUTMZone(&north) != 0) {
        zone = srs.GetUTMZone(&north);
        info.isUTM = true;
        info.zone = zone;
        info.hemisphere = (north != 0) ? "N" : "S";
    }

    return info;
}

// bool reproject_tif(const std::string &inFile,
//                    const std::string &outFile,
//                    const std::string &s_srs,
//                    const std::string &t_srs)
// {
//     // WORKS!
//     GDALAllRegister();

//     GDALDataset *srcDS = (GDALDataset *)GDALOpen(inFile.c_str(), GA_ReadOnly);
//     if (!srcDS)
//     {
//         std::cerr << "Failed to open source file: " << inFile << std::endl;
//         return false;
//     }

//     // Create warped VRT dataset (like gdalwarp in memory)
//     GDALDataset *warpedDS = (GDALDataset *)GDALAutoCreateWarpedVRT(
//         srcDS,
//         s_srs.c_str(),   // source SRS (can be nullptr if already in dataset)
//         t_srs.c_str(),   // target SRS
//         GRA_Bilinear,    // resampling (nearest, bilinear, cubic...)
//         0.0,             // error threshold
//         nullptr);        // warp options

//     if (!warpedDS)
//     {
//         std::cerr << "Warp VRT creation failed." << std::endl;
//         GDALClose(srcDS);
//         return false;
//     }

//     // Save warped dataset to GeoTIFF
//     GDALDriver *gtiffDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
//     if (!gtiffDriver)
//     {
//         std::cerr << "GTiff driver not available." << std::endl;
//         GDALClose(srcDS);
//         GDALClose(warpedDS);
//         return false;
//     }

//     GDALDataset *dstDS = gtiffDriver->CreateCopy(outFile.c_str(), warpedDS, FALSE, nullptr, nullptr, nullptr);
//     if (!dstDS)
//     {
//         std::cerr << "Failed to create output file: " << outFile << std::endl;
//         GDALClose(srcDS);
//         GDALClose(warpedDS);
//         return false;
//     }

//     // Cleanup
//     GDALClose(dstDS);
//     GDALClose(warpedDS);
//     GDALClose(srcDS);

//     return true;
// }

bool reproject_tif(const std::string &inFile,
                   const std::string &outFile,
                   const std::string &s_srs,
                   const std::string &t_srs,
                   const std::string &resampling) // default
{
    GDALAllRegister();

    // Open source file
    GDALDataset *srcDS = (GDALDataset *)GDALOpen(inFile.c_str(), GA_ReadOnly);
    if (!srcDS)
    {
        std::cerr << "Failed to open source file: " << inFile << std::endl;
        return false;
    }

    // Pick resampling algorithm
    GDALResampleAlg resampleAlg = GRA_NearestNeighbour;
    if (resampling == "bilinear") resampleAlg = GRA_Bilinear;
    else if (resampling == "cubic") resampleAlg = GRA_Cubic;
    else if (resampling == "cubicspline") resampleAlg = GRA_CubicSpline;
    else if (resampling == "lanczos") resampleAlg = GRA_Lanczos;
    else if (resampling == "average") resampleAlg = GRA_Average;
    else if (resampling == "mode") resampleAlg = GRA_Mode;
    else if (resampling == "max") resampleAlg = GRA_Max;
    else if (resampling == "min") resampleAlg = GRA_Min;
    else if (resampling == "med") resampleAlg = GRA_Med;
    else if (resampling == "q1") resampleAlg = GRA_Q1;
    else if (resampling == "q3") resampleAlg = GRA_Q3;

    // Create warped VRT (like gdalwarp in memory)
    GDALDataset *warpedDS = (GDALDataset *)GDALAutoCreateWarpedVRT(
        srcDS,
        s_srs.empty() ? nullptr : s_srs.c_str(), // can pass nullptr if dataset already has CRS
        t_srs.c_str(),
        resampleAlg,
        0.0,      // warp error threshold
        nullptr); // warp options

    if (!warpedDS)
    {
        std::cerr << "Warp VRT creation failed." << std::endl;
        GDALClose(srcDS);
        return false;
    }

    // Save warped dataset to GeoTIFF
    GDALDriver *gtiffDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    if (!gtiffDriver)
    {
        std::cerr << "GTiff driver not available." << std::endl;
        GDALClose(srcDS);
        GDALClose(warpedDS);
        return false;
    }

    GDALDataset *dstDS = gtiffDriver->CreateCopy(outFile.c_str(), warpedDS, FALSE, nullptr, nullptr, nullptr);
    if (!dstDS)
    {
        std::cerr << "Failed to create output file: " << outFile << std::endl;
        GDALClose(srcDS);
        GDALClose(warpedDS);
        return false;
    }

    // Cleanup
    GDALClose(dstDS);
    GDALClose(warpedDS);
    GDALClose(srcDS);

    return true;
}

bool reproject_tif_ds(const std::string &inFile,
                   const std::string &outFile,
                   std::vector<int16_t> &outDataset,
                   int &rows, int &cols,
                   const std::string &s_srs,
                   const std::string &t_srs,
                   const std::string &resampling) // default
{
    GDALAllRegister();

    // Open source file
    GDALDataset *srcDS = (GDALDataset *)GDALOpen(inFile.c_str(), GA_ReadOnly);
    if (!srcDS)
    {
        std::cerr << "Failed to open source file: " << inFile << std::endl;
        return false;
    }

    // Pick resampling algorithm
    GDALResampleAlg resampleAlg = GRA_NearestNeighbour;
    if (resampling == "bilinear") resampleAlg = GRA_Bilinear;
    else if (resampling == "cubic") resampleAlg = GRA_Cubic;
    else if (resampling == "cubicspline") resampleAlg = GRA_CubicSpline;
    else if (resampling == "lanczos") resampleAlg = GRA_Lanczos;
    else if (resampling == "average") resampleAlg = GRA_Average;
    else if (resampling == "mode") resampleAlg = GRA_Mode;
    else if (resampling == "max") resampleAlg = GRA_Max;
    else if (resampling == "min") resampleAlg = GRA_Min;
    else if (resampling == "med") resampleAlg = GRA_Med;
    else if (resampling == "q1") resampleAlg = GRA_Q1;
    else if (resampling == "q3") resampleAlg = GRA_Q3;

    // Create warped VRT (like gdalwarp in memory)
    GDALDataset *warpedDS = (GDALDataset *)GDALAutoCreateWarpedVRT(
        srcDS,
        s_srs.empty() ? nullptr : s_srs.c_str(), // can pass nullptr if dataset already has CRS
        t_srs.c_str(),
        resampleAlg,
        0.0,      // warp error threshold
        nullptr); // warp options

    if (!warpedDS)
    {
        std::cerr << "Warp VRT creation failed." << std::endl;
        GDALClose(srcDS);
        return false;
    }

    // Save warped dataset to GeoTIFF
    GDALDriver *gtiffDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    if (!gtiffDriver)
    {
        std::cerr << "GTiff driver not available." << std::endl;
        GDALClose(srcDS);
        GDALClose(warpedDS);
        return false;
    }

    GDALDataset *dstDS = gtiffDriver->CreateCopy(outFile.c_str(), warpedDS, FALSE, nullptr, nullptr, nullptr);
    if (!dstDS)
    {
        std::cerr << "Failed to create output file: " << outFile << std::endl;
        GDALClose(srcDS);
        GDALClose(warpedDS);
        return false;
    }

    // Allocate output buffer
    // int cols = warpedDS->GetRasterXSize();
    // int rows = warpedDS->GetRasterYSize();
    cols = warpedDS->GetRasterXSize();
    rows = warpedDS->GetRasterYSize();

    std::cout << "Warped dataset size: " 
              << cols << " x " << rows << std::endl;
    outDataset.resize(cols * rows);

    GDALRasterBand *band = warpedDS->GetRasterBand(1);
    if (band->RasterIO(GF_Read, 0, 0, cols, rows,
                       outDataset.data(), cols, rows,
                       GDT_Int16, 0, 0) != CE_None)
    {
        std::cerr << "ERROR: Failed to read warped raster into memory." << std::endl;
        GDALClose(srcDS);
        GDALClose(warpedDS);
        return false;
    }

    // Cleanup
    GDALClose(dstDS);
    GDALClose(warpedDS);
    GDALClose(srcDS);

    return true;
}

bool reproject_tif2ds(const std::string &inFile,
                   const std::string &outFile,
                   std::vector<int16_t> &outDataset,
                   double outGeoTransform[6],
                   GDALWarpOpts &opts)
{
    GDALAllRegister();

    // Open source file
    GDALDataset *srcDS = (GDALDataset *)GDALOpen(inFile.c_str(), GA_ReadOnly);
    if (!srcDS)
    {
        std::cerr << "Failed to open source file: " << inFile << std::endl;
        return false;
    }

    // Pick resampling algorithm
    GDALResampleAlg resampleAlg = GRA_NearestNeighbour;
    if (opts.resampling == "bilinear") resampleAlg = GRA_Bilinear;
    else if (opts.resampling == "cubic") resampleAlg = GRA_Cubic;
    else if (opts.resampling == "cubicspline") resampleAlg = GRA_CubicSpline;
    else if (opts.resampling == "lanczos") resampleAlg = GRA_Lanczos;
    else if (opts.resampling == "average") resampleAlg = GRA_Average;
    else if (opts.resampling == "mode") resampleAlg = GRA_Mode;
    else if (opts.resampling == "max") resampleAlg = GRA_Max;
    else if (opts.resampling == "min") resampleAlg = GRA_Min;
    else if (opts.resampling == "med") resampleAlg = GRA_Med;
    else if (opts.resampling == "q1") resampleAlg = GRA_Q1;
    else if (opts.resampling == "q3") resampleAlg = GRA_Q3;

    // This will automatically set the number of output rows/columns
    // // Create warped VRT (like gdalwarp in memory)
    // GDALDataset *warpedDS = (GDALDataset *)GDALAutoCreateWarpedVRT(
    //     srcDS,
    //     opts.s_srs.empty() ? nullptr : opts.s_srs.c_str(), // can pass nullptr if dataset already has CRS
    //     opts.t_srs.c_str(),
    //     resampleAlg,
    //     0.0,      // warp error threshold
    //     nullptr); // warp options

    // Create warp options
    GDALWarpOptions *psWarpOptions = GDALCreateWarpOptions();

    psWarpOptions->hSrcDS = srcDS;
    // psWarpOptions->hDstDS = dstDS;

    // Source projection: EPSG:32613
    OGRSpatialReference oSrcSRS;
    // oSrcSRS.importFromEPSG(32613);
    oSrcSRS.SetFromUserInput(opts.s_srs.c_str());
    char *pszSrcWKT = nullptr;
    oSrcSRS.exportToWkt(&pszSrcWKT);
    // psWarpOptions->pszSrcWKT = pszSrcWKT; // GDAL takes ownership

    // Target projection: EPSG:4326
    OGRSpatialReference oDstSRS;
    // oDstSRS.importFromEPSG(4326);
    oDstSRS.SetFromUserInput(opts.t_srs.c_str());
    char *pszDstWKT = nullptr;
    oDstSRS.exportToWkt(&pszDstWKT);
    // psWarpOptions->pszDstWKT = pszDstWKT; // GDAL takes ownership

    // Create transformer
    void *hTransformArg = GDALCreateGenImgProjTransformer(
                                        srcDS, pszSrcWKT,
                                        NULL, pszDstWKT,
                                        TRUE, 0.0, 1);

    psWarpOptions->pfnTransformer = GDALGenImgProjTransform;
    psWarpOptions->pTransformerArg = hTransformArg;

    // Set resampling algorithm
    psWarpOptions->eResampleAlg = resampleAlg;

    // Setup band mapping (input to output)
    int nBands = srcDS->GetRasterCount();
    psWarpOptions->nBandCount = nBands;
    psWarpOptions->panSrcBands = (int *) CPLMalloc(sizeof(int) * nBands);
    psWarpOptions->panDstBands = (int *) CPLMalloc(sizeof(int) * nBands);
    psWarpOptions->padfDstNoDataReal = (double*) CPLMalloc(sizeof(double)*nBands);
    for (int i = 0; i < nBands; i++) {
        psWarpOptions->panSrcBands[i] = i + 1;
        psWarpOptions->panDstBands[i] = i + 1;
        psWarpOptions->padfDstNoDataReal[i] = opts.nodata;
    }

    // psWarpOptions->padfDstNoDataReal = (double*) CPLMalloc(sizeof(double)*nBands);
    // for (int i = 0; i < nBands; i++) {
    //     psWarpOptions->padfDstNoDataReal[i] = opts.nodata;
    // }

    // Create the warped VRT
    GDALDataset *warpedDS = (GDALDataset *) GDALCreateWarpedVRT(
        srcDS,
        opts.out_cols,
        opts.out_rows,
        outGeoTransform,
        psWarpOptions
    );

    if (!warpedDS)
    {
        std::cerr << "Warp VRT creation failed." << std::endl;
        GDALClose(srcDS);
        return false;
    }
    std::cout << "Warped VRT created with size "
              << warpedDS->GetRasterXSize() << " x "
              << warpedDS->GetRasterYSize() << std::endl;

    // Save warped dataset to GeoTIFF
    GDALDriver *gtiffDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    if (!gtiffDriver)
    {
        std::cerr << "GTiff driver not available." << std::endl;
        GDALDestroyWarpOptions(psWarpOptions);
        GDALClose(srcDS);
        GDALClose(warpedDS);
        return false;
    }

    GDALDataset *dstDS = gtiffDriver->CreateCopy(outFile.c_str(), warpedDS, FALSE, nullptr, nullptr, nullptr);
    if (!dstDS)
    {
        std::cerr << "Failed to create output file: " << outFile << std::endl;
        GDALDestroyWarpOptions(psWarpOptions);
        GDALClose(srcDS);
        GDALClose(warpedDS);
        return false;
    }

    // Allocate output buffer
    // cols = warpedDS->GetRasterXSize();
    // rows = warpedDS->GetRasterYSize();

    std::cout << "Warped dataset size: " 
              << warpedDS->GetRasterYSize() << "x" << warpedDS->GetRasterXSize()
              << ", " << opts.out_rows << " x " << opts.out_cols << std::endl;
    outDataset.resize(opts.out_cols * opts.out_rows);

    GDALRasterBand *band = warpedDS->GetRasterBand(1);
    if (band->RasterIO(GF_Read, 0, 0, opts.out_cols, opts.out_rows,
                       outDataset.data(), opts.out_cols, opts.out_rows,
                       GDT_Int16, 0, 0) != CE_None)
    {
        std::cerr << "ERROR: Failed to read warped raster into memory." << std::endl;
        GDALDestroyWarpOptions(psWarpOptions);
        GDALClose(srcDS);
        GDALClose(warpedDS);
        return false;
    }

    // Cleanup
    GDALDestroyWarpOptions(psWarpOptions);
    GDALClose(dstDS);
    GDALClose(warpedDS);
    GDALClose(srcDS);

    return true;
}

bool reproject_to_memory(const std::string &inFile,
                         std::vector<float> &outData,
                         double geoTransform[6],
                         const GDALWarpOpts &opts)
{
    GDALAllRegister();

    // Open source dataset
    GDALDataset *srcDS = (GDALDataset *)GDALOpen(inFile.c_str(), GA_ReadOnly);
    if (!srcDS)
    {
        std::cerr << "ERROR: Failed to open source file: " << inFile << std::endl;
        return false;
    }

    // Pick resampling algorithm
    GDALResampleAlg resampleAlg = GRA_Average;
    if (opts.resampling == "nearest") resampleAlg = GRA_NearestNeighbour;
    else if (opts.resampling == "bilinear") resampleAlg = GRA_Bilinear;
    else if (opts.resampling == "cubic") resampleAlg = GRA_Cubic;
    else if (opts.resampling == "cubicspline") resampleAlg = GRA_CubicSpline;
    else if (opts.resampling == "lanczos") resampleAlg = GRA_Lanczos;
    else if (opts.resampling == "mode") resampleAlg = GRA_Mode;
    else if (opts.resampling == "max") resampleAlg = GRA_Max;
    else if (opts.resampling == "min") resampleAlg = GRA_Min;
    else if (opts.resampling == "med") resampleAlg = GRA_Med;
    else if (opts.resampling == "q1") resampleAlg = GRA_Q1;
    else if (opts.resampling == "q3") resampleAlg = GRA_Q3;

    // Create warped VRT
    GDALDataset *warpedDS = (GDALDataset *)GDALAutoCreateWarpedVRT(
        srcDS,
        opts.s_srs.empty() ? nullptr : opts.s_srs.c_str(),
        opts.t_srs.c_str(),
        resampleAlg,
        0.0,
        nullptr);

    if (!warpedDS)
    {
        std::cerr << "Warp VRT creation failed." << std::endl;
        GDALClose(srcDS);
        return false;
    }

    // Always set no-data value
    for (int i = 1; i <= warpedDS->GetRasterCount(); i++)
    {
        warpedDS->GetRasterBand(i)->SetNoDataValue(opts.nodata);
    }

    // Read output geotransform
    if (warpedDS->GetGeoTransform(geoTransform) != CE_None)
    {
        std::cerr << "Failed to get GeoTransform from warped dataset." << std::endl;
        GDALClose(srcDS);
        GDALClose(warpedDS);
        return false;
    }

    int xSize = warpedDS->GetRasterXSize();
    int ySize = warpedDS->GetRasterYSize();

    // Override output size if requested
    if (opts.out_cols > 0 && opts.out_rows > 0)
    {
        xSize = opts.out_cols;
        ySize = opts.out_rows;
    }

    // Allocate output buffer
    outData.resize(xSize * ySize);

    GDALRasterBand *band = warpedDS->GetRasterBand(1);
    if (band->RasterIO(GF_Read, 0, 0, xSize, ySize,
                       outData.data(), xSize, ySize,
                       GDT_Float32, 0, 0) != CE_None)
    {
        std::cerr << "Failed to read warped raster into memory." << std::endl;
        GDALClose(srcDS);
        GDALClose(warpedDS);
        return false;
    }

    GDALClose(srcDS);
    GDALClose(warpedDS);

    return true;
}

// Explicit template instantiations (so linker finds them)
template bool resample_dataset<int8_t>(const std::vector<int8_t>&, std::vector<int8_t>&, const GDALTranslateOpts&);
template bool resample_dataset<uint8_t>(const std::vector<uint8_t>&, std::vector<uint8_t>&, const GDALTranslateOpts&);
template bool resample_dataset<int16_t>(const std::vector<int16_t>&, std::vector<int16_t>&, const GDALTranslateOpts&);
template bool resample_dataset<uint16_t>(const std::vector<uint16_t>&, std::vector<uint16_t>&, const GDALTranslateOpts&);
template bool resample_dataset<int32_t>(const std::vector<int32_t>&, std::vector<int32_t>&, const GDALTranslateOpts&);
template bool resample_dataset<uint32_t>(const std::vector<uint32_t>&, std::vector<uint32_t>&, const GDALTranslateOpts&);
template bool resample_dataset<float>(const std::vector<float>&, std::vector<float>&, const GDALTranslateOpts&);
template bool resample_dataset<double>(const std::vector<double>&, std::vector<double>&, const GDALTranslateOpts&);

template bool saveToTiff<int8_t>(const std::string&, const std::vector<int8_t>&, int, int, const double[6], const std::string&, double);
template bool saveToTiff<uint8_t>(const std::string&, const std::vector<uint8_t>&, int, int, const double[6], const std::string&, double);
template bool saveToTiff<int16_t>(const std::string&, const std::vector<int16_t>&, int, int, const double[6], const std::string&, double);
template bool saveToTiff<uint16_t>(const std::string&, const std::vector<uint16_t>&, int, int, const double[6], const std::string&, double);
template bool saveToTiff<int32_t>(const std::string&, const std::vector<int32_t>&, int, int, const double[6], const std::string&, double);
template bool saveToTiff<uint32_t>(const std::string&, const std::vector<uint32_t>&, int, int, const double[6], const std::string&, double);
template bool saveToTiff<float>(const std::string&, const std::vector<float>&, int, int, const double[6], const std::string&, double);
template bool saveToTiff<double>(const std::string&, const std::vector<double>&, int, int, const double[6], const std::string&, double);

template bool reproject_dataset<int8_t>(const std::vector<int8_t> &,
                                     std::vector<int8_t> &,
                                     double[6],
                                     const GDALWarpOpts &);

template bool reproject_dataset<int16_t>(const std::vector<int16_t> &,
                                      std::vector<int16_t> &,
                                      double[6],
                                      const GDALWarpOpts &);

template bool reproject_dataset<int32_t>(const std::vector<int32_t> &,
                                      std::vector<int32_t> &,
                                      double[6],
                                      const GDALWarpOpts &);


template bool reproject_dataset<float>(const std::vector<float> &,
                                    std::vector<float> &,
                                    double[6],
                                    const GDALWarpOpts &);

template bool reproject_dataset<double>(const std::vector<double> &,
                                     std::vector<double> &,
                                     double[6],
                                     const GDALWarpOpts &);