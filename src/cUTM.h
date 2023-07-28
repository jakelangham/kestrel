// This file is part of the Kestrel software for simulations
// of sediment-laden Earth surface flows.
//
// Version 1.0
//
// Copyright 2023 Mark J. Woodhouse, Jake Langham, (University of Bristol).
//
// This program is free software: you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free 
// Software Foundation, either version 3 of the License, or (at your option) 
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT 
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
// more details.
//
// You should have received a copy of the GNU General Public License along with 
// this program. If not, see <https://www.gnu.org/licenses/>. 


// cUTM.h
#ifndef CUTM_H
#define CUTM_H
#include <stdexcept>
#include <math.h>
#include <limits>
#include <iostream>
#include <sstream>
#include <vector>
#include "proj.h"

const std::string ZONE_LETTERS = "CDEFGHJKLMNPQRSTUVWXX";

/**
 * @brief Transforms coordinates between WGS84 and WGS84 / UTM code
 * 
 * @note Creating a proj_transformer object using a WGS84 / UTM EPSG code
 * @note allows forward (to WGS84) and backward (from WGS84) transformations
 * @note to be performed.
 * 
 * @cite https://spatialreference.org/
 */
class proj_transformer {
    public:
        /**
         * @brief Construct a new proj transformer object from a EPSG utm code
         * 
         * @note EPSG codes for WGS84 / UTM zone are five-digit integers.
         * @note Zones in the northern hemisphere have the form 326**.
         * @note Zones in the southern hemisphere have the form 327**
         * 
         * @param[in] utm_code A five digit EPSG code for WGS84 / UTM zone
         */
        proj_transformer(int utm_code);

        /**
         * @brief Destroy the proj transformer object
         * 
         */
        ~proj_transformer();
        
        /**
         * @brief The EPSG code
         * This is a five-digit integer
         * in the range [32601,32660] in Northern hemisphere
         * or the range [32701,32760] in Southern hemisphere
         */
        int utmEPSG;
        
        /**
         * @brief Convert a WGS84 coordinate (latitude, longitude) to the required WGS84 / UTM zone
         * 
         * @param[in] latitude WGS84 latitude in decimal degrees
         * @param[in] longitude WGS84 longitude in decimal degrees
         * @return std::vector<double> UTM coordinates as a vector with (Easting, Northing)
         */
        std::vector<double> wgs84_to_utm(double latitude, double longitude) const;

        /**
         * @brief Convert a UTM coordinate (easting, northing) to WGS84 (latitude, longitude).
         * 
         * @param[in] easting UTM easting in metres
         * @param[in] northing UTM northing in metres
         * @return std::vector<double>  WGS84 coordinates as a vector with (longitude, latitude)
         */
        std::vector<double> utm_to_wgs84(double easting, double northing) const;
    private:
        PJ_CONTEXT *C;
        PJ *P;
        PJ *norm;
};

/**
 * @brief Determine if a integer is a valid EPSG code for a WGS84 / UTM zone
 * 
 * @param[in] utm_code An integer code
 * @return true if utm_code is a valid EPSG code for a WGS84 / UTM zone
 * @return false otherwise
 */
bool valid_EPSG_utm(int utm_code);

extern "C"
{

    /**
     * @brief C-compatible proj_transformer
     * 
     */
    class proj_transformer;
    typedef proj_transformer PT;

    /**
     * @brief C-compatible constructor
     * 
     * @param utm_code A five digit EPSG code for WGS84 / UTM zone
     * @return PT* pointer to proj_transformer object
     */
    PT* proj_transformer__new(int utm_code);

    /**
     * @brief C-compatible destructor
     * 
     * @param self proj_transformer to destroy
     */
    void proj_transformer__delete(PT* self);

    /**
     * @brief C-compatible interface to proj_transformer method wgs84_to_utm 
     * 
     * @param self proj_transformer object
     * @param latitude WGS84 latitude in decimal degrees
     * @param longitude WGS84 longitude in decimal degrees
     * @return double* two-element array containing UTM coordinates as a vector with (Easting, Northing)
     */
    double* proj_transformer__wgs84_to_utm(const PT* self, double latitude, double longitude);

    /**
     * @brief C-compatible interface to proj_transformer method utm_to_wgs84
     * 
     * @param self proj_transformer object
     * @param easting UTM easting in metres
     * @param northing UTM northing in metres
     * @return double* two-element array containing WGS84 coordinates as a vector with (longitude, latitude)
     */
    double* proj_transformer__utm_to_wgs84(const PT* self, double easting, double northing);

    /**
     * @brief Get UTM zone number from WGS84 latitude, longitude coordinates
     * 
     * @param[in] latitude WGS84 latitude in decimal degrees
     * @param[in] longitude WGS84 longitude in decimal degrees
     * @return int UTM zone number
     */
    int latlon_to_zone_number(double latitude, double longitude);

    /**
     * @brief Get central longitude from WGS84 / UTM zone number
     * 
     * @param[in] zone_number WGS84 / UTM zone number
     * @return int central longitude in decimal degrees
     */
    int zone_number_to_central_longitude(int zone_number);

    /**
     * @brief Get EPSG code for WGS84 / UTM zone from WGS84 latitude, longitude coordinates
     * 
     * @param[in] latitude WGS84 latitude in decimal degrees
     * @param[in] longitude WGS84 longitude in decimal degrees
     * @return int UTM zone number
     */
    int latlon_to_utm_epsg(double latitude, double longitude);
}

#endif
