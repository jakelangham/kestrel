// This file is part of the Kestrel software for simulations
// of sediment-laden Earth surface flows.
//
// Version v1.1.1
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


// This is a cpp version of the python utm package
// Main use in LaharFlow is to determine utm zone from
// latitude, longitude of the domain.
#include "cUTM.h"
#include <stdexcept>
#include <stdlib.h>
#include <math.h>
#include <limits>
#include <iostream>
#include <sstream>
#include <vector>

/* Constructor */
proj_transformer::proj_transformer(int utm_code) {

    if (!valid_EPSG_utm(utm_code)) {
        std::cerr << "Failed to create transformation object." << std::endl;
        std::cerr << "UTM code " << utm_code << " is not valid." << std::endl;
        exit (EXIT_FAILURE);
    }
    
    utmEPSG = utm_code;

    C = proj_context_create();

    std::string source_crs = "EPSG:4326";
    std::string target_crs = "EPSG:" + std::to_string(utm_code);

    P = proj_create_crs_to_crs(C, source_crs.c_str(), target_crs.c_str(), NULL);
    if (0 == P) {
        std::cerr << "Failed to create transformation object." << std::endl;
        exit (EXIT_FAILURE);
    }

    /* This will ensure that the order of coordinates for the input CRS */
    /* will be longitude, latitude, whereas EPSG:4326 mandates latitude, longitude.*/
    /* More importantly, I think it means output will always be easting, northing.*/
    norm = proj_normalize_for_visualization(C, P);
    if (0 == norm) {
        std::cerr << "Failed to normalize transformation object." << std::endl;
        exit (EXIT_FAILURE);
    }
    proj_destroy(P);
    P = norm;
}

/* Destructor */
proj_transformer::~proj_transformer() {
    proj_destroy(P);
    proj_context_destroy(C);
}

extern "C" {

// C constructor
PT* proj_transformer__new(int utm_code) {
    return new proj_transformer(utm_code);
}

// C destructor
void proj_transformer__delete(PT* self) {
    delete self;
}

}

/* Method: transform wgs84 latitude--longitude to UTM Easting--Northing*/
std::vector<double> proj_transformer::wgs84_to_utm(double latitude, double longitude) const{

    /* Change input lat-lon to a proj_coord.
       Note, as we have used proj_normalize_for_visualization(), the order
       of coordinates is longitude, latitude, for values in degrees. */
    PJ_COORD a = proj_coord(longitude, latitude, 0, 0);

    /* transform to required utm zone */
    PJ_COORD b = proj_trans(P, PJ_FWD, a);

    std::vector<double> EN;
    EN = {b.enu.e, b.enu.n};

    return EN;
}
// C interface
double* proj_transformer__wgs84_to_utm(const PT* self, double latitude, double longitude) {
    std::vector<double> EN;
    EN = self->wgs84_to_utm(latitude, longitude);
    static double a[2];
    a[0] = EN[0];
    a[1] = EN[1];
    return a;
}

/* Method: transform UTM Easting--Northing to wgs84 latitude--longitude */
std::vector<double> proj_transformer::utm_to_wgs84(double easting, double northing) const{
    
    /* Change input Easting-Northing to a proj_coord.
       Note, as we have used proj_normalize_for_visualization(), the order
       of coordinates is Easting, Northing, for values in meters. */
    PJ_COORD a = proj_coord(easting, northing, 0, 0);

    /* transform to required utm zone */
    PJ_COORD b = proj_trans(P, PJ_INV, a);

    std::vector<double> latlon;
    latlon = {b.lp.phi, b.lp.lam};
    
    return latlon;
}
// C interface
double* proj_transformer__utm_to_wgs84(const PT* self, double easting, double northing) {
    std::vector<double> latlon;
    latlon = self->utm_to_wgs84(easting, northing);
    static double a[2];
    a[0] = latlon[0];
    a[1] = latlon[1];
    return a;
}

int latlon_to_zone_number(double latitude, double longitude) {
    if ((56.0 <= latitude) && (latitude < 64.0) && (3.0 <= longitude) && (longitude < 12.0)) { return 32; }
    
    if ((72.0 <= latitude) && (latitude <= 84.0) && (longitude >= 0.0)) {
        if (longitude < 9.0)  { return 31; }
        if (longitude < 21.0) { return 33; }
        if (longitude < 33.0) { return 35; }
        if (longitude < 42.0) { return 37; }
    }
      
    int zone_number = (int) (((longitude + 180.0) / 6.0) + 1);
    return zone_number;
}


int zone_number_to_central_longitude(int zone_number) {
    return (zone_number - 1) * 6 - 180 + 3;
}

int latlon_to_utm_epsg(double latitude, double longitude) {
      int offset = (int) round((183.0 + longitude)/6.0);
      if (latitude > 0) { return 32600+offset; }
      else { return 32700+offset; }
}

bool valid_EPSG_utm(int utm_code) {
    if (((unsigned) (utm_code-32601) <= (32660-32601)) || ((unsigned) (utm_code-32701) <= (32760-32701))) { return true; }
    else { return false; }
}
