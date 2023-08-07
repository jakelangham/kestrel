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


#ifndef RASTERDATA_H
#define RASTERDATA_H

#include <algorithm>
#include <string>
#include <iostream>
#include <vector>

// We prefer to include std::filesystem, but if it's not present
// try looking for boost::filesystem.
#if HAVE_FILESYSTEM
#include <filesystem>
#else 
#ifdef HAVE_BOOST_FILESYSTEM
#include "boost/filesystem.hpp"
#endif
#endif

#include "gdal_proxy.h"
#include "vrtdataset.h"
#include "cUTM.h"

// Cf #include <filesystem> above
#if HAVE_FILESYSTEM
namespace fs = std::filesystem;
#else 
#ifdef HAVE_BOOST_FILESYSTEM
namespace fs = boost::filesystem;
#endif
#endif

/**
 * @brief A container for data read from a geotiff file using GDAL.
 * 
 * @note RasterData contains data values and other essential data from the geotiff file, including the geotransform.
 * A geotransform defines an affine transformation from the 
 * raster image coordinate space (row, column), also known as (pixel, line),
 * to the georeferenced coordinate space (projected or geographic coordinates).
 * The pixel/line coordinates are from (0,0) at the top-left corner of the top-left pixel
 * to (x_size, y_size) at the bottom-right corner of the bottom-left pixel.
 * The pixel/line coordinate (0, 0) defines the origin of the raster.
 * The pixel/line location of the centre of the top-left pixel would therefore be (0.5, 0.5).
 * The x-axis (columns) is west--east, giving e.g. Easting or Longitude.
 * The y-axis (rows) is north--south for a north-up image, giving e.g. Northing or Latitude.
 * 
 * A transformation from image coordinate space (I, J) to georeferenced coordinate space (X, Y) is:
 *      X = x_origin + I * pixel_width + J * row_rotation
 *      Y = y_origin + I * column_rotation + J * pixel_height
 * 
 * @cite https://gdal.org/tutorials/geotransforms_tut.html
 */
class RasterData
{
    public:
        RasterData(const char *filename); // Constructor from filename
        RasterData(const char *filename, int xsize, int ysize, int xoff, int yoff); // Constructor from filename, size and offset (for reading section of file)
        bool containsNoData(void); // Check if Raster values contains a no data pixel.

        int read_success; ///< Flag for success read of the geotiff file: read_success = 1 if read is ok, otherwise 0.

        int32_t x_size;         ///< number of columns (size of the x-axis)
	    int32_t y_size;         ///< number of rows (size of the y-axis)
	    int32_t size;           ///< total number of elements ( = x_size * y_size)

	    double pixel_width;     ///< physical width of a pixel (x-resolution)
	    double pixel_height;    ///< physical height of a pixel (y-resolution)

	    double x_origin;        ///< physical x-coordinate of the top-left corner
	    double y_origin;        ///< physical y-coordinate of the top-left corner

        double row_rotation;    ///< typically zero for raster data with row/columns aligned with Easting/Longitude and Northing/Latitude
        double column_rotation; ///< typically zero for raster data with row/columns aligned with Easting/Longitude and Northing/Latitude

        int EPSG_code; ///< EPSG code determined from CRS of raster
        bool latlon; ///< Does the raster CRS uses lat-long coordinates?  True if yes, else false

        /// For convenience, the corner coordinates are stored as arrays (Easting, Northing) or (Longitude, Latitude)
	    double *NW;             ///< Northwest corner (top-left)
	    double *NE;             ///< Northeast corner (top-right)
	    double *SE;             ///< Southeast corner (bottom-right)
	    double *SW;             ///< Southwest corner (bottom-left)

        /// Raster data values as an array with size (x_size, y_size)
	    double *values;

        // Value assigned in the geotiff file for no data pixels.
	    double no_data;
};

/// @brief Store for useful information required by GDAL API
typedef struct
{
	int isFileOK;
	int nRasterXSize;
	int nRasterYSize;
	double adfGeoTransform[6];
	int nBlockXSize;
	int nBlockYSize;
	GDALDataType firstBandType;
	int *panHasNoData;
	double *padfNoDataValues;
} DatasetProperty;

/// @brief Store for useful information required by GDAL API
typedef struct
{
	GDALColorInterp colorInterpretation;
	GDALDataType dataType;
	GDALColorTableH colorTable;
	int bHasNoData;
	double noDataValue;
} BandProperty;

// Enumeration for GDAL geotransform attributes
#define GEOTRSFRM_TOPLEFT_X 0
#define GEOTRSFRM_WE_RES 1
#define GEOTRSFRM_ROTATION_PARAM1 2
#define GEOTRSFRM_TOPLEFT_Y 3
#define GEOTRSFRM_ROTATION_PARAM2 4
#define GEOTRSFRM_NS_RES 5

/**
 * @brief Get the Spatial Reference object for a GDAL Dataset
 * 
 * @param poDataset An open GDAL Dataset
 * @return OGRSpatialReference object
 * 
 * @cite https://gdal.org/tutorials/osr_api_tut.html
 */
OGRSpatialReference GetSpatialReference(GDALDataset *poDataset);

/**
 * @brief Create a collection of SRTM tiles covering a domain (passed in utm coordinates)
 * 
 * @param path full path of directory to contain output SRTM.vrt
 * @param SRTMpath full path of directory containing SRTM tiles
 * @param utmEPSG EPSG code of the UTM coordinate system
 * @param minE minimum Easting (i.e. most westward extent of domain) in UTM coordinates
 * @param maxE maximum Easting (i.e. most eastward extent of domain) in UTM coordinates
 * @param minN minimum Northing (i.e. most southward extent of domain) in UTM coordinates
 * @param maxN maximum Northing (i.e. most northward extent of domain) in UTM coordinates
 * 
 * @return std::vector<std::string> Filenames of SRTM tiles
 * 
 * @note The directory structure of SRTMpath must be in the form of a SRTM archive 
 * @note (e.g. SRTMpath/N01/, SRTMpath/N02/, ..., SRTMpath/S01/, etc.)
 * 
*/
std::vector<std::string> GatherSRTMtiles(const char *path, const char *SRTMpath, int *utmEPSG,
                     double *minE, double *maxE, double *minN, double *maxN);

/**
 * @brief Create a virtual geospatial file SRTM.vrt containing SRTM tiles from latitude--longitude bounding box
 * 
 * @param path full path of the directory to contain output SRTM.vrt
 * @param InputFilenames list of SRTM filenames (including full path) e.g. as created by GatherSRTMtiles
 */
void BuildSRTMVRT(const char *path, std::vector<std::string> InputFilenames);

/// The following function are bound to Fortran.  This dictates the input types passed by pointers
extern "C"
{
    /**
     * @brief Retrieve information about a raster from a geotiff file using GDAL
     * 
     * Using the GDAL library, a specified geotiff file is opened and attributes
     * of the RasterData structure are read.
     * 
     * @note For binding to Fortran, const char pointers are used.
     * 
     * @param cpath Full path of directory containing the geotiff file
     * @param ctifname Name of the geotiff file
     * @param raster RasterData class containing data from the geotiff
     * 
     */
	void GeoTiffInfo(const char *ctifname, RasterData *raster);

    /**
     * @brief Read a section of a geotiff file into an array using GDAL
     * 
     * @param cpath Full path of directory containing the geotiff file
     * @param ctifname Name of the geotiff file
     * @param raster RasterData struct containing data from the geotiff
     * @param xoff Offset (pixel defining top-left corner of array section) in x-direction 
     * @param yoff Offset (pixel defining top-left corner of array section) in y-direction 
     * @param xsize Width (number of pixels) of the array section
     * @param ysize Height (number of pixels) of the array section
     */
	void GeoTiffArraySectionRead(const char *cpath, const char *ctifname, RasterData *raster, int *xoff, int *yoff, int *xsize, int *ysize);

    /**
     * @brief Transform geotif file to a specified EPSG code using a prescribed EPSG code
     * 
     * @param path Full path of directory containing the SRTM file
     * @param ctifname Name of the geotiff file
     * @param epsgcode EPSG code, as an integer (e.g. 4326 for WGS84)
     * @param interp Interpolation choice, options are:
     *  0/default: nearest neighbour
     *          1: bilinear;
	 *          2: cubic;
	 *          3: cubic spline;
	 *          4: Lanczos;
	 *          5: average;
	 *          6: mode;
     */

    /**
     * @brief Create a vrt file containing a DEM and (optionally) SRTM tiles for required domain and resolution in WGS84/UTM coordinates
     * 
     * @param path Full path of directory that will contain the output files
     * @param SRTMpath Full path of directory containing the SRTM files
     * @param DEMfile Full path of the DEM file
     * @param utmEPSG EPSG code of the WGS84/UTM zone. This should be a five digit integer of the form 326** for northern hemisphere and 367** for southern hemisphere
     * @param embed Should the DEM file be embedded within SRTM data?  If true, an SRTM archive containing the required tiles is needed and SRTM.vrt is created to hold the SRTM data.
     * @param minx Minimum Easting (i.e. most westerly point) of the domain extent in WGS84/UTM coordinates corresponding to utmEPSG
     * @param maxx Maximum Easting (i.e. most easterly point) of the domain extent in WGS84/UTM coordinates corresponding to utmEPSG
     * @param miny Minimum Northing (i.e. most southerly point) of the domain extent in WGS84/UTM coordinates corresponding to utmEPSG
     * @param maxy Maximum Northing (i.e. most northerly point) of the domain extent in WGS84/UTM coordinates corresponding to utmEPSG
     * @param xres Resolution in Easting
     * @param yres Resolution in Northing
     * 
     * @note BuildDEMVRT_raster attempts to create the file DEM.vrt in directory path.
     * @note If required to embed DEM in SRTM tiles (i.e. embed = true) then a number of other vrt files are created that hold SRTM data
     * @note at intermediate steps.
     * @note The files are created using system calls to gdalwarp.  If these fail, this is a fatal error and the program exits with reports of
     * @note the failed gdalwarp call reported to std_err.
     * 
     * @cite https://gdal.org/programs/gdalwarp.html
     * @cite https://spatialreference.org/
     */
    void BuildDEMVRT_raster(const char *path,
				 const char *SRTMpath,
				 const char *DEMfile,
                 int *utmEPSG,
				 bool *embed,
				 double *minx, double *maxx, double *miny, double *maxy,
				 double *xres, double *yres);
    
    /**
     * @brief Create a vrt file containing a SRTM tiles for required domain and resolution in WGS84/UTM coordinates
     * 
     * @param path Full path of directory that will contain the output files
     * @param SRTMpath Full path of directory containing the SRTM files
     * @param utmEPSG EPSG code of the WGS84/UTM zone. This should be a five digit integer of the form 326** for northern hemisphere and 367** for southern hemisphere
     * @param embed Should the DEM file be embedded within SRTM data?  If true, an SRTM archive containing the required tiles is needed and SRTM.vrt is created to hold the SRTM data.
     * @param minx Minimum Easting (i.e. most westerly point) of the domain extent in WGS84/UTM coordinates corresponding to utmEPSG
     * @param maxx Maximum Easting (i.e. most easterly point) of the domain extent in WGS84/UTM coordinates corresponding to utmEPSG
     * @param miny Minimum Northing (i.e. most southerly point) of the domain extent in WGS84/UTM coordinates corresponding to utmEPSG
     * @param maxy Maximum Northing (i.e. most northerly point) of the domain extent in WGS84/UTM coordinates corresponding to utmEPSG
     * @param xres Resolution in Easting
     * @param yres Resolution in Northing
     * 
     * @note BuildDEMVRT_srtm attempts to create the file DEM.vrt in directory path.
     * @note A number of other vrt files are created that hold SRTM data at intermediate steps.
     * @note The files are created using system calls to gdalwarp.  If these fail, this is a fatal error and the program exits with reports of
     * @note the failed gdalwarp call reported to std_err.
     * 
     * @cite https://gdal.org/programs/gdalwarp.html
     * @cite https://spatialreference.org/
     * @cite https://lpdaac.usgs.gov/products/srtmgl1v003/
     */
    void BuildDEMVRT_srtm(const char *path,
				 const char *SRTMpath,
                 int *utmEPSG,
				 double *minx, double *maxx, double *miny, double *maxy,
				 double *xres, double *yres);
}

/**
 * @brief Check that the destination bounding box intersects the source bounding box
 * 
 * @param[in] psDP Information on source data contained in a DatasetProperty object
 * @param[in] we_res Resolution in Easting
 * @param[in] ns_res Resolution in Northing
 * @param[in] minX Minimum Easting (i.e. most westerly point)
 * @param[in] minY Minimum Northing (i.e. most southerly point)
 * @param[in] maxX Maximum Easting (i.e. most easterly point)
 * @param[in] maxY Maximum Northing (i.e. most northerly point)
 * @param[out] pnSrcXOff Offset in Easting of source dataset
 * @param[out] pnSrcYOff Offset in Northing of source dataset
 * @param[out] pnSrcXSize Size in Easting of source dataset
 * @param[out] pnSrcYSize Size in Easting of source dataset
 * @param[out] pnDstXOff Offset in Easting of destination dataset
 * @param[out] pnDstYOff Offset in Northing of destination dataset
 * @param[out] pnDstXSize Size in Easting of destination dataset
 * @param[out] pnDstYSize Size in Easting of destination dataset
 * @return true if the destination bounding box intersects the source bounding box
 * @return false if destination bounding box does not intersect the source bounding box
 * 
 * @note This is a utility function within BuildVRT
 */
bool GetSrcDstWin(DatasetProperty *psDP,
				 double we_res, double ns_res,
				 double minX, double minY, double maxX, double maxY,
				 int *pnSrcXOff, int *pnSrcYOff, int *pnSrcXSize, int *pnSrcYSize,
				 int *pnDstXOff, int *pnDstYOff, int *pnDstXSize, int *pnDstYSize);


/**
 * @brief Build a virtual file (vrt) of specified extend from list of files
 * 
 * @param[in] pszOutputFilename Name of output file
 * @param[in] InputFilenames List of input file names
 * @param[in] minX Minimum Easting (i.e. most westerly point)
 * @param[in] minY Minimum Northing (i.e. most southerly point)
 * @param[in] maxX Maximum Easting (i.e. most easterly point)
 * @param[in] maxY Maximum Northing (i.e. most northerly point)
 * @param pfnProgress GDALProgressFunc or NULL
 * @param pProgressData 
 * @return CPLErr GDAL error indicator
 * 
 * @cite https://gdal.org/drivers/raster/vrt.html
 */
CPLErr BuildVRT(const char *pszOutputFilename,
					std::vector<std::string> InputFilenames,
					double minX, double minY, double maxX, double maxY,
					GDALProgressFunc pfnProgress, void *pProgressData);

/**
 * @brief Create a zero-padded string from integer
 * 
 * @param[in] i Integer to convert to zero-padded string
 * @param[in] len Length of output string
 * @return std::string Zero-padded integer as a string
 */
std::string zero_pad(int i, int len);

/**
 * @brief Get name of longitude tile from integer (e.g. E012)
 * 
 * @param[in] i Integer tile longitude
 * @param[in] CAPS If true the E/W character is in caps, otherwise not
 * @return std::string tile longitude string
 */
std::string LongName(int i, bool CAPS=true);

/**
 * @brief Get name of latitude tile from integer (e.g. N01)
 * 
 * @param[in] i Integer tile latitude
 * @param[in] CAPS If true the E/W character is in caps, otherwise not
 * @return std::string tile latitude string
 */
std::string LatName(int i, bool CAPS=true);
#endif
