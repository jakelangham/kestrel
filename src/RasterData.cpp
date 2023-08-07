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


#include "RasterData.h"

// We prefer to include std::filesystem, but if it's not present
// try looking for boost::filesystem.
#if HAVE_FILESYSTEM
#include <filesystem>
#else 
#ifdef HAVE_BOOST_FILESYSTEM
#include "boost/filesystem.hpp"
#endif
#endif

RasterData::RasterData(const char *filename) {
	// Register GDAL drivers
	GDALDataset *poDataset;
	GDALAllRegister();

	// Open file using GDAL in read-only mode
	poDataset = (GDALDataset *)GDALOpen(filename, GA_ReadOnly);

	// Check if file opens successfully
	if (poDataset == NULL) // File could not be read
	{
        read_success = 0;
	}
	else // Successful file open
	{
        read_success = 1; // Report successful read.

        x_size = poDataset->GetRasterXSize();  	// Read raster x-size using GDAL
	    y_size = poDataset->GetRasterYSize();  	// Read raster y-size using GDAL
	    size = x_size * y_size; // Total number of data values = x_size * y_size

        // Get EPSG projection reference
	    OGRSpatialReference sref = GetSpatialReference(poDataset);
        latlon = sref.EPSGTreatsAsLatLong(); // Check whether CRS uses lat-long coordinates
		
	    // Retrieve EPSG code and store
        try {
	        if (strcmp(sref.GetAuthorityName(NULL), "EPSG")==0) {
	    	    EPSG_code = atoi(sref.GetAuthorityCode(NULL));
	        }
            else {
                throw (poDataset->GetFileList());
            }
        }
        catch (char **fileList) {
            std::cerr << "Error: no EPSG authority in " << fileList << std::endl;
            std::cerr << "Check files and ensure they are properly georeferenced" << std::endl;
            exit(EXIT_FAILURE);
        }

        try {
	        // Get geotransform using GDAL
	        double adfGeoTransform[6]; // Create array to contain geotransform
	        CPLErr err = poDataset->GetGeoTransform(adfGeoTransform);
	        if (err == CE_None) // Read geotransform data using GDAL is successful
	        {
		        x_origin = adfGeoTransform[GEOTRSFRM_TOPLEFT_X]; // x-coordinate of upper-left corner or the upper-left pixel
                pixel_width = adfGeoTransform[GEOTRSFRM_WE_RES]; // west--east pixel resolution
                row_rotation = adfGeoTransform[GEOTRSFRM_ROTATION_PARAM1]; // row rotation, typically zero
		        y_origin = adfGeoTransform[GEOTRSFRM_TOPLEFT_Y]; // y-coordinate of upper-left corner or the upper-left pixel
		        column_rotation = adfGeoTransform[GEOTRSFRM_ROTATION_PARAM2]; // column rotation, typically zero
		        pixel_height = adfGeoTransform[GEOTRSFRM_NS_RES]; // north--south pixel resolution (negative value for a north-up image)

		        // Add corner coordinates
                NW = new double[2];
                NW[0] = y_origin;
                NW[1] = x_origin;

		        NE = new double[2];
                NE[0] = y_origin;
                NE[1] = x_origin + x_size * pixel_width;

                SW = new double[2];
                SW[0] = y_origin + y_size * pixel_height;
                SW[1] = x_origin;

                SE = new double[2];
	            SE[0] = y_origin + y_size * pixel_height;
	            SE[1] = x_origin + x_size * pixel_width;

                if (latlon) { // if lat-long coordinates are used, make sure the longitude is in [0, 360]
		            NW[1] = fmod(360. + NW[1], 360.);
		            NE[1] = fmod(360. + NE[1], 360.);
		            SW[1] = fmod(360. + SW[1], 360.);
		            SE[1] = fmod(360. + SE[1], 360.);
                }
            }
            else {
                throw (err);
            }
        }
        catch (CPLErr err) {
            std::cerr << "Error: could not get geotransform.  GDAL error " << err << std::endl;
            std::cerr << "Check files and ensure they are properly georeferenced" << std::endl;
            exit(EXIT_FAILURE);
        }

        try {
            GDALRasterBand *poBand;
		    poBand = poDataset->GetRasterBand(1);

            // Get no-data value
            int pbSuccess = 0;
		    no_data = poBand->GetNoDataValue(&pbSuccess);
            if (pbSuccess == 0) {
                throw (pbSuccess);
            }
        }
        catch (int pbSuccess) {
            std::cerr << "Error: could not get no data value for raster band." << std::endl;
            std::cerr << "Check files and ensure they have a no data value correctly set." << std::endl;
            exit(EXIT_FAILURE);
        }

	}
	GDALClose(poDataset);	// Close the geotiff
	return;
}

RasterData::RasterData(const char *filename, int xsize, int ysize, int xoff, int yoff) {
	// Register GDAL drivers
	GDALDataset *poDataset;
	GDALAllRegister();

	// Open file using GDAL in read-only mode
	poDataset = (GDALDataset *)GDALOpen(filename, GA_ReadOnly);

	// Check if file opens successfully
	if (poDataset == NULL) // File could not be read
	{
        read_success = 0;
	}
	else // Successful file open
	{
        read_success = 1; // Report successful read.

        x_size = xsize;  	// Set x-size from input
	    y_size = ysize;  	// Set y-size from input
	    size = x_size * y_size; // Total number of data values = x_size * y_size

        // Get EPSG projection reference
	    OGRSpatialReference sref = GetSpatialReference(poDataset);
        latlon = sref.EPSGTreatsAsLatLong(); // Check whether CRS uses lat-long coordinates
		
	    // Retrieve EPSG code and store
        try {
	        if (strcmp(sref.GetAuthorityName(NULL), "EPSG")==0) {
	    	    EPSG_code = atoi(sref.GetAuthorityCode(NULL));
	        }
            else {
                throw (poDataset->GetFileList());
            }
        }
        catch (char **fileList) {
            std::cerr << "Error: no EPSG authority in " << fileList << std::endl;
            std::cerr << "Check files and ensure they are properly georeferenced" << std::endl;
            exit(EXIT_FAILURE);
        }

        try {
	        // Get geotransform using GDAL
	        double adfGeoTransform[6]; // Create array to contain geotransform
	        CPLErr err = poDataset->GetGeoTransform(adfGeoTransform);
	        if (err == CE_None) // Read geotransform data using GDAL is successful
	        {
		        x_origin = adfGeoTransform[GEOTRSFRM_TOPLEFT_X] + xoff * adfGeoTransform[GEOTRSFRM_WE_RES]; // x-coordinate of upper-left corner or the upper-left pixel
                pixel_width = adfGeoTransform[GEOTRSFRM_WE_RES]; // west--east pixel resolution
                row_rotation = adfGeoTransform[GEOTRSFRM_ROTATION_PARAM1]; // row rotation, typically zero
		        y_origin = adfGeoTransform[GEOTRSFRM_TOPLEFT_Y] + yoff * adfGeoTransform[GEOTRSFRM_NS_RES]; // y-coordinate of upper-left corner or the upper-left pixel
		        column_rotation = adfGeoTransform[GEOTRSFRM_ROTATION_PARAM2]; // column rotation, typically zero
		        pixel_height = adfGeoTransform[GEOTRSFRM_NS_RES]; // north--south pixel resolution (negative value for a north-up image)

		        // Add corner coordinates
                NW = new double[2];
                NW[0] = y_origin;
                NW[1] = x_origin;

		        NE = new double[2];
                NE[0] = y_origin;
                NE[1] = x_origin + x_size * pixel_width;

                SW = new double[2];
                SW[0] = y_origin + y_size * pixel_height;
                SW[1] = x_origin;

                SE = new double[2];
	            SE[0] = y_origin + y_size * pixel_height;
	            SE[1] = x_origin + x_size * pixel_width;

                if (latlon) { // if lat-long coordinates are used, make sure the longitude is in [0, 360]
		            NW[1] = fmod(360. + NW[1], 360.);
		            NE[1] = fmod(360. + NE[1], 360.);
		            SW[1] = fmod(360. + SW[1], 360.);
		            SE[1] = fmod(360. + SE[1], 360.);
                }
            }
            else {
                throw (err);
            }
        }
        catch (CPLErr err) {
            std::cerr << "Error: could not get geotransform.  GDAL error " << err << std::endl;
            std::cerr << "Check files and ensure they are properly georeferenced" << std::endl;
            exit(EXIT_FAILURE);
        }

        try {
            GDALRasterBand *poBand;
		    poBand = poDataset->GetRasterBand(1);

            // Get no-data value
            int pbSuccess = 0;
		    no_data = poBand->GetNoDataValue(&pbSuccess);
            if (pbSuccess == 0) {
                throw (pbSuccess);
            }
        }
        catch (int pbSuccess) {
            std::cerr << "Error: could not get no data value for raster band." << std::endl;
            std::cerr << "Check files and ensure they have a no data value correctly set." << std::endl;
            exit(EXIT_FAILURE);
        }

	}
	GDALClose(poDataset);	// Close the geotiff
	return;
}

OGRSpatialReference GetSpatialReference(GDALDataset *poDataset)
{
	// Get projection reference
	const char *proj = GDALGetProjectionRef(poDataset);
	// Initialize OGRSpatialReference using projection reference
	OGRSpatialReference sref;
	OGRErr serr = sref.importFromWkt(proj);
	// Report error
	// TODO: should this kill?
	if (serr != OGRERR_NONE) { // Error importing spatial reference
		std::cerr << "GDAL error getting spatial reference " << serr << std::endl;
	}
	return sref;
}

void GeoTiffInfo(const char *ctifname, RasterData *raster)
{

    *raster = RasterData(ctifname);
	return;
}

void GeoTiffArraySectionRead(const char *cpath, const char *ctifname, RasterData *raster, int *xoff, int *yoff, int *xsize, int *ysize)
{
	// Convert input char to string for ease of concatenation
	std::string pathStr(cpath);
	std::string ctifnameStr(ctifname);

	// Register drivers
	GDALAllRegister();

	// Create full path of file to open
    std::string sourceFilename = pathStr + ctifnameStr;
	fs::path sourceFile { sourceFilename };
	char *inFile = new char[sourceFilename.length() + 1];
	strcpy(inFile, sourceFilename.c_str()); // GDAL required char for file name

    *raster = RasterData(inFile, *xsize, *ysize, *xoff, *yoff);

    if (raster->read_success) {
        
        try {
			GDALDataset *poDataset;
            // Open file
	        poDataset = (GDALDataset *)GDALOpen(inFile, GA_ReadOnly);

	        if (poDataset == NULL) // Failed to open file
	        {
                throw;
	        }

			// Fetch a raster band
			GDALRasterBand *poBand;
			poBand = poDataset->GetRasterBand(1);

			// Retrieve block sizes
			int nBlockXSize, nBlockYSize;
			poBand->GetBlockSize(&nBlockXSize, &nBlockYSize);

			// Allocate array for values
			raster->values = (double *)CPLMalloc(sizeof(double) * (*xsize) * (*ysize));

			// Get no-data value
			raster->no_data = poBand->GetNoDataValue();

			// Read values in the array section with origin at (xoff, yoff) and with size (xsize, ysize)
			CPLErr read_err = poBand->RasterIO(GF_Read, *xoff, *yoff, *xsize, *ysize, raster->values, *xsize, *ysize, GDT_Float64, 0, 0);

			if (read_err != CE_None) // Error reading
			{
				CPLError(read_err, CPLE_AppDefined, "Error reading file %s.", inFile);
				raster->read_success = 0; // Report failed read -- note this exception is handled on the Fortran side.
			}
			GDALClose(poDataset);
        }
        catch (...) {
            std::cerr << "Error: could not open raster " << sourceFilename << std::endl;
            exit(EXIT_FAILURE);
        }		
	}
	return;
}

std::vector<std::string> GatherSRTMtiles(const char *path, const char *SRTMpath, int *utmEPSG,
                     double *minE, double *maxE, double *minN, double *maxN)
{

    proj_transformer PT = proj_transformer(*utmEPSG);

    /* Get latitude longitude of corners */
    /* Note: these are ordered as latitude-longitude. */
    std::vector<double> NE = PT.utm_to_wgs84(*maxE, *maxN);
    std::vector<double> SE = PT.utm_to_wgs84(*maxE, *minN);
    std::vector<double> SW = PT.utm_to_wgs84(*minE, *minN);
    std::vector<double> NW = PT.utm_to_wgs84(*minE, *maxN);
    
    /* Get extent in latitude--longitude coordinates */
    double minlat = std::min(std::min(SW[0], SE[0]), std::min(NW[0], NE[0]));
    double maxlat = std::max(std::max(SW[0], SE[0]), std::max(NW[0], NE[0]));
    double minlon = std::min(std::min(SW[1], NW[1]), std::min(SE[1], NE[1]));
    double maxlon = std::max(std::max(SW[1], NW[1]), std::max(SE[1], NE[1]));
    
    /* Cast to integer for tile labels */
    int minLatFloor = static_cast<int>(floor(minlat));
    int maxLatFloor = static_cast<int>(floor(maxlat));
    int minLonFloor = static_cast<int>(floor(minlon));
    int maxLonFloor = static_cast<int>(floor(maxlon));

	std::string SRTMdir(SRTMpath); // Store path to SRTM files

    std::vector<std::string> srtm_ext = {".hgt", ".tif"};

	/* Create container for SRTM tile filenames */
	std::vector<std::string> InputFilenames;

	int latCount = 1; // Counter for latitude tiles
	int tileCount = 0; // Counter for total number of tiles
	for (int latIndex = minLatFloor; latIndex < maxLatFloor + 1; latIndex++) // loop over latitude tiles
	{
		int longIndex = minLonFloor; // index for longitude
		int longCount = 1; // Counter for longitude tiles
		while ((longIndex % 360) != ((maxLonFloor + 1) % 360)) // loop until we reach maximum longitude, but take care about wrapping
		{
			std::string latID = LatName(latIndex);		// ID of latitude e.g. N01
			std::string lonID = LongName(longIndex);	// ID of longitude e.g. E123

			// Construct SRTM filename
			std::string SRTMfile_base = SRTMdir + latID + "/" + latID + lonID; // + ".hgt";
            std::string SRTMfile;

			// Look for file
            std::cout << " Looking for ";
            bool found_file = false;
            for (int i=0; i<2; i++) {
                if (i==1) {
                    std::string latID = LatName(latIndex, false);
				    std::string lonID = LongName(longIndex, false);
                    SRTMfile_base = SRTMdir + latID + "/" + latID + lonID;
                }
                for (auto& ext : srtm_ext) {
                    SRTMfile = SRTMfile_base + ext;
                    std::cout << SRTMfile << " ... ";
                    if (fs::exists(SRTMfile)) { 
                        std::cout << "found it!" << std::endl;
                        found_file = true;
                        break;
                    }
                }
                if (found_file) break;
            }
            if (! found_file)
			{

				std::cout << " did not find SRTM file" << std::endl;
                std::cerr << " SRTM file not found" << std::endl;
				exit(EXIT_FAILURE);
			}

			// If we got here, the SRTM file was found, so add to list of input files
			InputFilenames.push_back(SRTMfile);

			// Advance counters
			longIndex++;
			longIndex = longIndex % 360; // Take care of wrapping
			longCount++;
			tileCount++;
		}
		latCount++;
	}
	return InputFilenames;
}

void BuildSRTMVRT(const char *path, std::vector<std::string> InputFilenames)
{
	// initialize Ouput file, progress bar, extents
	const char *pszOutputFilename = NULL;
	GDALProgressFunc pfnProgress = NULL;
	double xmin = 0, ymin = 0, xmax = 0, ymax = 0;

	// Register gdal drivers
	GDALAllRegister();

	// Convert input paths to strings
	std::string OutputDir(path);

	// Create full path of output file
	std::string OutputFilename = OutputDir + "/SRTM.vrt";
	pszOutputFilename = OutputFilename.c_str(); // Convert to character string

	// Build vrt file using GDAL
	BuildVRT(pszOutputFilename, InputFilenames,
				 xmin, ymin, xmax, ymax,
				 pfnProgress, NULL);

	return;
}

void BuildDEMVRT_raster(const char *path,
				 const char *SRTMpath,
				 const char *DEMfile,
                 int *utmEPSG,
				 bool *embed,
				 double *minE, double *maxE, double *minN, double *maxN,
				 double *xres, double *yres)
{
	fs::path OutputDir = path;
	fs::path SRTMdir = SRTMpath;
	fs::path DEMfullfile = DEMfile;

	fs::path OutputFile = OutputDir / fs::path("DEM.vrt");

	RasterData input_dem = RasterData(DEMfile);
	// Get info on input DEM
	// GeoTiffInfo(DEMfile, &input_dem);

    proj_transformer PT = proj_transformer(*utmEPSG);

    std::string dem_to_warp;

    if (*embed) {
        // Create SRTM.vrt that gathers the required SRTM tiles
        std::vector<std::string> InputSRTMfiles = GatherSRTMtiles(path, SRTMpath, utmEPSG, minE, maxE, minN, maxN);
		BuildSRTMVRT(path, InputSRTMfiles);

        // Name virtual file for SRTM_EPSG_{X}.vrt
        std::string vrt_srtm_name = std::string("SRTM_EPSG_") + std::to_string(input_dem.EPSG_code) + std::string(".vrt");
        fs::path vrt_srtm_path = OutputDir / fs::path(vrt_srtm_name);

        // Warp SRTM_vrt to the same CRS as input DEM
        std::ostringstream gdal_exec_str_SRTM_wrp;
        gdal_exec_str_SRTM_wrp << "gdalwarp -overwrite -of VRT -et 0 -r cubic -ot Float64 "; // Start creation of warp string
        // set resolution
        std::ostringstream tr_str_SRTM_wrp; 
        tr_str_SRTM_wrp.precision(17); // hold resolution at sufficient precision
        tr_str_SRTM_wrp << std::fixed;
        tr_str_SRTM_wrp << "-tr " << input_dem.pixel_width << " " << input_dem.pixel_height; // set resolution to match DEM file
        gdal_exec_str_SRTM_wrp << tr_str_SRTM_wrp.str(); // Append resolution to gdal executable
        gdal_exec_str_SRTM_wrp << " -t_srs EPSG:" << std::to_string(input_dem.EPSG_code) << " "; // set CRS to match DEM
        gdal_exec_str_SRTM_wrp << " -q "; // Quiet gdalwarp
        gdal_exec_str_SRTM_wrp << path << "SRTM.vrt " << vrt_srtm_path.string(); // Create warped output file
        // Report gdal commandline statement
        
        // Run gdal command, report if fail
        int gdal_ret_SRTM_wrp = system((gdal_exec_str_SRTM_wrp).str().c_str());
        if (gdal_ret_SRTM_wrp != 0) {
            std::cerr << "Error: could not create SRTM_utm.vrt" << std::endl;
            std::cerr << "The failed call to gdalwarp is:" << std::endl;
            std::cerr << gdal_exec_str_SRTM_wrp.str() << std::endl;
            exit(EXIT_FAILURE);
        }

        /// Create a temporary VRT file containing DEM, called DEM_EPSG_{X}.vrt
	    /// We do this to ensure the settings are the same as vrt_srtm_name
	    /// so we can mosiac and reproject
	    std::string vrt_dem_name = fs::path(DEMfile).stem().string() + std::string(".vrt");
        fs::path vrt_dem_path = OutputDir / fs::path(vrt_dem_name);
	    std::vector<std::string> dem_list = {std::string(DEMfile)};
	    BuildVRT(vrt_dem_path.c_str(), dem_list,
	    	input_dem.SW[0], input_dem.SW[1], input_dem.NE[0], input_dem.NE[1],
	    	NULL, NULL
	    );

        /// Now we create a VRT file containing vrt_srtm_path and vrt_dem_path
        /// These are in the same CRS so can be mosaiced in a warp
        std::string vrt_combined_name = std::string("DEM_EPSG_") + std::to_string(input_dem.EPSG_code) + std::string(".vrt");
	    fs::path vrt_combined_path = OutputDir / fs::path(vrt_combined_name);
	    std::vector<std::string> vrt_list = {vrt_srtm_path.string(), vrt_dem_path.string()};
	    BuildVRT(vrt_combined_path.c_str(), vrt_list,
		    *minE, *minN, *maxE, *maxN,
		    NULL, NULL
	    );

        /// The VRT created above is the file that will be warped.
        dem_to_warp = vrt_combined_path.string();

    }
    else {
        /// We can directly warp the input DEM geotif file into a VRT in the required UTM coordinate system
        dem_to_warp = DEMfullfile.string();
    }

	/// Finally, we warp the dem_to_warp file
	/// to the required CRS and resolution.
	/// The ouput is DEM.vrt.
	std::ostringstream gdal_exec_str;
	gdal_exec_str << "gdalwarp -overwrite -of VRT -et 0 -r cubic -ot Float64 "; // use cubic interpolation and set output format to 64-bit float
	// set resolution
	std::ostringstream tr_str;
	tr_str.precision(17); // Store at sufficient resolution
	tr_str << std::fixed;
	tr_str << "-tr " << *xres << " " << *yres; // set resolution to match required grid resolution
	// Set extent
	std::ostringstream te_str;
	te_str.precision(17); // Store at sufficient resolution
	te_str << " -te " << *minE << " " << *minN << " " << *maxE << " " << *maxN << " "; // set extent to match required domain
	// Append to gdal command line string
	gdal_exec_str << tr_str.str();
	gdal_exec_str << te_str.str();
	gdal_exec_str << " -t_srs EPSG:" << std::to_string(*utmEPSG) << " "; // Add CRS from required utmEPSG
    gdal_exec_str << " -q "; // Quiet gdalwarp
	gdal_exec_str << dem_to_warp << " "; // Add DEM file
	gdal_exec_str << OutputFile.string(); // Set output name
	// Execute gdal command, report if fail
	int gdal_ret = system((gdal_exec_str).str().c_str());
	if (gdal_ret != 0) {
		std::cerr << "Error: could not create DEM.vrt" << std::endl;
        std::cerr << "The failed call to gdalwarp is:" << std::endl;
        std::cerr << gdal_exec_str.str() << std::endl;
		exit(EXIT_FAILURE);
	}

    // Get info for DEM.vrt
    // RasterData out_dem = RasterData(OutputFile.c_str());
    // GeoTiffInfo(OutputFile.c_str(), &out_dem);

    // double dem_minE = std::min(out_dem.SW[1], out_dem.NW[1]);
    // double dem_maxE = std::max(out_dem.SE[1], out_dem.NE[1]);
    // double dem_minN = std::min(out_dem.SW[0], out_dem.SE[0]);
    // double dem_maxN = std::max(out_dem.NW[0], out_dem.NE[0]);

    // // Check DEM.vrt covers the domain of the 
    // if ((*minE < dem_minE) or (*maxE > dem_maxE) or (*minN < dem_minN) or (*maxN > dem_maxN))
    // { 
    //     std::cerr << "Error: domain extends beyond DEM.vrt" << std::endl;
    //     std::cerr << "either decrease domain size or embed raster in SRTM" << std::endl;
    //     exit(EXIT_FAILURE);
    // }
	return;
}

void BuildDEMVRT_srtm(const char *path,
				 const char *SRTMpath,
                 int *utmEPSG,
				 double *minE, double *maxE, double *minN, double *maxN,
				 double *xres, double *yres)
{
	fs::path OutputDir = path;
	fs::path SRTMdir = SRTMpath;

	fs::path OutputFile = OutputDir / fs::path("DEM.vrt");

    std::string dem_to_warp;

    // Create SRTM.vrt that gathers the required SRTM tiles
    std::vector<std::string> InputSRTMfiles = GatherSRTMtiles(path, SRTMpath, utmEPSG, minE, maxE, minN, maxN);
	BuildSRTMVRT(path, InputSRTMfiles);

    fs::path srtm_vrt = OutputDir / fs::path("SRTM.vrt");

	/// Warp the SRTM.vrt file
	/// to the required CRS and resolution.
	/// The ouput is DEM.vrt.
	std::ostringstream gdal_exec_str;
	gdal_exec_str << "gdalwarp -overwrite -of VRT -et 0 -r cubic -ot Float64 ";
	// set resolution
	std::ostringstream tr_str;
	tr_str.precision(17); // Store at sufficient resolution
	tr_str << std::fixed;
	tr_str << "-tr " << *xres << " " << *yres; // Set resolution to match required grid resolution
	// Set extent
	std::ostringstream te_str;
	te_str.precision(17); // Store at sufficient resolution
	te_str << " -te " << *minE << " " << *minN << " " << *maxE << " " << *maxN << " "; // set extent to match required domain
	// Append to gdal command line string
	gdal_exec_str << tr_str.str();
	gdal_exec_str << te_str.str();
	
	gdal_exec_str << " -t_srs EPSG:" << std::to_string(*utmEPSG) << " "; // Add CRS from required utmEPSG
    gdal_exec_str << " -q "; // Quiet gdalwarp
	gdal_exec_str << srtm_vrt.string() << " "; // Add SRTM.vrt
	gdal_exec_str << OutputFile.string(); // Set output name
	
	// Execute gdal command, report if fail
	int gdal_ret = system((gdal_exec_str).str().c_str());
	if (gdal_ret != 0) {
		std::cerr << "Error: could not create DEM.vrt" << std::endl;
        std::cerr << "The failed call to gdalwarp is:" << std::endl;
        std::cerr << gdal_exec_str.str() << std::endl;
		exit(EXIT_FAILURE);
	}
	return;
}

bool GetSrcDstWin(DatasetProperty *psDP,
				 double we_res, double ns_res,
				 double minX, double minY, double maxX, double maxY,
				 int *pnSrcXOff, int *pnSrcYOff, int *pnSrcXSize, int *pnSrcYSize,
				 int *pnDstXOff, int *pnDstYOff, int *pnDstXSize, int *pnDstYSize)
{
	/* Check that the destination bounding box intersects the source bounding box */
	if (psDP->adfGeoTransform[GEOTRSFRM_TOPLEFT_X] + psDP->nRasterXSize * psDP->adfGeoTransform[GEOTRSFRM_WE_RES] < minX)
	{
		return FALSE;
	}
	if (psDP->adfGeoTransform[GEOTRSFRM_TOPLEFT_X] > maxX)
	{
		return FALSE;
	}
	if (psDP->adfGeoTransform[GEOTRSFRM_TOPLEFT_Y] + psDP->nRasterYSize * psDP->adfGeoTransform[GEOTRSFRM_NS_RES] > maxY)
	{
		return FALSE;
	}
	if (psDP->adfGeoTransform[GEOTRSFRM_TOPLEFT_Y] < minY)
	{
		return FALSE;
	}

	*pnSrcXSize = psDP->nRasterXSize;
	*pnSrcYSize = psDP->nRasterYSize;
	if (psDP->adfGeoTransform[GEOTRSFRM_TOPLEFT_X] < minX)
	{
		*pnSrcXOff = (int)((minX - psDP->adfGeoTransform[GEOTRSFRM_TOPLEFT_X]) / psDP->adfGeoTransform[GEOTRSFRM_WE_RES] + 0.5);
		*pnDstXOff = 0;
	}
	else
	{
		*pnSrcXOff = 0;
		*pnDstXOff = (int)(0.5 + (psDP->adfGeoTransform[GEOTRSFRM_TOPLEFT_X] - minX) / we_res);
	}
	if (maxY < psDP->adfGeoTransform[GEOTRSFRM_TOPLEFT_Y])
	{
		*pnSrcYOff = (int)((psDP->adfGeoTransform[GEOTRSFRM_TOPLEFT_Y] - maxY) / -psDP->adfGeoTransform[GEOTRSFRM_NS_RES] + 0.5);
		*pnDstYOff = 0;
	}
	else
	{
		*pnSrcYOff = 0;
		*pnDstYOff = (int)(0.5 + (maxY - psDP->adfGeoTransform[GEOTRSFRM_TOPLEFT_Y]) / -ns_res);
	}
	*pnDstXSize = (int)(0.5 + psDP->nRasterXSize * psDP->adfGeoTransform[GEOTRSFRM_WE_RES] / we_res);
	*pnDstYSize = (int)(0.5 + psDP->nRasterYSize * psDP->adfGeoTransform[GEOTRSFRM_NS_RES] / ns_res);

	return TRUE;
}

CPLErr BuildVRT(const char *pszOutputFilename,
					std::vector<std::string> InputFilenames,
					double minX, double minY, double maxX, double maxY,
					GDALProgressFunc pfnProgress, void *pProgressData)
{
	char *projectionRef = NULL;
	int nBands = 0;
	BandProperty *bandProperties = NULL;
	int i, j;
	int rasterXSize;
	int rasterYSize;
	int nCount = 0;
	int bFirst = TRUE;
	VRTDatasetH hVRTDS = NULL;
	CPLErr eErr = CE_None;

	if (pfnProgress == NULL) { pfnProgress = GDALDummyProgress; }

	double we_res = 0;
	double ns_res = 0;

	int nInputFiles = InputFilenames.size();

	DatasetProperty *psDatasetProperties = (DatasetProperty *)CPLCalloc(InputFilenames.size(), sizeof(DatasetProperty));

	int bAllowSrcNoData = TRUE;
	double *padfSrcNoData = NULL;

	int bAllowVRTNoData = TRUE;
	double *padfVRTNoData = NULL;

	for (std::vector<std::string>::size_type i = 0; i < InputFilenames.size(); i++)
	{
		const char *dsFileName = InputFilenames[i].c_str();
		// const char *dsFileName = InputFile.c_str();

		if (!pfnProgress(1.0 * (i + 1) / InputFilenames.size(), NULL, pProgressData))
		{
			eErr = CE_Failure;
			goto end;
		}

		// GDALDatasetH hDS = GDALOpen(InputFilenames[i].c_str(), GA_ReadOnly);
		GDALDatasetH hDS = GDALOpen(dsFileName, GA_ReadOnly);
		psDatasetProperties[i].isFileOK = FALSE;

		if (hDS)
		{
			const char *proj = GDALGetProjectionRef(hDS);
			GDALGetGeoTransform(hDS, psDatasetProperties[i].adfGeoTransform);

			psDatasetProperties[i].nRasterXSize = GDALGetRasterXSize(hDS);
			psDatasetProperties[i].nRasterYSize = GDALGetRasterYSize(hDS);
			double product_minX = psDatasetProperties[i].adfGeoTransform[GEOTRSFRM_TOPLEFT_X];
			double product_maxY = psDatasetProperties[i].adfGeoTransform[GEOTRSFRM_TOPLEFT_Y];
			double product_maxX = product_minX + GDALGetRasterXSize(hDS) * psDatasetProperties[i].adfGeoTransform[GEOTRSFRM_WE_RES];
			double product_minY = product_maxY + GDALGetRasterYSize(hDS) * psDatasetProperties[i].adfGeoTransform[GEOTRSFRM_NS_RES];

			GDALGetBlockSize(GDALGetRasterBand(hDS, 1), &psDatasetProperties[i].nBlockXSize, &psDatasetProperties[i].nBlockYSize);

			int _nBands = GDALGetRasterCount(hDS);

			psDatasetProperties[i].firstBandType = GDALGetRasterDataType(GDALGetRasterBand(hDS, 1));

			psDatasetProperties[i].padfNoDataValues = (double *)CPLMalloc(sizeof(double) * _nBands);
			psDatasetProperties[i].panHasNoData = (int *)CPLMalloc(sizeof(int) * _nBands);
			for (j = 0; j < _nBands; j++)
			{
				psDatasetProperties[i].padfNoDataValues[j] = GDALGetRasterNoDataValue(GDALGetRasterBand(hDS, j + 1), &psDatasetProperties[i].panHasNoData[j]);
			}

			if (bFirst)
			{
				if (proj)
					projectionRef = CPLStrdup(proj);
				minX = product_minX;
				minY = product_minY;
				maxX = product_maxX;
				maxY = product_maxY;
				nBands = _nBands;

				bandProperties = (BandProperty *)CPLMalloc(nBands * sizeof(BandProperty));
				for (j = 0; j < nBands; j++)
				{
					GDALRasterBandH hRasterBand = GDALGetRasterBand(hDS, j + 1);
					// bandProperties[j].colorInterpretation = GDALGetRasterColorInterpretation(hRasterBand);
					bandProperties[j].colorInterpretation = GCI_GrayIndex;
					bandProperties[j].dataType = GDALGetRasterDataType(hRasterBand);
					bandProperties[j].colorTable = 0;
					bandProperties[j].noDataValue = GDALGetRasterNoDataValue(hRasterBand, &bandProperties[j].bHasNoData);
				}
			}
			else
			{
				for (j = 0; j < nBands; j++)
				{
					GDALRasterBandH hRasterBand = GDALGetRasterBand(hDS, j + 1);
					(void)hRasterBand;
				}
				if (j != nBands)
					continue;

				if (product_minX < minX)
					minX = product_minX;
				if (product_minY < minY)
					minY = product_minY;
				if (product_maxX > maxX)
					maxX = product_maxX;
				if (product_maxY > maxY)
					maxY = product_maxY;
			}

			we_res += psDatasetProperties[i].adfGeoTransform[GEOTRSFRM_WE_RES];
			ns_res += psDatasetProperties[i].adfGeoTransform[GEOTRSFRM_NS_RES];

			psDatasetProperties[i].isFileOK = 1;
			nCount++;
			bFirst = FALSE;
			GDALClose(hDS);
		}
		else
		{
			CPLError(CE_Warning, CPLE_AppDefined, "Warning : can't open %s. Skipping it", dsFileName);
		}
	}

	if (nCount == 0)
		goto end;

	we_res /= nCount;
	ns_res /= nCount;

	rasterXSize = (int)(0.5 + (maxX - minX) / we_res);
	rasterYSize = (int)(0.5 + (maxY - minY) / -ns_res);

	if (rasterXSize == 0 || rasterYSize == 0)
	{
		CPLError(CE_Failure, CPLE_AppDefined, "Computed VRT dimension is invalid. You've probably specified unappropriate resolution.");
		goto end;
	}

	hVRTDS = VRTCreate(rasterXSize, rasterYSize);
	GDALSetDescription(hVRTDS, pszOutputFilename);

	if (projectionRef)
	{
		GDALSetProjection(hVRTDS, projectionRef);
	}

	double adfGeoTransform[6];
	adfGeoTransform[GEOTRSFRM_TOPLEFT_X] = minX;
	adfGeoTransform[GEOTRSFRM_WE_RES] = we_res;
	adfGeoTransform[GEOTRSFRM_ROTATION_PARAM1] = 0;
	adfGeoTransform[GEOTRSFRM_TOPLEFT_Y] = maxY;
	adfGeoTransform[GEOTRSFRM_ROTATION_PARAM2] = 0;
	adfGeoTransform[GEOTRSFRM_NS_RES] = ns_res;
	GDALSetGeoTransform(hVRTDS, adfGeoTransform);

	for (j = 0; j < nBands; j++)
	{
		GDALRasterBandH hBand;
		GDALAddBand(hVRTDS, bandProperties[j].dataType, NULL);
		hBand = GDALGetRasterBand(hVRTDS, j + 1);
		GDALSetRasterColorInterpretation(hBand, bandProperties[j].colorInterpretation);
		if (bandProperties[j].colorInterpretation == GCI_PaletteIndex)
		{
			GDALSetRasterColorTable(hBand, bandProperties[j].colorTable);
		}
		if (bAllowVRTNoData && bandProperties[j].bHasNoData)
		{
			GDALSetRasterNoDataValue(hBand, bandProperties[j].noDataValue);
		}

		for (i = 0; i < nInputFiles; i++)
		{
			if (psDatasetProperties[i].isFileOK == 0)
				continue;

			int nSrcXOff, nSrcYOff, nSrcXSize, nSrcYSize, nDstXOff, nDstYOff, nDstXSize, nDstYSize;
			if (!GetSrcDstWin(&psDatasetProperties[i],
							  we_res, ns_res, minX, minY, maxX, maxY,
							  &nSrcXOff, &nSrcYOff, &nSrcXSize, &nSrcYSize,
							  &nDstXOff, &nDstYOff, &nDstXSize, &nDstYSize))
				continue;

			const char *dsFileName = InputFilenames[i].c_str();

			GDALProxyPoolDatasetH hProxyDS = GDALProxyPoolDatasetCreate(dsFileName,
																		psDatasetProperties[i].nRasterXSize,
																		psDatasetProperties[i].nRasterYSize,
																		GA_ReadOnly, TRUE, projectionRef,
																		psDatasetProperties[i].adfGeoTransform);

			for (j = 0; j < nBands; j++)
			{
				GDALProxyPoolDatasetAddSrcBandDescription(hProxyDS,
														  bandProperties[j].dataType,
														  psDatasetProperties[i].nBlockXSize,
														  psDatasetProperties[i].nBlockYSize);
			}

			for (j = 0; j < nBands; j++)
			{
				VRTSourcedRasterBandH hVRTBand = (VRTSourcedRasterBandH)GDALGetRasterBand(hVRTDS, j + 1);

				/* Place the raster band at the right position in the VRT */
				if (bAllowSrcNoData && psDatasetProperties[i].panHasNoData[j])
				{
					VRTAddComplexSource(hVRTBand, GDALGetRasterBand((GDALDatasetH)hProxyDS, j + 1),
										nSrcXOff, nSrcYOff, nSrcXSize, nSrcYSize,
										nDstXOff, nDstYOff, nDstXSize, nDstYSize,
										0, 1, psDatasetProperties[i].padfNoDataValues[j]);
				}
				else
				{
					VRTAddSimpleSource(hVRTBand, GDALGetRasterBand((GDALDatasetH)hProxyDS, j + 1),
									   nSrcXOff, nSrcYOff, nSrcXSize, nSrcYSize,
									   nDstXOff, nDstYOff, nDstXSize, nDstYSize,
									   "near", VRT_NODATA_UNSET);
				}
			}

			GDALDereferenceDataset(hProxyDS);
		}
	}
	GDALClose(hVRTDS);

end:
	for (i = 0; i < nInputFiles; i++)
	{
		CPLFree(psDatasetProperties[i].padfNoDataValues);
		CPLFree(psDatasetProperties[i].panHasNoData);
	}
	CPLFree(psDatasetProperties);
	if (bandProperties != NULL)
	{
		for (j = 0; j < nBands; j++)
		{
			GDALDestroyColorTable(bandProperties[j].colorTable);
		}
	}
	CPLFree(bandProperties);
	CPLFree(projectionRef);
	CPLFree(padfSrcNoData);
	CPLFree(padfVRTNoData);

	return eErr;
}

std::string zero_pad(int i, int len) {
	std::string s = std::to_string(i);
	int s_len = static_cast<int>(s.size());
	int p = len - std::min(len, s_len);
	return std::string(p, '0').append(s);
}

std::string LongName(int i, bool CAPS)
{
	std::string Name;

	std::string H; // Hemisphere label, can be E or W}
	if (0<= i && i <= 180) { H = (CAPS) ? "E" : "e"; }
	else { 
		H = (CAPS) ? "W" : "w";
		if (i > 180) { i = 360 - i; }
		else { i = -i; }
	}

	Name = H.append(zero_pad(i, 3));

	return Name;
}

std::string LatName(int i, bool CAPS)
{

	std::string Name;

	std::string H; // Hemisphere label, can be N or S
	if (i >= 0) { H = (CAPS) ? "N" : "n"; }
	else { H = (CAPS) ? "S" : "s"; }

	Name = H.append(zero_pad(abs(i), 2));
	
	return Name;
}
