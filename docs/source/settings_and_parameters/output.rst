Output
------

The *Output* block sets up output from Kestrel. It is identified using the block
keyword :code:`Output:`.

The only **required** setting in the *Output* block is:

    :code:`N out = 10`

        The number of output files to be produced. These record simulation data
        at equally spaced time intervals over the duration of the simulation.

        .. note::

            The initial condition is also recorded and does not count towards
            this total.

The **optional** settings of the output block are:

    :code:`base path = ./`

        A path to a base directory to hold the output directory.  Default is
        current working directory.  This is created if it does not exist,
        providing that permissions allow.

    :code:`directory = results/`

        A directory to store the results.  This is created if it does not exist
        and permissions allow.

    :code:`format = txt`

        Format of the output files.  Options are:

            :code:`txt` (default) -- column-headed, comma-delimited text files.

            :code:`nc` or :code:`netcdf` -- NetCDF files.  Requires compilation with NetCDF4 (see :ref:`installation`)

            :code:`kml` -- KML files.  Requires simulation on georeferenced topography.
        
        .. note::

            Multiple output formats can be specified as a comma-separated list (e.g. :code:`format = txt, nc`).
        
        .. note::

            KML output feature is currently limited.

    :code:`info filename = RunInfo.txt`

        Name of a text file to contain data on the simulation.

The following **conditionally optional** setting can be given if :code:`format` value includes :code:`txt`, :code:`nc` or :code:`netcdf`:

    :code:`maximums filename = Maximums`

        The name of a file to contain aggregated data over the duration of the simulation.
        For a text file, the extension .txt if added.
        For a NetCDF file, the extension .nc is added.

The following **conditionally optional** setting can be given if the :code:`format` value includes :code:`txt`:

    :code:`compression = off`

        Compress text files using tar.

The following **conditionally optional** setting can be given if :code:`format` value includes :code:`kml`:

    :code:`kml height = 0.1`

        A threshold depth, in metres, for writing flow depth and flow speed outputs to KML files.  Cells with flow depths below :code:`kml height` are ignored when writing to KML.

    :code:`inundation time filename = InundationTime`

        The name of a KML file to store the first time of inundation of points in the domain.

    :code:`max height filename = MaxHeights`

        The name of a KML file to store the maximum flow depth, and the time of this maximum, for points in the domain.

    :code:`max speed filename = MaxSpeeds`

        The name of a KML file to store the maximum flow speed, and the time of this maximum, for points in the domain.

