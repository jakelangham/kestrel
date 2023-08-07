.. _output:

Output
======

Kestrel simulation results are contained in a number of output files.
The main numerical results of the simulation are cell-averaged values 
of quantities at locations within the flow domain.  There are two types
of numerical results files:

    - *snapshot* files that contain data produced by Kestrel at a single
    point in time. The files have sequential numeric names, e.g.
    :code:`000001.nc`, :code:`000002.nc`, ....  The snapshot file
    :code:`000000.nc` contains any initial conditions of the flow.
    
    - *aggregated* files that contain data that summarize the flow up to
    the point at which they are written.  The aggregated results file
    update during the simulation.

The numerical results are only stored in *active* regions of the domain,
and therefore, the data volume in later snapshot files is larger than in
earlier files, and aggregated results files increase in size as simulations
progress.

The results files can be provided in two main formats: as plain text in 
column-headed comma-separated values (with a :code:`.txt` extension),
or as NetCDF files (with a :code:`.nc` extension).  Here we first describe
the output files in these two formats.  Additionally a RunInfo.txt file is
created that summarizes the settings used in the model, and files related to the 
topographic data are created.  We describe these later.

.. _output_netcdf:

NetCDF files
------------

NetCDF files are typically smaller in size and contain data at higher
numerical precision than Kestrel's text file outputs, and are
georeferenced when simulations use topographic data.  Therefore, NetCDF are our
preferred format, but require a compatible NetCDF library for installation
(see :ref:`installation`_) and post-processing tools that are able to read NetCDF
files.  Many GIS software packages are capable of reading NetCDF, and packages are
available in Julia, MatLab, Python, and other languages.

Aggregated results are stored in the single :code:`Maximums.nc` file for NetCDF outputs.

The NetCDF files contain global attributes that contain some of the important settings
of the model.  These include the attribute :code:`time` that gives the simulated time, 
in seconds, of the snapshot data or latest write for aggregated results.

For georeferenced simulations, the global attribute :code:`crs_epsg` gives the EPSG code
for the UTM coordinate reference system.

There are four NetCDF dimensions defined:
    - :code:`x` - the Easting dimension, for cell centres;
    - :code:`y` - the Northing dimension for cell centres;
    - :code:`x_vertex` - the Easting dimension for cell vertices;
    - :code:`y_vertex` - the Northing dimension for cell vertices.

.. note::
    
    The vertex data is required when restarting simulations, but should not be used for post-processing.
    Therefore, the vertex data is not georeferenced.

The following variables in the NetCDF file (given with their dimensions) define the location of the data:

    :code:`x(x)`
        The Easting coordinate of cell centres, in metres.
        When georeferenced the coordinate is in the UTM coordinate reference system
        given by :code:`crs_epsg`.
    

    :code:`y(y)`
        The Northing coordinate of cell centres, in metres.
        When georeferenced the coordinate is in the UTM coordinate reference system
        given by :code:`crs_epsg`.

    :code:`x_vertex(x_vertex)`
        The Easting coordinate of cell vertices, in metres.
        These values are not georeferenced, so are relative to the domain centre.

    :code:`y_vertex(y_vertex)`
        The Northing coordinate of cell vertices, in metres.
        These values are not georeferenced, so are relative to the domain centre.
    
    :code:`crs()`
        A NetCDF variable to store the following attributes defining the coordinate reference system:
            - :code:`grid_mapping_name` taking the value :code:`transverse_mercator`.
            - :code:`false_easting` taking the value :code:`500000`.
            - :code:`false_northing` taking the value :code:`0` for domain centres in the Northern hemisphere
            and :code:`10000000` for domain centres in the Southern hemisphere.
            - :code:`latitude_of_projection_origin` taking the value :code:`0`.
            - :code:`longitude_of_central_meridian` with value dependent on the central longitude of the domain.
            - :code:`scale_factor_at_central_meridian` taking the value :code:`0.9996`.


The remaining variables (given with their units in brackets) in *snapshot* result files are:
    
    :code:`flow_depth(y, x)`
        The flow depth, :math:`H` [m].

    :code:`x_velocity(y, x)`
        The :math:`x`-component of the flow velocity (m/s), :math:`\bar{u}` [m/s].
    
    :code:`y_velocity(y, x)`
        The :math:`y`-component of the flow velocity (m/s), :math:`\bar{v}` [m/s]
    
    :code:`flow_speed(y, x)`
        The flow speed, :math:`\sqrt{\bar{u}^{2} + \bar{v}^{2}}` [m/s]
    
    :code:`density(y, x)`
        The density of the mixture, :math:`\bar{\rho}` [kg/m:superscript:`3`]

    :code:`solids_fraction(y, x)`
        The volume fraction of solids in the mixture, :math:`\bar{\psi}` [*dimensionless*]
    
    :code:`x_flux(y, x)`
        The :math:`x`-component of the volumetric flux per unit area, :math:`H\bar{u}` [m:superscript:`3`/s]
    
    :code:`y_flux(y, x)`
        The :math:`y`-component of the volumetric flux per unit area, :math:`H\bar{v}` [m:superscript:`3`/s]
    
    :code:`Hnpsi(y, x)`
        The volume of solids per unit area, :math:`H\bar{\psi}` [m] 
    
    :code:`base_elevation(y, x)`
        The initial topographic elevation, :math:`b_{0} = b(x,y,0)` [m]
    
    :code:`elevation_change(y, x)`
        The change in topographic elevation, :math:`\delta b_{t} = b(x,y,t) - b(x,y,0)` [m]
    
    :code:`x_slope(y, x)`
        The topographic slope along the :math:`x` coordinate, :math:`\partial b/\partial x` [*dimensionless*]
    
    :code:`y_slope(y, x)`
        The topographic slope along the :math:`y` coordinate, :math:`\partial b/\partial y` [*dimensionless*]
    
    :code:`B0_vertex(y_vertex, x_vertex)`
        The initial topographic elevation at cell vertices, :math:`b_{0} = b(x,y,0)` [m]
        This data is required when restarting simulations, but should not be used for post-processing.
        It is not georeferenced.
    
    :code:`Bt_vertex(y_vertex, x_vertex)`
        The change in topographic elevation at cell vertices, :math:`\delta b_{t} = b(x,y,t) - b(x,y,0)` [m]
        This data is required when restarting simulations, but should not be used for post-processing.
        It is not georeferenced.

    :code:`w(y, x)`
        The conserved quantity :math:`w = H/\gamma + b` that is computed in the model.
        This is required for restarting simulations but should not be used for post-processing.
    
    

.. warning::
    TODO
