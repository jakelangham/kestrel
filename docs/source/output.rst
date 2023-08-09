.. _output:

Output
======

Kestrel simulation results are contained in a number of output files.
The main numerical results of the simulation are cell-averaged values 
of quantities at locations within the flow domain.  There are two types
of numerical results files:

    * *snapshot* files that contain data produced by Kestrel at a single
      point in time. The files have sequential numeric names, e.g.
      *000001.nc*, *000002.nc*, .... The snapshot file
      *000000.nc* contains any initial conditions of the flow.

    * *aggregated* files that contain data that summarize the flow up to
      the point at which they are written.  The aggregated results file
      update during the simulation.  Aggregated results are stored in a single
      file with the filename set using the optional input
      :code:`maximums filename` with default value :code:`Maximums`.

The numerical results are only stored in *active* regions of the domain,
and therefore, the data volume in later snapshot files is larger than in
earlier files, and aggregated results files increase in size as simulations
progress.

The results files can be provided in two main formats: as plain text in 
column-headed comma-separated values (with a *.txt* extension),
or as NetCDF files (with a *.nc* extension).  Here we first describe
the output files in these two formats.  Additionally, Kestrel produces:

    * a *RunInfo.txt* file that summarizes the settings used in the model;
    * a *Volume.txt* file that summarizes the time series of volumes of material;
    * files related to the topographic data are created.

.. _output_netcdf:

NetCDF files
------------

NetCDF files are typically smaller in size and contain data at higher
numerical precision than Kestrel's text file outputs, and are
georeferenced when simulations use topographic data.  Therefore, NetCDF are our
preferred format, but require a compatible NetCDF library for installation
(see :ref:`installation`) and post-processing tools that are able to read NetCDF
files.  Many GIS software packages are capable of reading NetCDF, and packages are
available in Julia, MatLab, Python, and other languages.

Computed solution fields are stored as Float64 variables.

The NetCDF files contain global attributes that contain some of the important settings
of the model.  These include the attribute :code:`time` that gives the simulated time, 
in seconds, of the snapshot data or latest write for aggregated results.

For georeferenced simulations, the global attribute :code:`crs_epsg` gives the EPSG code
for the UTM coordinate reference system.

There are four NetCDF dimensions defined:
    * :code:`x` - the Easting dimension, for cell centres;
    * :code:`y` - the Northing dimension for cell centres;
    * :code:`x_vertex` - the Easting dimension for cell vertices (only for *snapshot* files);
    * :code:`y_vertex` - the Northing dimension for cell vertices (only for *snapshot* files).

.. note::

    For one-dimensional simulations, :code:`y` has size 1, and :code:`y_vertex` has size 2.

.. note::
    
    The vertex data is required when restarting simulations but should not be used for post-processing.
    Therefore, the vertex data is not georeferenced.

For georeferenced simulations, the NetCDF files contain a variable

    :code:`crs`
        A NetCDF variable to store the following attributes defining the coordinate reference system:
            - :code:`grid_mapping_name` taking the value :code:`transverse_mercator`.
            - :code:`false_easting` taking the value :code:`500000`.
            - :code:`false_northing` taking the value :code:`0` for domain centres in the Northern hemisphere and :code:`10000000` for domain centres in the Southern hemisphere.
            - :code:`latitude_of_projection_origin` taking the value :code:`0`.
            - :code:`longitude_of_central_meridian` with value dependent on the central longitude of the domain.
            - :code:`scale_factor_at_central_meridian` taking the value :code:`0.9996`.

.. warning::
    The aggregated results NetCDF file is overwritten during the simulation. This will fail is the NetCDF file is open
    during an attempt to write.  This will result in an error and termination of the simulation.

    It is recommended to copy the file before opening to perform analysis and processing while a Kestrel simulation is in progress.

.. _output_txt:

Text files
----------

The text format output files produced by Kestrel contain column-headed, comma-separated numerical data.
The text files contain the same data fields as NetCDF files, but at a reduced numerical precision.
Additionally, to aid post-processing of these files, when simulations are georeferenced there are
additional columns providing the cell-centre latitude and longitude in the WGS84 coordinate reference system.
Furthermore, to reduce file size, fields that are appropriate only for two-dimensional simulations are not
recorded when one-dimensional simulations are performed.

Computed solution fields are stored in scientific form (:code:`ES18.10E3`).

When restarting simulations, topographic data is required at cell vertices as well as at cell centres.  These data are
stored in separate files with extension *.txt_topo*.

The first column of text data files is named :code:`tile` and refers an internal integer identifier of domain tiles
within Kestrel, and is required for restarting simulations from an existing result file.

.. _output_fields:

Output fields
-------------

The *snapshot* and *aggregated* result files contain the following fields for location:

    :code:`x`
        The Easting coordinate of cell centres [m].
        
        When georeferenced, the coordinate is in the UTM coordinate reference system
        given by :code:`crs_epsg`.
    
    :code:`y`
        The Northing coordinate of cell centres [m].
        
        When georeferenced the coordinate is in the UTM coordinate reference system
        given by :code:`crs_epsg`.
        
        For one-dimensional simulations :code:`y = 0` in NetCDF files and is absent
        from text files.

    :code:`x_vertex`
        The Easting coordinate of cell vertices [m].

        These values are not georeferenced, so are relative to the domain centre.
        
        Only defined in NetCDF *snapshot* files.

    :code:`y_vertex`
        The Northing coordinate of cell vertices [m].
        
        These values are not georeferenced, so are relative to the domain centre.
        
        Only defined in NetCDF *snapshot* files.
        
        For one-dimensional simulations :code:`y_vertex = -0.5, 0.5`.

The *snapshot* files contain the following solution fields:
    
    :code:`flow_depth`
        The flow depth, :math:`H` [m].

    :code:`flow_speed`
        The flow slope-aligned speed, :math:`\left|\mathbf{u}\right|` [m/s]

    :code:`x_velocity`
        The :math:`x`-component of the flow velocity, :math:`\bar{u}` [m/s].
    
    :code:`y_velocity`
        The :math:`y`-component of the flow velocity, :math:`\bar{v}` [m/s]
    
    :code:`density`
        The density of the mixture, :math:`\bar{\rho}` [kg/m\ :sup:`3`\ ]

    :code:`solids_fraction`
        The volume fraction of solids in the mixture, :math:`\bar{\psi}` [*dimensionless*]
    
    :code:`x_flux`
        The :math:`x`-component of the mass flux per unit area, :math:`\bar{\rho} H\bar{u}` [m\ :sup:`3`\ /s]
    
    :code:`y_flux`
        The :math:`y`-component of the volumetric flux per unit area, :math:`\bar{\rho} H\bar{v}` [m\ :sup:`3`\ /s]
    
    :code:`Hnpsi`
        The volume of solids per unit area, :math:`H\bar{\psi}` [m] 
    
    :code:`base_elevation`
        The initial topographic elevation, :math:`b_{0} = b(x,y,0)` [m]
    
    :code:`elevation_change`
        The change in topographic elevation, :math:`\delta b_{t} = b(x,y,t) - b(x,y,0)` [m]
    
    :code:`x_slope`
        The topographic slope along the :math:`x` coordinate, :math:`\partial b/\partial x` [*dimensionless*]
    
    :code:`y_slope`
        The topographic slope along the :math:`y` coordinate, :math:`\partial b/\partial y` [*dimensionless*]
    
    :code:`B0_vertex`
        The initial topographic elevation at cell vertices, :math:`b_{0} = b(x,y,0)` [m]
        This data is required when restarting simulations, but should not be used for post-processing.
        It is not georeferenced.
        Stored in :code:`.txt_topo` file for text output.
    
    :code:`Bt_vertex`
        The change in topographic elevation at cell vertices, :math:`\delta b_{t} = b(x,y,t) - b(x,y,0)` [m]
        This data is required when restarting simulations, but should not be used for post-processing.
        It is not georeferenced.
        Stored in :code:`.txt_topo` file for text output.

    :code:`w`
        The conserved quantity :math:`w = H/\gamma + b` that is computed in the model.
        This is required for restarting simulations but should not be used for post-processing.
    
The *aggregated* result files contain the following solution fields:

    :code:`max_depth`
        The maximum flow depth that has occurred in each cell of the domain [m].

    :code:`t_max_depth`
        The time at which the maximum flow depth occurred in each cell of the domain [s].

    :code:`max_speed`
        The maximum flow speed that has occurred in each cell of the domain [m/s].

    :code:`t_max_speed`
        The time at which the maximum flow depth occurred in each cell of the domain [s].

    :code:`max_erosion`
        The maximum depth of erosion that has occurred in each cell of the domain [m].

    :code:`t_max_erosion`
        The time at which the maximum depth of erosion occurred in each cell of the domain [s].

    :code:`max_deposit`
        The maximum depth of deposition that has occurred in each cell of the domain [m].

    :code:`t_max_deposit`
        The time at which the maximum depth of deposition occurred in each cell of the domain [s].

    :code:`max_solids_frac`
        The maximum solids volume fraction that has occurred in each cell of the domain [*dimensionless*].

    :code:`t_max_solids_frac`
        The time at which the maximum volume fraction occurred in each cell of the domain [s].

    :code:`inundation_time`
        The time at which flow material first reaches each cell of the domain [s].

.. _output_volume:

Flow Volumes
------------

The file :code:`Volume.txt` contains time series recording the evolution of flow volumes and masses
during the simulation.  These are summary values, integrated over the full simulation domain.

The data are stored as column-headed, comma-separated values with the first row 
recording the time, and subsequent columns recording quantities calculated from computed fields.
The calculated values are written at high precision (sixteen figures) as they can be used to
verify accurate computation of conserved quantities.

The following columns are stored:

    :code:`volume`
        The total volume of material in the flow domain.

        This volume includes material added in initial conditions,
        material added by flux sources, and material entrained into
        flows by erosion, and is given by

        .. math::
            V_{total} = \int_{A} H\gamma\ \mathrm{d}A.



    :code:`total_bed_volume`
        The total volume of material derived from the bed.

        This volume is the difference of material eroded from the bed
        from that deposited to the bed, and is given by

        .. math::
            V_{bed} = \int_{A} \left(b(\mathbf{x},t) - b(\mathbf{x},0)\right)\ \mathrm{d}A
        
        where :math:`A` is the area of the flow.

    :code:`total_mass`
        The total mass of material in the flow domain.

        This mass includes material added in initial conditions,
        material added by flux sources, and material entrained into
        flows by erosion, and is given by

        .. math::
            M_{total} = \int_{A} \bar{\rho}H\gamma\ \mathrm{d}A.
        
    :code:`bed_mass`
        The total mass of material derived from the bed.

        This mass is the difference of material eroded from the bed
        from that deposited to the bed, and is given by

        .. math::
            M_{bed} = \rho_{b}V_{bed}

        where :math:`\rho_{b}` is the density of bed material.

    :code:`total_solids_mass`
        The total mass of solids in the flow domain.

        This mass includes solids added in initial conditions,
        solids added by flux sources, and solids entrained into
        flows by erosion, and is given by

        .. math::
            M_{solids} = \int_{A} \rho_{s}H\bar{\psi}\gamma\ \mathrm{d}A.

    :code:`bed solids mass`
        The total mass of solids derived from the bed.

        This mass is the difference of solid mass added by erosion from the bed
        from the solid mass deposited to the bed, and is given by

        .. math::
            M_{bed solids} = \rho_{s}V_{bed}p
        
        where :math:`p` is the bed porosity.

.. warning::
    TODO RunInfo.txt
