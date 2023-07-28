.. _settings_and_parameters:

Settings and parameters
=======================

Kestrel input files are divided into different blocks, which may be specified in
any order. We document the options for each below.

The settings are given in the form :code:`keyword = value`.  Here we provide example values for each keyword.

Some settings are **required** and the simulation will not commence if values are not given.

Other settings are **optional** and default values will be used if they are not given in the input file.  In this case the default value is used as the example.

There are also some settings that are required in only some circumstances (conditionally required).

.. _set_domain:

Domain
------

The *Domain* block specifies the location, size, resolution and boundary conditions that are applied at domain edges.  The domain block is identified using the block keyword :code:`Domain:`.

The domain is an interval in one-dimensional simulations and a rectangular region for two-dimensional simulations.  In one-dimension, the interval is taken to be along the :math:`x` axis.  In two-dimensions, the coordinates at :math:`x` and :math:`y`.  For georeferenced simulations on topography, the domain is rectangular in the WGS84 UTM zone corresponding to the centre of the domain, with :math:`x` corresponding to the Easting coordinate and :math:`y` corresponding to the Northing coordinate.

The domain is divided into equal sized tiles, with the tiles divided into equal sized cells; this allows for efficient computation for flows following irregular paths along topography.

The following **required** settings are specified in the *Domain* block:

    :code:`Lat = 13.248745`

        The latitude of the center of the domain in decimal degrees in WGS84 coordinates, for georeferenced simulations.
        For simulations using an analytical surface, this should be set to the required central point for y-coordinate.

    :code:`Lon = 123.686939`

        The longitude of the center of the domain in decimal degrees in WGS84 coordinates, for georeferenced simulations.
        For simulations using an analytical surface, this should be set to the required central point for x-coordinate.

    :code:`nXtiles = 40`
    
        Number of tiles in the x-direction.

    :code:`nYtiles = 10`

        Number of tiles in the y-direction.

    :code:`nXpertile = 25`

        Number of cells in the x-direction in each tile.

    :code:`nYpertile = 25`

        Number of cells in the y-direction in each tile.

    :code:`Xtilesize = 50`

        Length of a tile in the x-direction in metres.  (Currently this also sets the length of a tile in the y-direction.)

    .. note::

        The numerical resolution of the simulation is determined by the cell size, given by :code:`Xtilesize` / :code:`nXpertile`.

The following **optional** settings can be specified (here given with their default value if left unset):

    :code:`Boundary conditions = halt`

        Sets the conditions that apply at the edge of the domain.
        Options are: 
    
            :code:`halt` (default) -- simulation terminates if flow reaches domain boundary;

            :code:`sponge` -- flow quantities are forced to decay near to the domain boundary;

            :code:`periodic` -- flow quantities at a boundary are transferred to opposite boundary;

            :code:`dirichlet` -- flow quantities at a boundary are specified; values for :math:`H`, :math:`\bar{u}`, :math:`\bar{v}` and :math:`\bar{\psi}` are needed, and are specified using the domain block variables :code:`Boundary Hn`, :code:`Boundary U`, :code:`Boundary V`, and :code:`Boundary psi`, respectively.

The following settings are **conditionally required** when the optional :code:`Boundary conditions = dirichlet`:

    :code:`Boundary Hn`

        Value of flow depth, :math:`H`, to impose at a domain edge if :code:`Boundary conditions = dirichlet`.
        Unused by default.

    :code:`Boundary U`

        Value of flow velocity component in the x-direction, :math:`\bar{u}`, to impose at a domain edge if :code:`Boundary conditions = dirichlet`.
        Unused by default.

    :code:`Boundary V`

        Value of flow velocity component in the y-direction, :math:`\bar{v}`, to impose at a domain edge if :code:`Boundary conditions = dirichlet`.
        Unused by default.

    :code:`Boundary psi`

        Value of flow solids concentration, :math:`\bar{\psi}`, to impose at a domain edge if :code:`Boundary conditions = dirichlet`.
        Unused by default.


Output
------

TODO

Parameters
----------

TODO

Solver
------

TODO

Source
------

A *Source* block specifies conditions for a release of material onto the domain through a time series (referred to as a *flux source*).  A source block is identified using the block keyword :code:`Source:`.

Multiple flux sources can be added through additional Source blocks.

The flux source is modelled as a circular area through which material is added to the domain at a specified volumetric flux and with a specified solids fraction. The flux source requires a location, size and time series for the volumetric flux and solids fraction.

The location of the source can be specified by giving *either*
    
    - the latitude (:code:`sourceLat`) and longitude (:code:`sourceLon`) of the centre of the source;

*or* 

    - the offset of the source centre from the centre of the domain (:code:`sourceX`, :code:`sourceY`), in metres.

.. note::

    If using an artificial analytical topographic surface, the location must be set using :code:`sourceX`, :code:`sourceY`.

These **required** location specifies give:

    :code:`sourceLat = 13.248745`

        The latitude of the center of the flux source in decimal degrees in WGS84 coordinates.

    :code:`sourceLon = 123.686939`

        The longitude of the center of the flux source in decimal degrees in WGS84 coordinates.

    :code:`sourceX = 100`

        The offset of the center of the flux source along the :math:`x` axis in metres from the centre of the domain.

    :code:`sourceY = -50`

        The offset of the center of the flux source along the :math:`y` axis in metres from the centre of the domain.

The following are the additional **required** settings for a source block:

    :code:`sourceRadius = 5`

        The radius of the circular flux source, in metres.
        .. note::

            The radius should be large enough to ensure that the source can be represented on the numerical grid.

    :code:`sourceTime = (  0, 360, 720)`
    
        A list of times for which the volumetric flux and solids fraction are given. 
        This takes the form :code:`sourceTime = (t0, t1, t2, ..., tN)` with ascending times and can contain as many increments as needed.

    :code:`sourceFlux = (5.0, 7.0, 0.0)` 
    
        A list of the volume flux (m:sup:`3`/s) at the times given in :code:`sourceTime`, and takes the form :code:`sourceFlux = (Q0, Q1, Q2, ..., QN)`.

    :code:`sourceConc = (0.0, 0.0, 0.0)`
    
        A list of the solids concentration at the times given in :code:`sourceTime`, and takes the form :code:`sourceConc = (psi0, psi1, psi2, ..., psiN)`.

    .. note::
    
        Each of :code:`sourceTime`, :code:`sourceFlux` and :code:`sourceConc` must contain the same number of points.
        
        For times t<t0 and t>tN, Q=0, psi=0.

        Between the given time increments, the flux and concentration are linearly interpolated.


TODO

Topog
-----

TODO
