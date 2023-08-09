Domain
------

The *Domain* block specifies the location, size, resolution and boundary
conditions that are applied at domain edges.  The domain block is identified
using the block keyword :code:`Domain:`.

The domain is an interval in one-dimensional simulations and a rectangular
region for two-dimensional simulations.  In one dimension, the interval is taken
to be along the :math:`x` axis.  In two dimensions, the coordinates are :math:`x`
and :math:`y`.  For georeferenced simulations on topography, the domain is
rectangular in the WGS84 UTM zone corresponding to the centre of the domain,
with :math:`x` corresponding to the Easting coordinate and :math:`y`
corresponding to the Northing coordinate.

The domain is divided into equal sized tiles, with the tiles divided into equal
sized cells; this allows for efficient computation for flows following irregular
paths along topography.

The following **required** settings are specified in the *Domain* block:

    :code:`Lat = 13.248745`

        The latitude of the center of the domain in decimal degrees in WGS84
        coordinates, for georeferenced simulations.  For simulations using an
        analytical surface, this should be set to the required central point for
        :math:`y`-coordinate.

    :code:`Lon = 123.686939`

        The longitude of the center of the domain in decimal degrees in WGS84
        coordinates, for georeferenced simulations.  For simulations using an
        analytical surface, this should be set to the required central point for
        :math:`x`-coordinate.

    :code:`nXtiles = 40`
    
        Number of tiles in the :math:`x`-direction.

    :code:`nYtiles = 10`

        Number of tiles in the :math:`y`-direction.

    :code:`nXpertile = 25`

        Number of finite volume cells in the :math:`x`-direction in each tile.

    :code:`nYpertile = 25`

        Number of finite volume cells in the :math:`y`-direction in each tile.

    :code:`Xtilesize = 50`

        Length of a tile in the :math:`x`-direction in metres.  (Currently this
        also sets the length of a tile in the :math:`y`-direction.)

    .. note::

        The numerical resolution of the simulation is determined by the cell
        size, given by :code:`Xtilesize` / :code:`nXpertile`.

    .. note::

        To perform one-dimensional simulations, set :code:`nYtiles = 1` and
        :code:`nYpertile = 1`.

The following **optional** settings can be specified (here given with their default value if left unset):

    :code:`Boundary conditions = halt`

        Sets the conditions that apply at the edge of the domain.
        Options are: 
    
            :code:`halt` (default) -- simulation terminates if flow reaches domain boundary;

            :code:`sponge` -- flow quantities are forced to decay near to the domain boundary;

            :code:`periodic` -- flow quantities at a boundary are transferred to opposite boundary;

            .. warning::
                
                Periodic boundary conditions *do not* enforce periodicity in the
                initial bed elevation. This means that applying them for
                non-periodic topographies (such as constant slopes) will lead to
                step changes in :math:`b` across the boundaries.

            :code:`dirichlet` -- flow quantities at a boundary are specified;
            values for :math:`H`, :math:`\bar{u}`, :math:`\bar{v}` and
            :math:`\bar{\psi}` are needed, and are specified using the domain
            block variables :code:`Boundary Hn`, :code:`Boundary U`,
            :code:`Boundary V`, and :code:`Boundary psi`, respectively.

The following settings are **conditionally required** when the optional
:code:`Boundary conditions = dirichlet`:

    :code:`Boundary Hn`

        Value of flow depth, :math:`H`, to impose at a domain edge if
        :code:`Boundary conditions = dirichlet`.  Unused by default.

    :code:`Boundary U`

        Value of flow velocity component in the :math:`x`-direction,
        :math:`\bar{u}`, to impose at a domain edge if :code:`Boundary
        conditions = dirichlet`.  Unused by default.

    :code:`Boundary V`

        Value of flow velocity component in the :math:`y`-direction,
        :math:`\bar{v}`, to impose at a domain edge if :code:`Boundary
        conditions = dirichlet`.  Unused by default.

    :code:`Boundary psi`

        Value of flow solids concentration, :math:`\bar{\psi}`, to impose at a
        domain edge if :code:`Boundary conditions = dirichlet`.  Unused by
        default.

    .. note:: 

        Dirchlet conditions only affect parts of the simulation containing flow.
        Therefore, the typical way to make use of them is to provide some
        flowing material via a `Cap` block (see below) that touches the edge of
        the domain.
