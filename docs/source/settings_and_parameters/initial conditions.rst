Initial and source conditions
-----------------------------

Cap
^^^

.. warning::
    TODO

Cube
^^^^

.. warning::
    TODO

Source
^^^^^^

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

