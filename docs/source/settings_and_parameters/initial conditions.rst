Initial and source conditions
-----------------------------

Flowing material may be input into simulations via ``Cap``, ``Cube`` and
``Source`` blocks. Any number of these may be added to the input file.

.. note::
   The ``Initial condition`` declaration, which is part of the ``Solver`` block
   (described below) may also be used to specify arbitrary initial conditions by
   providing valid solution data file(s).

Cap
^^^

A *Cap* block places a volume of flowing material with a circular base into the
simulation at the initial time. The following declarations are **required**:

The location of the cap can be specified by giving *either*
    
    - the latitude (:code:`capLat`) and longitude (:code:`capLon`) of its centre

*or* 

    - the offset of the cap centre from the centre of the domain (:code:`capX`,
      :code:`capY`), in metres.

These **required** specifiers give:

    :code:`capLat = 13.248745`

        The latitude of the center of the volume in decimal degrees in WGS84
        coordinates.

    :code:`capLat = 123.686939`

        The longitude of the center of the volume in decimal degrees in WGS84
        coordinates.

    :code:`capX = 100`

        The offset of the center of the volume along the :math:`x` axis in
        metres from the centre of the domain.

    :code:`capY = -50`

        The offset of the center of the volume along the :math:`y` axis in
        metres from the centre of the domain.

To set the volume of the release, exactly two out of three of the following
settings are **required**:

    :code:`capRadius`

        Sets the radius of the circular base in metres.

        .. note::

            The radius should be large enough to ensure that the volume can be
            represented on the numerical grid.

    :code:`capHeight`

        Sets the height of the volume in metres. N.b. the precise meaning of
        this depends on the value of :code:`capShape` (see below).

    :code:`capVolume`

        Sets the volume in metres.

From the two specified values, Kestrel determines the third.

The following **optional** parameters set additional attributes for the volume
release:

    :code:`capU = 0`

        Sets the initial velocity :math:`\bar{u}` of the release in the
        :math:`x`-direction (m/s).

    :code:`capV = 0`

        Sets the initial velocity :math:`\bar{v}` of the release in the
        :math:`y`-direction (m/s).

    :code:`capConc = 0`

        Sets the volumetric fraction :math:`\bar{\psi}` of solids present in the
        release.

    :code:`capShape = para`

        This selects the geometry of the volume. There are three options to
        choose from

            :code:`flat`

                Selects a volume with a flat top with height ``capHeight``, i.e.
                a cylinder.

            :code:`level`

                Selects a volume whose top surface lies at constant vertical
                elevation. The value of this elevation is set to ``capHeight``.

                .. note::

                    In this case only, both :code:`capRadius` and
                    :code:`capHeight` must be set.

            :code:`para` (default)

                Selects a parabolic volume with formula

                    :math:`H(\mathbf{x}) = C_H (1 - |\mathbf{x} -
                    \mathbf{x}_c|^2/C_R)`,

                where :math:`C_H` is ``capHeight``, :math:`C_R` is ``capRadius``
                and :math:`\mathbf{x}_c` is the centre point of the volume's
                base at ``capX``, ``capY``.

Cube
^^^^

A *Cube* block places a volume of flowing material with a rectangular base into
the simulation at the initial time. The following declarations are **required**:

The location of the cube can be specified by giving *either*
    
    - the latitude (:code:`cubeLat`) and longitude (:code:`cubeLon`) of its centre

*or* 

    - the offset of the cube centre from the centre of the domain (:code:`cubeX`,
      :code:`cubeY`), in metres.

These **required** specifiers give:

    :code:`cubeLat = 13.248745`

        The latitude of the center of the volume in decimal degrees in WGS84
        coordinates.

    :code:`cubeLat = 123.686939`

        The longitude of the center of the volume in decimal degrees in WGS84
        coordinates.

    :code:`cubeX = 100`

        The offset of the center of the volume along the :math:`x` axis in
        metres from the centre of the domain.

    :code:`cubeY = -50`

        The offset of the center of the volume along the :math:`y` axis in
        metres from the centre of the domain.

To set the volume of the release, the following settings are **required**:

    :code:`cubeHeight`

        Sets the height of the volume in metres. N.b. the precise meaning of
        this depends on the value of :code:`cubeShape` (see below).

    :code:`cubeLength`

        Sets the extent of the volume in the :math:`x`-direction.

    :code:`cubeWidth`

        Sets the extent of the volume in the :math:`y`-direction.

    .. note::

        Both ``cubeLength`` and ``cubeWidth`` should be large enough to ensure
        that the volume can be represented on the numerical grid.

The following **optional** parameters set additional attributes for the volume
release:

    :code:`cubeU = 0`

        Sets the initial velocity :math:`\bar{u}` of the release in the
        :math:`x`-direction (m/s).

    :code:`cubeV = 0`

        Sets the initial velocity :math:`\bar{v}` of the release in the
        :math:`y`-direction (m/s).

    :code:`cubeConc = 0`

        Sets the volumetric fraction :math:`\bar{\psi}` of solids present in the
        release.

    :code:`cubeShape = flat`

        This selects the geometry of the volume. There are three options to
        choose from

            :code:`flat (default)`

                Selects a volume with a flat top with height ``cubeHeight``,
                i.e. a cuboid.

            :code:`level`

                Selects a volume whose top surface lies at constant vertical
                elevation. The value of this elevation is set to ``capHeight``.

Source
^^^^^^

A *Source* block specifies conditions for a release of material onto the domain
through a time series (referred to as a *flux source*).  A source block is
identified using the block keyword :code:`Source:`.

Multiple flux sources can be added through additional Source blocks.

The flux source is modelled as a circular area through which material is added
to the domain at a specified volumetric flux and with a specified solids
fraction. The flux source requires a location, size and time series for the
volumetric flux and solids fraction.

The location of the source can be specified by giving *either*
    
    - the latitude (:code:`sourceLat`) and longitude (:code:`sourceLon`) of the centre of the source;

*or* 

    - the offset of the source centre from the centre of the domain (:code:`sourceX`, :code:`sourceY`), in metres.

.. note::

    If using an artificial analytical topographic surface, the location must be set using :code:`sourceX`, :code:`sourceY`.

These **required** specifiers give:

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
    
        A list of the volume flux (m\ :sup:`3`/s) at the times given in :code:`sourceTime`, and takes the form :code:`sourceFlux = (Q0, Q1, Q2, ..., QN)`.

    :code:`sourceConc = (0.0, 0.0, 0.0)`
    
        A list of the solids concentration at the times given in :code:`sourceTime`, and takes the form :code:`sourceConc = (psi0, psi1, psi2, ..., psiN)`.

    .. note::
    
        Each of :code:`sourceTime`, :code:`sourceFlux` and :code:`sourceConc` must contain the same number of points.
        
        For times t<t0 and t>tN, Q=0, psi=0.

        Between the given time increments, the flux and concentration are linearly interpolated.

