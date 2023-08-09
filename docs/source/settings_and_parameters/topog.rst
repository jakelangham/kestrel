Topog
-----

The *Topog* block specifies the topography to be used for the simulation.
Topography can be given either as a Digital Elevation Model (DEM) or through a
parameterized analytic function.  The *Topog* block is identified using the block
keyword :code:`Topog:`.

There is a single **required** parameter in the *Topog* block:

    :code:`Type`

        The type of topography to use. Options are:

        :code:`DEM` or :code:`Raster`

            Use a georeferenced DEM contained in a geotiff file.
        
        :code:`SRTM`

            Use a georeferenced SRTM file, contained in an SRTM archive.

        :code:`Function`

            Use a parameterized analytic function defining a surface.

The following statment is **conditionally required** when :code:`Type = DEM` or
:code:`Type = Raster`:

    :code:`raster file`

        The name of the georeferenced raster file containing the DEM.
    
The following **conditionally optional** parameters may be used when :code:`Type
= DEM` or :code:`Type = Raster`:

    :code:`dem directory`

        The directory containing :code:`raster file`.  Default is the current
        working directory.

    :code:`embed raster`

        Should the raster DEM contained in :code:`raster file` be embedded in
        SRTM topography?  Options are:

          - :code:`Yes` or :code:`On` -- embed the DEM within SRTM data.
            
            .. warning::
                If :code:`embed raster = on` then flows will transition between
                the DEM topography to SRTM at the edge of the DEM. While this
                can be useful for simulations where the available DEM coverage
                is partial, any offsets between the DEM and SRTM surface maps
                can cause problems.

          - :code:`No` or :code:`Off` -- do no embed the DEM.

            .. warning::
                If :code:`embed raster = off` then error will occur if the flow
                reaches the edge of the DEM.

The following parameter is **conditionally required** when :code:`Type = DEM` or
:code:`Type = Raster` *and* :code:`embed raster = on`. Additionally, it is
**conditionally required** when :code:`Type = SRTM`.

    :code:`srtm directory`

        The directory containing the SRTM data. This should be in the form of
        an SRTM archive, i.e. with a structure N01/N01*.tif, N02/N02*.tif, ...

        .. note::
            SRTM files may be either geotiff or USGS hgt format (included zipped
            hgt files with extension :code:`.SRTMGL1.hgt.zip`).

There is one **conditionally required** setting needed when :code:`Type =
Function` -- :code:`Topog function` giving the choice of function :math:`b_0` to
be used for the initial bed elevation.  This is accompanied by a **conditionally
optional** setting -- :code:`Topog params`, which specifies a list of parameters
to use with the function.  The :code:`Topog params` variable takes the form
:code:`Topog params = (a, b, c, ...)`. Its values convey a particular meaning
for each function choice, as described below for the currently implemented
functions:

    :code:`Topog function = flat`

        A flat surface, :math:`b_0(x,y) = 0`.

        :code:`Topog params` not used.

    :code:`Topog function = xslope`

        A constant slope along the :math:`x`-coordinate, :math:`b_0(x,y) = \alpha x`.

        :code:`Topog params = (alpha)`
            - :code:`alpha` -- The slope of the plane.
    
    :code:`Topog function = yslope`

        A constant slope along the :math:`y`-coordinate, :math:`b_0(x,y) = \beta y`.

        :code:`Topog params = (beta)`

            - :code:`beta` -- The slope of the plane.

    :code:`Topog function = xyslope`

        A general plane inclined along both the :math:`x`- and
        :math:`y`-coordinates, :math:`b_0(x,y) = \alpha x + \beta y`.

        :code:`Topog params = (alpha, beta)`

            - :code:`alpha` -- The slope of the plane along :math:`x`.
            - :code:`beta` -- The slope of the plane along :math:`y`.

    :code:`Topog function = x2slopes`

        A surface with two slopes in the :math:`x`-direction, connected by a circular arc near :math:`x = 0`.

        :code:`Topog params = (alpha, beta, R)`

            - :code:`alpha` -- slope on left-hand-side.
            - :code:`beta` -- slope on right-hand-side.
            - :code:`R` -- radius of connecting cicular arc.

    :code:`Topog function = xBislope`

        A surface with two slopes in the limits :math:`x\to\pm\infty`, connected
        by a smooth transition.  The surface has the form 

        :math:`b_{0}(x,y) = -\tfrac{1}{2}\left(\tan\phi_{1} +
        \tan\phi_{2}\right)x + \tfrac{1}{2}\left(\tan\phi_{1} -
        \tan\phi_{2}\right)\lambda\log\left[\cosh\left(x/\lambda\right)\right].`
        
        :code:`Topog params = (phi1, phi2, lambda)`

            - :code:`phi1` -- the slope angle for :math:`x\to -\infty`, in degrees.  A positive value corresponds to an elevation decreasing from left to right.
            - :code:`phi2` -- the slope angle for :math:`x\to +\infty`, in degrees.  A positive value corresponds to an elevation decreasing from left to right.
            - :code:`lambda` -- the characteristic length scale of the smooth transition region.

    :code:`Topog function = USGS`

        Parameterization of the USGS flume.  This has slope of 31° for
        :math:`x<0`, and slope 2.4° for :math:`x>x_{1}>0` that are connected
        by a smooth :math:`\cosh` curve section.  Note :math:`x_{1}` is
        determined to ensure smooth connection of the slope elements.  The flume
        is confined by walls for :math:`x<8.5` m, that are represented as
        :math:`\tanh` profile humps.  See `Iverson et al. (2010)
        <https://doi.org/10.1029/2009JF001514>`_ for details.

        :code:`Topog params = (wallH, sigma)`

            - :code:`wallH` -- the height of the sidewalls of the flume.
            - :code:`sigma` -- the width of the sidewalls of the flume.

    :code:`Topog function = xsinslope`

        One-dimensional sinusoidal variation along the x-direction, with one
        complete period in the specified domain. Letting :math:`L_{x}` denote
        the domain length in :math:`x`, the surface is
        :math:`b_{0}(x,y) = \epsilon \sin(2\pi x / L_{x}).`

        :code:`Topog params = (epsilon)`

            - :code:`epsilon` -- the amplitude of the sinusoidal variation.

    :code:`Topog function = xysinslope`

        Two-dimensional sinusoidal variation, with one complete period in the
        specified domain. Letting :math:`L_{x}` and :math:`L_{y}` denote the
        domain lengths in :math:`x` and :math:`y` respectively, the suface is
        :math:`b_{0}(x,y) = \epsilon \sin(2\pi x / L_{x}) \sin(2\pi y / L_{y}).`

        :code:`Topog params = (epsilon)`

            - :code:`epsilon` -- the amplitude of the sinusoidal variation.

    :code:`Topog function = xhump`

        One-dimensional cosine hump on a flat topography, 

            :math:`b_{0}(x,y) = \tfrac{1}{2} A \left(1 + \cos(\pi x/L)\right)`,

        for :math:`-L \le x \le L`.

        :code:`Topog params = (A, L)`

            - :code:`A` -- the amplitude of the hump.
            - :code:`L` -- the half-length of the hump.

    :code:`Topog function = xtanh`

        One-dimensional :math:`\tanh` surface,

            :math:`b_{0}(x,y) = A\left[ 1 + \tanh\left((x-x_{0})/L\right) \right]`
        
        :code:`Topog params = (x0, A, L)`

            - :code:`x0` -- the centre of the tanh profile.
            - :code:`A` -- the amplitude of the hump.
            - :code:`L` -- the half-length of the hump.

    :code:`Topog function = xparab`

        One-dimensional parabolic surface,

            :math:`b_{0}(x,y) = Ax^{2}`
        
        :code:`Topog params = (A)`

            - :code:`A` -- coefficient of the parabola.
    
    :code:`Topog function = xyparab`

        Two-dimensional parabolic surface,

            :math:`b_{0}(x,y) = Ax^{2} + By^{2}`
        
        :code:`Topog params = (A, B)`

            - :code:`A` -- coefficient of :math:`x^{2}` for the parabola.
            - :code:`B` -- coefficient of :math:`y^{2}` for the parabola.
