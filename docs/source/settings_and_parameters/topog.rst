Topog
-----

The *Topog* block specifies the topography to be used for the simulation.
Topography can be given either as a Digital Elevation Model (DEM) or through a
parameterised analytic function.  The *Topog* block is identified using the block
keyword :code:`Topog:`.

There is a single **required** parameter in the *Topog* block:

    :code:`Type`

        The type of topography to use. Options are:

        :code:`DEM` or :code:`Raster`

            Use a georeferenced DEM contained in a geotiff file.
        
        :code:`SRTM`

            Use a georeferenced SRTM file, contained in an SRTM archive.

        :code:`Function`

            Use a parameterised analytic function defining a surface.

The following statement is **conditionally required** when :code:`Type = DEM` or
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
        an SRTM archive, i.e. either:
        - containing zip files with extension :code:`.SRTMGL1.hgt.zip`
        - a directory with a structure N01/N01*.tif, N02/N02*.tif, ...

        .. note::
            SRTM files may be either geotiff or USGS hgt format (including zipped
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
            - :code:`R` -- radius of connecting circular arc.

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

    :code:`Topog function = xTrislope`

        A surface with three constant slopes connected piecewise continuously 
        by trigonometric functions over a characteristic length scale
        :math:`\lambda`.

        :code:`Topog params = (phi1, phi2, phi3, lambda, x1, x2)`

            - :code:`phi1` -- the slope angle :math:`\phi_1` for :math:`x<(x_1-\lambda/2)`, in degrees.  A positive value corresponds to an elevation decreasing from left to right.
            - :code:`phi2` -- the slope angle :math:`\phi_2` for :math:`(x_1-\lambda/2)<x<(x_2-\lambda/2)`, in degrees.  A positive value corresponds to an elevation decreasing from left to right.
            - :code:`phi3` -- the slope angle :math:`\phi_3` for :math:`(x_2-\lambda/2)<x`, in degrees.  A positive value corresponds to an elevation decreasing from left to right.
            - :code:`lambda` -- the length scale :math:`\lambda` over which transitions between
              constant slope values occur.
            - :code:`x1` -- the point :math:`x_1` where the transition between
              :math:`\phi_1` and :math:`\phi_2` is centred.
            - :code:`x2` -- the point :math:`x_2` where the transition between
              :math:`\phi_2` and :math:`\phi_3` is centred.

    :code:`Topog function = USGS`

        Parametrisation of the USGS flume.  This has slope of 31° for
        :math:`x<0`, and slope 2.4° for :math:`x>x_{1}>0` that are connected
        by a smooth :math:`\cosh` curve section.  Note :math:`x_{1}` is
        determined to ensure smooth connection of the slope elements.  The flume
        is confined by walls for :math:`x<8.5` m, that are represented as
        :math:`\tanh` profile humps.  See `Iverson et al. (2010)
        <https://doi.org/10.1029/2009JF001514>`_ for details.

        :code:`Topog params = (wallH, sigma)`

            - :code:`wallH` -- the height of the sidewalls of the flume.
            - :code:`sigma` -- the characteristic lateral slope of the channel
              banks, which are centred at :math:`y=\pm 3/2`. 
              Note that sufficiently low :code:`sigma` values also
              control the effective width of the channel.

    :code:`Topog function = flume`

        Generalised USGS flume geometry that requires six parameters to specify
        its features. It defines a channel at slope angle :math:`\theta_0` for
        :math:`x<0`, confined by walls for :math:`x<x_{\mathrm{wall}}` that
        smoothly connects to an unconfined runout plane of slope angle
        :math:`\theta_1` for :math:`x>x_1>0` (note that :math:`x_1` is
        determined to ensure the smooth connection).
        The confining walls are constructed by :math:`\tanh` bumps, centred at
        :math:`y=\pm W` with characteristic width :math:`W`.

        :code:`Topog params = (theta0, theta1, xwall, wallW, wallH, sigma)`

            - :code:`theta0` -- the slope (in degrees) for :math:`x < 0`.
            - :code:`theta1` -- the slope (in degrees) for :math:`x > x_1`.
            - :code:`xwall` -- :math:`x_{\mathrm{wall}}` coordinate defining the
              transition between confined and unconfined parts of the topography.
            - :code:`wallW` -- the characteristic width :math:`W` of the channel.
            - :code:`wallH` -- the height of the sidewalls of the flume.
            - :code:`sigma` -- the characteristic lateral slope of the channel
              banks, which are centred at :math:`y=\pm W`.
              Note that sufficiently low :code:`sigma` values also
              control the effective width of the channel.

    :code:`Topog function = channel power law`
    
        Channel with constant slope in :math:`x` and banks defined by a power law as so

        :math:`b(x, y) = Sx + \cos(\theta)|y/W|^\alpha`,

        where :math:`S = \tan(\theta)` is the slope in :math:`x` and :math:`W`, 
        :math:`\alpha` are parameters defining a power law cross-section.

        :code:`Topog params = (slope, W, alpha)`

            - :code:`slope` -- the slope :math:`S` in :math:`x`.
            - :code:`W` -- the characteristic width :math:`W` of the power law
              cross-section.
            - :code:`alpha` -- index :math:`\alpha` of the power law.

    :code:`Topog function = channel trapezium`

        Channel with constant slope in :math:`x` and banks defined by a trapezium as so

        :math:`b(x, y) = Sx + \cos(\theta)\max\{0, S_b (|y| - W/2)\}`,

        where :math:`S = \tan(\theta)` is the slope in :math:`x`, :math:`W` is
        the width of the trapezoid cross-section and :math:`S_b` is the slope of
        the banks in a frame oriented along the channel.

        :code:`Topog params = (slope, W, Sb)`

            - :code:`slope` -- the slope :math:`S` in :math:`x`.
            - :code:`W` -- the width of the channel base.
            - :code:`Sb` -- the gradient of the banks in slope-aligned
              coordinates.

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
        domain lengths in :math:`x` and :math:`y` respectively, the surface is
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
