Solver
------

The *Solver* block sets up the core numerical solver in Kestrel.  The solver block is identified using the block keyword :code:`Solver:`.

The only **required** setting in the solver block is:

    :code:`T end`

        The time, in seconds, at which to end the simulation.


The **optional** settings of the solver block are:

    :code:`limiter = MinMod2`

        The slope limiter using to approximate gradients in the governing equations.  Options are:

            :code:`MinMod1`

                The minmod limiter, defined as

                    .. math::
                        \mathrm{MinMod}(a, b) = \tfrac{1}{2}(\mathrm{sgn}(a) + \mathrm{sgn}(b)).\min(\left|a\right|, \left|b\right|).

                for real numbers :math:`a` and :math:`b` and :math:`\mathrm{sgn}(a)` denotes the signum function.
            
            :code:`MinMod2` (default)

                The generalised minmod limiter, defined as

                    .. math::
                        \mathrm{MinMod}(a, b) = \begin{cases} \mathrm{sgn}(a) . \min(\theta\left|a\right|, \tfrac{1}{2}\left|a+b\right|, \theta\left|b\right|), & \text{if $\mathrm{sgn}(a)=\mathrm{sgn}(b)$}\\ 0, & \text{if $\mathrm{sgn}(a)+\mathrm{sgn}(b)$ = 0}. \end{cases}
                
                The parameter :math:`1\le \theta \le 2` and we impose :math:`\theta = 1.3` which we have found to work well (see also `Chertock et al. (2015) <https://doi.org/10.1002/fld.4023>`_).

                    .. note::
                        :code:`MinMod2` is strongly recommended.
            
            :code:`van albada` or :code:`albada`

                The van Albada limiter is a smooth function that tends to a
                central difference approximation in smooth regions, and across
                discontinuities the slope is biased to smallest value of the two
                one-sided slopes. In general, the van Albada limiter has a
                parameter, but here we use a version with parameter set to zero
                (see `Lu & Tadmor (2023)
                <https://doi.org/10.48550/arXiv.2304.00437>`_)

                    .. math::
                        \mathrm{vanAlbada}(a,b) = \frac{a^2 b + a b^2}{a^2 + b^2}.
                
            :code:`WENO`

                A weighted essentially non-oscillatory type limiter.

                    .. math::
                        \mathrm{WENO}(a,b) = \frac{w(a)a + w(b)b}{w(a) + w(b)} \text{ with weighting } w(x) = (x^2+\epsilon)^{-2}

                Kestrel imposes :math:`\epsilon = 10^{-6}`.
            
            :code:`None`

                No slope limiter used.  Gradients are approximated using central differences.

                .. note::
                        As the governing equations are hyperbolic, solutions can exhibit discontinuities.  Therefore, :code:`None` is **not** recommended.

    :code:`height threshold = 1e-6`

        A threshold on flow depths.  Flow depths below :code:`height threshold` are neglected.

    :code:`Tile buffer = 1`

        Neighbouring tiles are activated if flow reaches :code:`tile buffer` cells from a tile edge.  The default value :code:`tile buffer = 1` should ensure that a neighbouring tile is added when needed.

    :code:`CFL = 0.5` for 1D simulations; :code:`CFL = 0.25` for 2D simulations

        The Courant-Friedrichs-Lewy (CFL) number of the simulation.  This determines the maximum time step.

            - For 1D simulations the scheme requires :code:`CFL` :math:`\le 0.5`.
            - For 2D simulations the scheme requires :code:`CFL` :math:`\le 0.25`.

    :code:`max dt = HUGE(1.0)`

        The maximum time step. Note that the default value ensures that the
        maximum time step is the either the time step determined by the CFL
        condition, or the output time step.
    
    :code:`T start = 0`

        The start time of the simulation, in seconds.
    
    :code:`Restart = off`

        Should the simulation restart from a previous result? Options are:

            :code:`on` -- restart from a previous simulation.  This requires an
            existing Kestrel results directory containing the RunInfo.txt file,
            which is used to determine suitable simulation parameters and to
            locate the last output file to be used as the initial condition.

            :code:`off` -- start a new simulation.

            .. note::

                This feature is useful if a simulation is interrupted for some
                reason. By selecting ``on`` and re-running, the simulation will
                pick up from where it left off.
    
    :code:`Initial condition`

        Specifies the path to a Kestrel result file to be used as an initial
        condition. On start-up Kestrel loads the solution fields from this file
        and simulates forward from this point.

        If ``Restart = on``, then `RunInfo.txt` is used to determine the
        simulation parameters. Otherwise, (by default) they are read in from the
        usual input file given on the command line.

The following **conditionally optional** variables are used only if
:code:`Boundary conditions = sponge` in the *Domain* block:

    :code:`Sponge strength = 0.2`

        When using a sponge layer boundary condition, the solution's quantities
        are gradually damped on the tiles bordering the domain boundary.  The
        damping rate is set by the :code:`Sponge strength` settings.

        .. note::
            Care must be taken in setting :code:`Sponge strength`.  If the
            damping is too weak, flow quantities may be non-zero at the domain
            edge, causing errors.  If the damping is too strong, flow quantities
            in the interior can be influenced by those in the sponge layer
            tiles.  The flow in the sponge layer tiles should be discarded.
