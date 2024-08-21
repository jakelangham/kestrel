Parameters
----------

The *Parameters* block specifies the model closures and associated parameters
for the simulation. It is identified using the block keyword
:code:`Parameters:`.

There are two **required** settings in the *Parameters* block:

    :code:`drag = chezy`

    This sets the basal drag closure function :math:`\tau_b`, which we assume to
    be of the form :math:`\tau_b = \frac{\bar{\rho}
    \bar{\mathbf{u}}}{\mathbf{u}}\mathcal{F}`, where (as detailed in
    :ref:`physical_model`) :math:`\mathcal{F}` is a friction function that may
    depend on the local flow variables and any conditionally optional parameters
    relevant to the selected drag law. Options are:

        :code:`chezy`

            Sets the Chézy drag law

                :math:`\mathcal{F} = C_d |\mathbf{u}|^2`

            for modelling the basal resistance of a turbulent fluid.

            The **conditionally optional** setting

                :code:`chezy co = 0.01`

            is used to set the value of the constant drag coefficient
            :math:`C_d`. It is required to be strictly positive.

        :code:`coulomb`

            Sets the Coulomb drag law

                :math:`\mathcal{F} = \mu g_\perp H`

            for modelling granular friction. Here, :math:`g_\perp` denotes
            gravitational acceleration resolved perpendicular to the local
            slope.

            The **conditionally optional** setting

                :code:`coulomb co = 0.1`

            is used to set the value of the constant friction coefficient
            :math:`\mu`. It is required to be strictly positive.

        :code:`manning`

            Sets the Manning drag law
            
                :math:`\mathcal{F} = \frac{g_\perp n^2}{H^{1/3}}`,

            where :math:`g_\perp` is gravitational acceleration resolved
            perpendicular to the local slope.

            The **conditionally optional** setting

                :code:`manning co = 0.03`
     
            is used to set the *dimensional* constant coefficient :math:`n`.  It
            is required to be strictly positive.

        :code:`pouliquen`

            Sets Pouliquen and Forterre's granular drag law (see e.g. Pouliquen
            & Forterre, *J. Fluid Mech.*, **453**, `2002
            <https://doi.org/10.1017/S0022112001006796>`_)

                :math:`\mathcal{F} = \mu(I) g_\perp H`,

            where :math:`g_\perp` is gravitational acceleration resolved
            perpendicular to the local slope and :math:`I` is a dimensionless
            *inertial number*, defined by :math:`I = d |\mathbf{u}| /
            \sqrt{g_\perp H^3}`. (The parameter :math:`d` is the characteristic
            sediment diameter, which may be set via the ``solid diameter``
            option.) The friction coefficient is given by

                :math:`\mu(I) = \mu_1 + \frac{\mu_2 - \mu_1}{1 + \beta / I}`

            The following **conditionally optional** settings are used to set
            the various empirical coefficients in this model:

                :code:`pouliquen min = 0.1`

                Sets :math:`\mu_1`, which is required to be strictly positive.
            
                :code:`pouliquen max = 0.4`

                Sets :math:`\mu_2`, which is required to be strictly positive.
            
                :code:`pouliquen beta = 0.136`

                Sets :math:`\beta`, which is required to be strictly positive.

        :code:`edwards2019`

            Sets Edwards et al.'s extension to the Pouliquen law (see e.g. 
            Edwards, Russell, Johnson, Gray, *J. Fluid Mech.*, **875**, `2019
            <https://doi.org/10.1017/jfm.2019.517>`_)

                :math:`\mathcal{F} = \mu(Fr, H) g_\perp H`,

            where :math:`g_\perp` is gravitational acceleration resolved
            perpendicular to the local slope, :math:`Fr = |\mathbf{u}|/\sqrt{g_\perp H}` 
            is the local Froude number and the friction function is defined as

                :math:`\mu(Fr, H) = \mu_1 + \frac{\mu_2 - \mu_1}{1 + H \beta / (L(Fr + \Gamma))}`,

            when :math:`Fr > \beta_*` and

                :math:`\mu(Fr, H) = \left(\frac{Fr}{\beta_*}\right)^\kappa \left\{\mu_1 + \frac{\mu_2 - \mu_1}{1 + H \beta / (L(\beta_* + \Gamma))} - \mu_{\mathrm{start}(H)}\right\} + \mu_{\mathrm{start}(H)}`

            when :math:`Fr \leq \beta_*`, with :math:`\mu_{\mathrm{start}(H)} =
            \mu_3 + (\mu_2 - \mu_1) / (1 + H / L)`.

            The following **conditionally optional** settings are used to set
            the various empirical coefficients in this model:

                :code:`pouliquen min = 0.1`

                Sets :math:`\mu_1`, which is required to be strictly positive.
            
                :code:`pouliquen max = 0.4`

                Sets :math:`\mu_2`, which is required to be strictly positive.

                :code:`pouliquen intermediate = 0.2`
            
                Sets :math:`\mu_3`, which is required to be strictly positive.

                :code:`pouliquen beta = 0.136`

                Sets :math:`\beta`, which is required to be strictly positive.

                :code:`edwards2019 betastar = 0.136`
                
                Sets :math:`\beta_*`, which is required to be strictly positive.

                :code:`edwards2019 kappa = 1.0`

                Sets :math:`\kappa`, which is required to be strictly positive.

                :code:`edwards2019 gamma = 0.0`

                Sets :math:`\Gamma`.

        :code:`variable`

            Sets the following drag law, which interpolates between the Chézy
            and Pouliquen laws, depending on the solids fraction:

                :math:`\mathcal{F} = (1 - f(\bar{\psi})) C_d |\mathbf{u}|^2 +
                f(\bar{\psi}) \mu(I) g_\perp H`,

            where :math:`f` is a switching function, equal to

                :math:`f(\bar{\psi})=\frac{1}{2}[1+\tanh(V_R(\bar{\psi}-\psi^*))]`.

            Its parameters may be set via the **conditionally optional**
            statements

                * :code:`voellmy switch rate = 3.0`, which sets :math:`V_R`
                * :code:`voellmy switch value = 0.2`, which sets
                  :math:`\psi^*`

        :code:`voellmy`

            Sets the Voellmy drag closure, which is the sum of Chézy and Coulomb
            drag:

                :math:`\mathcal{F} = C_d |\mathbf{u}|^2 + \mu g_\perp H`.

            As above, :math:`C_d` and :math:`\mu` may be set by the
            **conditionally optional** statements :code:`chezy co` and
            :code:`coulomb co` respectively.

    :code:`erosion = on`

        This sets the closure function :math:`\mathcal{E}` for erosion (see
        :ref:`physical_model`. If set to anything other than :code:`erosion =
        off`, it also activates Kestrel's morphodynamic capabilities. The
        following options are available:

        :code:`fluid`

            Sets a 'fluid-like' erosion, with

                :math:`\mathcal{E} = \max\{ \varepsilon_f u_p (\theta - \theta_c), 0\}`,

            where :math:`u_p = \sqrt{g_\perp d (\rho_s/\rho_f - 1)}` with
            :math:`\rho_s, \rho_f` denoting sediment and fluid densities
            respectively.  :math:`\theta = C_d |\mathbf{u}|^2/u_p^2` denotes the
            Shields number for the flow. (N.b. :math:`g_\perp` and :math:`d` are
            defined in the discussion of the :code:`drag` setting.)

            Erosion occurs when :math:`\theta` exceeds the critical value
            :math:`\theta_c`, determined by the empirical closure

                :math:`\theta_c = \frac{0.3}{1 + 1.2 R} + 0.055[1-\exp(-0.02R)]`,

            where :math:`R = d [g (\rho_s/\rho_f - 1) / \nu_w^2]^{1/3}` and
            :math:`\nu_w = 1.2\times 10^{-6}\textrm{m}^2/\textrm{s}` is the
            kinematic viscosity of water. (cf. Soulsby, *Dynamics of Marine
            Sands*, 1997)

            Three **conditionally optional** statements affect this closure:
            
                * The constant coefficient :math:`\varepsilon_f` may be set via the

                    :code:`erosion rate = 1e-3`

                  whose value is required to be strictly positive.

                * Sediment and fluid densities may be set via

                    :code:`rhos = 2000`

                  and

                    :code:`rhow = 1000`

                  respectively. These are required to be strictly positive.

        :code:`granular`

            Sets a 'granular-like' erosion, with

                :math:`\mathcal{E} = \max\{ \varepsilon_g u_p [\mu(I)g_\perp H - \mu_n], 0\}`,

            where :math:`\mu(I)` is Pouliquen's friction coefficient (see the
            :code:`drag` discussion above) and 

                :math:`\mu_n = \mu_1 + \frac{\mu_s - \mu_1}{1 +
                \left(\frac{H}{25d}\right)^2}`

            with :math:`\mu_s = [\mu_1 + \tan(1^\circ)] / [1 - \mu_1
            \tan(1^\circ)]` (:math:`\mu_1` is set by :code:`pouliquen min`).

            The **conditionally optional** declaration

                :code:`granular erosion rate = 1e-3`

            may be used to set the constant coefficient :math:`\varepsilon_g`,
            which is required to be strictly positive.

        :code:`mixed` or :code:`on` (default)

            Sets the following erosion law that switches between fluid-like and
            granular-like erosion rates, depending on the solids fraction:

                :math:`\mathcal{E} = (1 - f(\bar{\psi})) \mathcal{E}_f +
                f(\bar{\psi}) \mathcal{E}_g`,

            where :math:`\mathcal{E}_f` and :math:`\mathcal{E}_g` are the
            corresponding erosion rates according to the :code:`fluid` and
            :code:`granular` closures respectively. The function :math:`f` is
            the same switching function as in the case of :code:`drag =
            variable`.

        :code:`off`

            Sets :math:`\mathcal{E} = 0`.

        :code:`simple`

            Sets a simple model for erosion based on the Shields number, with no
            critical value:

                :math:`\mathcal{E} = \varepsilon u_p \theta`.

            The constant coefficient :math:`\varepsilon` may be defined via the
            **conditionally optional** setting 

                :code:`erosion rate = 1e-3`

            which is required to be strictly positive.

The remaining settings in the *Parameters* block are **optional**. We list them
below:

    :code:`bed porosity = 0.35`

        Sets the bed porosity :math:`p`. This is related to the solid fraction
        :math:`\psi_b` of the bed by :math:`\psi_b = 1 - p` and as such, affects
        the rate of sediment transfer between the flow and bed (see
        :ref:`physical_model`). Kestrel requires :math:`0 < p \leq 1`.

        .. note::
            In most cases, it is prudent to have ``bed porosity`` equal to ``1 -
            maxPack``.
 
    :code:`deposition = Spearman Manning`

        This sets the deposition rate closure :math:`\mathcal{D}`. The following
        options are available:

            :code:`none`

                Sets :math:`\mathcal{D} = 0`.

            :code:`simple`

                Sets a simple quadratic hindered settling law of the form

                    :math:`\mathcal{D} = w_s \bar{\psi}(1 -
                    \bar{\psi}/\psi_{\max})`,

                where :math:`w_s` is characteristic sediment settling speed
                and :math:`\psi_{\max}` is the maximum volume fraction that the
                flowing sediment may be packed into. These constant coefficients
                may be set via the **conditionally optional** declarations:

                    :code:`settling speed = 1e-3`

                    sets :math:`w_s`. If not explicitly set, Kestrel defaults to
                    using an empirical law based on the solid diameter
                    :math:`d`:

                        :math:`w_s = \frac{\nu_w}{d}\left\{\sqrt{10.36^2 + 1.048R} - 10.36\right\}`,

                    where (as above) :math:`R = d [g (\rho_s/\rho_f - 1) / \nu_w^2]^{1/3}` and
                    :math:`\nu_w = 1.2\times 10^{-6}\textrm{m}^2/\textrm{s}` is the
                    kinematic viscosity of water.

                    :code:`maxPack = 0.65`

                    sets :math:`\psi_{\max}`.

            :code:`Spearman Manning`

                Sets an empirical hindered settling law due to Spearman &
                Manning (*Ocean Dynam.* **67(3)**, 2017):

                    :math:`\mathcal{D} = w_s \bar{\psi} (1 - \bar{\psi})^a (1 -
                    \bar{\psi}/\psi_{\max})^b`

                The exponents :math:`a` and :math:`b` are determined via the
                formulae :math:`a = 2.7 - 0.15 n` and :math:`b = 0.62n - 1.46`
                where

                    :math:`n = \frac{4.7 + 0.41 (u_p d / \nu_w)^{3/4}}{1 + 0.175
                    (u_p d / \nu_w)^{3/4}}`

                (Rowe, P. N, *Chem. Eng. Sci*, **42**, 1987). The constants
                :math:`w_s` and :math:`\psi_{\max}` may be set via
                :code:`settling speed` and :code:`maxPack` resp. (as above).

    :code:`eddy viscosity = 0.0`

        Sets the (constant) value of eddy viscosity :math:`\nu` in the model
        (see :ref:`physical_model`). This value is required to be non-negative.

        .. warning::
            If you want to simulate morphodynamics then eddy viscosity must be
            non-zero. Otherwise, the underlying governing equations are
            ill-posed as an initial value problem and Kestrel's numerical
            solutions will fail to converge as the grid resolution is refined.

    :code:`erosion critical height = 0.01`

        Sets a critical flow depth :math:`H_c` in metres, below which erosion is
        not permitted. This ensures that rapid thin flows do not unphysically
        erode the bed. It is recommended that this is at least equal to the
        characteristic solid diameter :math:`d`. It is required to be strictly
        positive.

        A phenomenological function :math:`\chi(H)` is pre-multiplied to the
        value of the morphodynamic transfer rates :math:`\mathcal{E}` and
        :math:`\mathcal{D}` to achieve this. It is user-selectable via the
        **conditionally optional** setting

        :code:`morphodynamic damping`

        This has options

            * ``off``. Sets :math:`\chi = 1`.
            * ``rat3``. Sets :math:`\chi = 0` if :math:`H < H_c`,
              :math:`\chi = 1` if :math:`H > 2 H_c`, :math:`\chi = (\frac{H}{H_c} - 1)^3 [(2 - \frac{H}{H_c})^{3} + (\frac{H}{H_c} - 1)]^{-1}` otherwise.
            * ``tanh`` (default). Sets :math:`\chi = \frac{1}{2}[1 + \tanh(10
              \log(H/H_c))]`.

    :code:`erosion depth = 1`

        Sets the depth in metres up to which erosion is permitted. In the
        notation of :ref:`physical_model`, this means setting :math:`\Delta
        b_{\max}`. It is required to be non-negative.

        A phenomenological function :math:`\Theta(\Delta b)` (where
        :math:`\Delta b \equiv b(\mathbf{x},t) - b_0(\mathbf{x},0)`) is
        pre-multiplied to the value of the erosion rate :math:`\mathcal{E}` to
        achieve this. It is user-selectable via the **conditionally optional**
        setting

        :code:`erosion transition`

        This has options

            * ``off``. Sets :math:`\Theta = 1`.
            * ``smooth`` (default). Sets :math:`\Theta = \frac{1}{2}[1 +
              \tanh(10^5(\Delta b + \Delta b_{\max}))]`.
            * ``step``. Sets :math:`\Theta = 0` if :math:`\Delta b < -\Delta
              b_{\max}`, :math:`\Theta = 1` otherwise.

    :code:`g = 9.81`

        Sets the gravitational acceleration :math:`g`.

    :code:`geometric factors = on`

        This option selects whether the model equations that Kestrel solves
        should consider geometric corrections that arise due to variations in
        the topographic surface. This is the default model, described in
        :ref:`physical_model`. If ``geometric factors = off``, the only effect
        of slope variation left in the model is that of gravitational forcing
        along the direction of steepest descent. This is equivalent to setting
        :math:`\gamma \equiv 1` in :ref:`physical_model`.

        For more information, see this `reference
        <https://arxiv.org/abs/2306.16185>`_.

        .. warning::
            Note that, if :code:`geometric factors = off`, then :math:`g_\perp
            \equiv g`.

    :code:`solid diameter = 1e-3`

        Sets the characteristic sediment diameter :math:`d`. This affects
        various optional closures, that model the physics of grains such as the
        settling speed, ``Pouliquen`` drag rule and ``granular`` erosion. It is
        required to be strictly positive.
 
