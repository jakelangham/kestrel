.. _set_params:

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

                Sets :math:`\mu_1`
            
                :code:`pouliquen max = 0.4`

                Sets :math:`\mu_2`
            
                :code:`pouliquen beta = 0.136`

                Sets :math:`\beta`

        :code:`variable`

        :code:`voellmy`

    :code:`erosion = on`

.. warning::
    TODO

.. warning::
    TODO: optional settings

    :code:`bed porosity = 0.35`
    
    :code:`eddy viscosity = 0.0`

    :code:`g = 9.81`

    :code:`geometric factors = on`

    .. warning::
        Note that, if this is off, :math:`g_\perp \equiv g`.

    :code:`maxPack = 0.65`

    :code:`rhos = 2000`

    :code:`rhow = 1000`

.. warning::
    TODO: conditionally optional settings

    :code:`erosion transition = smooth`

    :code:`erosion rate = 1e-3`

    :code:`granular erosion rate = 4.0`

    :code:`erosion depth = 1`

    :code:`hindered settling = Spearman Manning`

    :code:`morphodynamic damping = tanh`

    :code:`solid diameter = 1e-3`

    :code:`switch function = tanh`

    :code:`erosion critical height = 0.01`

    :code:`settling speed = 1e-3`

    :code:`voellmy switch rate = 3.0`

    :code:`voellmy switch value = 0.2`

