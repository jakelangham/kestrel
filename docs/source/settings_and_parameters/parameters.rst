.. _set_params:

Parameters
----------

The *Parameters* block specifies the model closures and associated parameters
for the simulation. It is identified using the block keyword
:code:`Parameters:`.

There are two **required** settings in the *Parameters* block:

    :code:`drag = chezy`

.. warning::
    TODO

and

    :code:`erosion = on`

.. warning::
    TODO

.. warning::
    TODO: optional settings

    :code:`bed porosity = 0.35`
    
    :code:`eddy viscosity = 0.0`

    :code:`g = 9.81`

    :code:`geometric factors = on`

    :code:`maxPack = 0.65`

    :code:`rhos = 2000`

    :code:`rhow = 1000`

.. warning::
    TODO: conditionally optional settings

    :code:`chezy co = 0.01`

    :code:`coulomb co = 0.1`

    :code:`manning co = 0.03`
     
    :code:`pouliquen min = 0.1`

    :code:`pouliquen max = 0.4`

    :code:`pouliquen beta = 0.136`

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