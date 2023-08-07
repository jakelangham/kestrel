.. _set_params:

Parameters
----------

The *Parameters* block sets up parameters and model closures to be used in Kestrel.  The parameters block is identified using the block keyword :code:`Parameters:`.

The model parameters are closely linked to the model closures used, and therefore many are conditionally required or conditionally optional.  Here we list the parameter inputs by collecting with their associated model closures.  See :ref:`_physical_model_closures` for further details.

.. _set_params_drag:

Basal friction
^^^^^^^^^^^^^^

The basal friction (also referred to as *drag*) is selected through the *required* block variable :code:`Drag`.  The currently implemented drag closures, and their corresponding conditional variables, are:

    :code:`Drag = Chezy`

        The basal friction is modelled using a Chezy form,
            :math:`\mathcal{F} = `