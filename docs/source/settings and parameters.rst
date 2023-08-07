.. _settings_and_parameters:

Settings and parameters
=======================

Kestrel input files are divided into different blocks, which may be specified in
any order. We document the options for each below.

The settings are given in the form :code:`keyword = value`.  Here we provide example values for each keyword.

Some settings are **required** and the simulation will not commence if values are not given.

Other settings are **optional** and default values will be used if they are not given in the input file.  In this case the default value is used as the example.

There are also some settings that are required in only some circumstances (**conditionally required**), and some that are optionally used in some circumstances (**conditionally optional**).

.. include:: settings_and_parameters/domain.rst

.. include:: settings_and_parameters/output.rst

.. include:: settings_and_parameters/parameters.rst

.. include:: settings_and_parameters/solver.rst

.. include:: settings_and_parameters/source.rst

.. include:: settings_and_parameters/topog.rst
