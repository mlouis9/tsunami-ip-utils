
.. api_type:: public
.. _public_api:

Public API Documentation
========================


Subpackages
-----------

.. toctree::
   :maxdepth: 4

   tsunami_ip_utils.viz

Submodules
----------

tsunami\_ip\_utils.comparisons module
-------------------------------------

.. automodule:: tsunami_ip_utils.comparisons
   :members:
   :undoc-members:
   :show-inheritance:

tsunami\_ip\_utils.config module
--------------------------------

.. automodule:: tsunami_ip_utils.config
   :members:
   :undoc-members:
   :show-inheritance:

tsunami\_ip\_utils.integral\_indices module
-------------------------------------------

.. automodule:: tsunami_ip_utils.integral_indices
   :members:
   :undoc-members:
   :show-inheritance:

tsunami\_ip\_utils.perturbations module
---------------------------------------

.. automodule:: tsunami_ip_utils.perturbations
   :members:
   :undoc-members:
   :show-inheritance:

tsunami\_ip\_utils.readers module
---------------------------------

.. automodule:: tsunami_ip_utils.readers
   :members:
   :undoc-members:
   :show-inheritance:

tsunami\_ip\_utils.utils module
-------------------------------

.. automodule:: tsunami_ip_utils.utils
   :members:
   :undoc-members:
   :show-inheritance:

tsunami\_ip\_utils.xs module
----------------------------

.. automodule:: tsunami_ip_utils.xs
   :members:
   :undoc-members:
   :show-inheritance:


.. _scale-template-input-files:

SCALE Template Input Files
--------------------------
These are the SCALE input files that are used to interface with SCALE via ``scalerte``. Each of the input file is meant to 
be read as a python `Template string <https://docs.python.org/3/library/string.html#template-strings>`_.

.. _mg_reader.inp:

Multigroup SCALE XS Library Reader
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This file is templated with the following arguments:

* multigroup_library_path
    Path to the SCALE multigroup library (e.g. /home/example_user/codes/SCALE-6.3.1/data/scale.rev08.252groupcov7.1).
* nuclide_zaid
    The ZAID of the nuclide to be extracted from the multigroup library ('-1' to select all nuclides).
* reaction_mt
    The reaction MT to be extracted for the selected nuclide(s) ('-1' to select all reactions).
* plot_option
    The charmin option for formatting the dump, e.g. 'fido'.
* output_file_path
    The path to direct the cross section dump to.

Which generates a SCALE input that uses the AMPX utilities tabasco and charmin to read the selected AMPX master library.
Further details on valid options can be found in `the AMPX manual <https://www.ornl.gov/sites/default/files/AMPX-6.pdf>`_.

.. literalinclude:: ../../../src/tsunami_ip_utils/input_files/MG_reader.inp
   :language: text
   :caption: :download:`Download this input file <../../../src/tsunami_ip_utils/input_files/MG_reader.inp>`.
   :linenos:

.. _tsunami_ip_uncertainty_contributions.inp:

Uncertainty Edit Generator
^^^^^^^^^^^^^^^^^^^^^^^^^^
This input file is used to generate the uncertainty edits (i.e. covariance-wise :math:`\Delta k/k` contributions) 
for a given set of SDFs using `TSUNAMI-IP <scale-manual.ornl.gov/tsunami-ip.html#user-input>`_. This file is templated with the following arguments:

* filenames
    A multiline string containing the paths to the SDF files for which the uncertainty edits are to be generated.
* first_file
    An arbitrary selected file to be chosen as an experiment in the TSUNAMI-IP input.

This uses the ``'uncert_long'`` option to generate the uncertainty edits of all of the selected applications (which is just
all of the files supplied).

.. literalinclude:: ../../../src/tsunami_ip_utils/input_files/tsunami_ip_uncertainty_contributions.inp
   :language: text
   :caption: :download:`Download this input file <../../../src/tsunami_ip_utils/input_files/tsunami_ip_uncertainty_contributions.inp>`.
   :linenos:

.. _generate_perturbed_library.inp:

Perturbed Library Generator
^^^^^^^^^^^^^^^^^^^^^^^^^^^
This input file is used for generating a randomly sampled/perturbed multigroup cross section library using the perturbation
factors built-in to SCALE (this should be compared to the analogous feature in 
`SAMPLER <https://scale-manual.ornl.gov/sampler.html#nuclear-data-perturbations-for-multigroup-calculations>`_, which generates
similar input files). This template takes the following arguments:

* base_library
    Path to the base SCALE multigroup library.
* perturbation_factors
    Path to the file containing the perturbation factors for the chosen base multigroup library (e.g. 
    /home/example_user/codes/SCALE-6.3.1/data/perturb/252n.v7.1/Sample20).
* sample_number
    The sample number to use. Note this **must** correspond with the chosen perturbation_factors or the library will not
    be generated correctly (e.g. the correct sample number, corresponding with the example perturbation factors is 20).
* output
    Path to copy the generated perturbed library to.

.. literalinclude:: ../../../src/tsunami_ip_utils/input_files/generate_perturbed_library.inp
   :language: text
   :caption: :download:`Download this input file <../../../src/tsunami_ip_utils/input_files/generate_perturbed_library.inp>`.
   :linenos:


Module contents
---------------

.. automodule:: tsunami_ip_utils
   :members:
   :undoc-members:
   :show-inheritance:
