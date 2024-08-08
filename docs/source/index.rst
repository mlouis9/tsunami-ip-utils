.. tsunami_ip_utils documentation master file, created by
   sphinx-quickstart on Sun Jul  7 12:00:07 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to tsunami_ip_utils's documentation!
============================================

This is the documentation for tsunami_ip_utils, a python package for providing enhanced visualization capabilities for criticality safety similarity indices calculated by
`TSUNAMI-IP <https://scale-manual.ornl.gov/tsunami-ip.html#tsunami-ip>`_: a code in the Standardized Computer Analysis for Licensing 
Evaluation (`SCALE <https://scale-manual.ornl.gov/>`_) suite. These methods are still under development and are *not* intended for
production use, but can be nice as an exploratory tool.

**NOTE**: This package requires SCALE to be installed and the binaries to be on the system path. SCALE is an export controlled code
and can be requested from `RSICC <https://rsicc.ornl.gov/>`_.

The documentation is organized into a public API (i.e. all functions that are intended to be used by users, and whose interface will
not change without a deprecation warning) and a private API (i.e. all functions that are not intended to be used by users, and whose
interface may change without a deprecation warning). The public API is the recommended way to interact with the package, while the
private API is intended for developers who wish to extend the package.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Public API <public_api/index>
   Private API <private_api/index>
   Paths <paths>
   Examples <auto_examples/index>
   Report <report/index>

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
