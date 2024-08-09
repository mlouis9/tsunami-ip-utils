# Introduction
This is a python package created to provide enhanced visualization capabilities for criticality safety similarity indices calculated by
[TSUNAMI-IP](https://scale-manual.ornl.gov/tsunami-ip.html#tsunami-ip): a code in the Standardized Computer Analysis for Licensing Evaluation ([SCALE](https://scale-manual.ornl.gov/)) suite. These methods are still under development and are *not* intended for
production use, but can be nice as an exploratory tool.

# Installation and Setup
This package is [published on PyPI](https://pypi.org/project/tsunami-ip-utils/), and so can be installed (along with all of the necessary dependencies) via `pip`
```
pip install tsunami_ip_utils
```
Note many of the features require a valid SCALE installation (SCALE is an export controlled code and can be requested from [RSICC](https://rsicc.ornl.gov/)), and all fo the SCALE binaries to be on PATH.

# Documentation and Examples
If you're interested in running some examples, just clone this repository via the instructions [here](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository#cloning-a-repository) and navigate to the `examples` directory.
More detailed instructions (in a human readable form) can be found in [the documentation](https://mlouis9.github.io/tsunami_ip_utils/),
along with a description of the theory behind the visualization tools, and a fully documented API.
