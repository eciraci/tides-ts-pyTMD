=======================================
Compute Tide Time Series with pyTMD
=======================================
|Language|
|License|

.. |Language| image:: https://img.shields.io/badge/python%20-3.7%2C%203.8%2C%203.9-brightgreen?style=plastic
   :target: https://www.python.org/

.. |License| image:: https://img.shields.io/badge/license-MIT-green.svg
   :target: https://github.com/eciraci/Download_ECMWF_Data/blob/main/LICENSE


Installation
=============
1. Setup minima conda installation using  `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_
2. Create Python Virtual Environment

    - Creating an environment with commands (`Link <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands>`_);
    - Creating an environment from an environment.yml file (`Link <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file>`_);

3. Install Python Dependencies

    .. code-block:: bash

       conda install -y pip numpy matplotlib cartopy pandas geopandas jupyter netcdf4 xarray
       conda install -c conda-forge pytmd

    or run *./install_dependencies.sh*
4. Test Python Dependencies Installation

    .. code-block:: bash

        python test_requirements_installation.py



**PYTHON DEPENDENCIES**:
#######
 - `numpy: The fundamental package for scientific computing with Python <https://numpy.org>`_
 - `xarray: xarray: N-D labeled arrays and datasets in Python <https://xarray.pydata.org/en/stable>`_
 - `pandas: Python Data Analysis Library <https://pandas.pydata.org>`_
 - `geopandas: Python tools for geographic data <https://geopandas.org/en/stable/>`_
 - `rasterio: access to geospatial raster data <https://rasterio.readthedocs.io>`_
 - `fiona: reads and writes geographic data files <https://fiona.readthedocs.io>`_
 - `shapely: Manipulation and analysis of geometric objects in the Cartesian plane. <https://shapely.readthedocs.io/en/stable>`_
 - `cartopy: Python package designed to produce maps and other geospatial data analyses. <https://scitools.org.uk/cartopy>`_
 - `matplotlib: Library for creating static, animated, and interactive visualizations in Python. <https://matplotlib.org>`_

\
\
License
#######

The content of this project is licensed under the
`Creative Commons Attribution 4.0 Attribution license <https://creativecommons.org/licenses/by/4.0/>`_
and the source code is licensed under the `MIT license <LICENSE>`_.