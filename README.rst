=======================================
Compute Tide Time Series with pyTMD
=======================================
|Language|
|License|

.. |Language| image:: https://img.shields.io/badge/python%20-3.7%2C%203.8%2C%203.9-brightgreen?style=plastic
   :target: https://www.python.org/

.. |License| image:: https://img.shields.io/badge/license-MIT-green.svg
   :target: https://github.com/eciraci/Download_ECMWF_Data/blob/main/LICENSE


This repository contains small collection of scripts that can be used to Compute
Tide Time Series by employing the python pyTMD module.

- PyTMD Project Homepage: https://pytmd.readthedocs.io/en/latest/

- PyTMD GitHub: https://github.com/tsutterley/pyTMD

*Other Resources*:

- `ESR Polar Tide Models <https://www.esr.org/research/polar-tide-models/list-of-polar-tide-models/>`_
- `Finite Element Solution (FES) tide models <https://www.aviso.altimetry.fr/en/data/products/auxiliary-products/global-tide-fes.html>`_


\
\


**Installation**:

1. Setup minimal **conda** installation using  `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_
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


\
\


**Execution**:

1. compute_pt_tide_time_series.py
    Compute Tide Time Series at a selected geographic location

    .. code-block:: bash

        python compute_pt_tide_time_series.py ./parameters.yml -68.72 -5.32

2. compute_shp_tide_time_series.py
    Compute Hourly Tide Time Series at the selected locations listed inside a ESRI shapefile

    .. code-block:: bash

        python compute_shp_tide_time_series.py ./parameters.yml

3. **parameters.yaml** - Parameters File Content.

.. code-block:: yaml

    model: Model - [CATS2008, FES2014, AOTIM5]
    path: Absolute Path to Model's data directory
    shp: ./data/Antarctic_GZ_tide/Antarctic_GZ_tide.shp
    out_path: Absolute Path to Output directory
    out_format: ascii/txt/netcdf
    date1: Initial Date MM/DD/YYYY
    date2: Final Date MM/DD/YYYY
    nproc: Number of Maximum simultaneous processes - integer



**PYTHON DEPENDENCIES**:
#######
 - `pyTMD: Python-based tidal prediction software that reads OTIS, GOT and FES formatted tidal solutions.  <https://github.com/tsutterley/pyTMD>`_
 - `numpy: The fundamental package for scientific computing with Python. <https://numpy.org>`_
 - `xarray: xarray: N-D labeled arrays and datasets in Python. <https://xarray.pydata.org/en/stable>`_
 - `pandas: Python Data Analysis Library. <https://pandas.pydata.org>`_
 - `geopandas: Python tools for geographic data. <https://geopandas.org/en/stable/>`_
 - `rasterio: access to geospatial raster data. <https://rasterio.readthedocs.io>`_
 - `fiona: reads and writes geographic data files. <https://fiona.readthedocs.io>`_
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