#!/usr/bin/env python
from __future__ import print_function
import os
import datetime
import numpy as np
import netCDF4 as nC4
import matplotlib.pyplot as plt
# - Test pyTMD Installation
from pyTMD.read_tide_model import extract_tidal_constants

print('# - All Dependencies Correctly Imported.')