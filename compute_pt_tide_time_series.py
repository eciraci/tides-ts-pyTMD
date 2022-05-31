#!/usr/bin/env python
u"""
compute_pt_tide_time_series.py
Written by Enrico Ciraci' (05/2022)

Compute Hourly Tide Time Series at the selected location.

usage: compute_pt_tide_time_series.py [-h] parameters latitude longitude

positional arguments:
  parameters  Processing Parameters File [yml - format].
  latitude    Point Latitude [deg].
  longitude   Point Longitude [deg].

optional arguments:
  -h, --help  show this help message and exit


PYTHON DEPENDENCIES:
    argparse: Parser for command-line options, arguments and sub-commands
           https://docs.python.org/3/library/argparse.html
    numpy: The fundamental package for scientific computing with Python
           https://numpy.org/
    matplotlib: library for creating static, animated, and interactive
           visualizations in Python.
           https://matplotlib.org
    datetime: Basic date and time types
           https://docs.python.org/3/library/datetime.html#module-datetime
    pyTMD: Python-based tidal prediction software that reads OTIS, GOT and FES
           https://github.com/tsutterley/pyTMD
    xarray: xarray: N-D labeled arrays and datasets in Python
           https://docs.xarray.dev/en/stable/
UPDATE HISTORY:
"""
# - Python Dependencies
from __future__ import print_function
import os
import argparse
import datetime
import yaml
import numpy as np
import xarray as xr
# - pyTMD
from pyTMD.read_tide_model import extract_tidal_constants
from pyTMD.infer_minor_corrections import infer_minor_corrections
from pyTMD.predict_tide import predict_tide
from pyTMD.read_FES_model import extract_FES_constants
# - Utility functions
from utils.create_dir import create_dir


def compute_pt_tide_ts(pt_lon: float, pt_lat: float,
                       date1: datetime.datetime, date2: datetime.datetime,
                       tide_model: str, tide_model_path: str,
                       verbose: bool = False) -> dict:
    """
    Compute Tide Time Series at the selected Geographic Coordinates.
    :param pt_lon: Point Longitude - deg
    :param pt_lat: Point Latitude - deg
    :param date1: Initial Time - datetime.datetime
    :param date2: Final Time - datetime.datetime
    :param tide_model: Tidal Model - ['CATS2008', 'AOTIM5', 'FES2014']
    :param tide_model_path: Absolute Path to Tide Model Data
    :param verbose: Print intermediate results on Standard Output.
    :return: Python Dictionary Containing the computed Time Series.
    """
    # - Define Model's Specific Processing Parameters
    if tide_model in ['CATS2008', 'AOTIM5']:
        if tide_model == 'CATS2008':
            # - Define path to CATS2008 data - See parameters.yml file
            grid_file = os.path.join(tide_model_path, 'grid_CATS2008')
            model_file = os.path.join(tide_model_path, 'hf.CATS2008.out')
            epsg_code = tide_model
        else:
            # - Define path to AOTIM5 data - See parameters.yml file
            grid_file = os.path.join(tide_model_path, 'grid_Arc5km2018')
            model_file = os.path.join(tide_model_path, 'h_Arc5km2018')
            epsg_code = 'PSNorth'

        # - pyTDM parameters
        model_format = 'OTIS'
        # - Variable to Read
        var_type = 'z'
        # -
        model_type = None
        model_version = None
        model_scale = None

    else:
        # - Define path to FES2014 data  - See parameters.yml file
        model_directory = os.path.join(tide_model_path,
                                       'fes2014_elevations_and_load',
                                       'fes2014b_elevations',
                                       'ocean_tide')

        model_files = ['2n2.nc', 'eps2.nc', 'j1.nc', 'k1.nc',
                       'k2.nc', 'l2.nc', 'la2.nc', 'm2.nc', 'm3.nc', 'm4.nc',
                       'm6.nc', 'm8.nc', 'mf.nc', 'mks2.nc', 'mm.nc',
                       'mn4.nc', 'ms4.nc', 'msf.nc', 'msqm.nc', 'mtm.nc',
                       'mu2.nc', 'n2.nc', 'n4.nc', 'nu2.nc', 'o1.nc', 'p1.nc',
                       'q1.nc', 'r2.nc', 's1.nc', 's2.nc', 's4.nc', 'sa.nc',
                       'ssa.nc', 't2.nc']
        model_file = [os.path.join(model_directory, x) for x in model_files]
        model_version = 'FES2014'
        model_format = 'FES'
        model_type = 'z'
        model_scale = 1.0 / 100.0
        # - pyTMD parameters
        epsg_code = None
        var_type = None
        grid_file = None

    # -- read tidal constants and interpolate to grid points
    if verbose:
        print('# - Loading Model Parameters.')
    if model_format in ('OTIS', 'ATLAS'):
        amp, ph, d, c = extract_tidal_constants(pt_lon, pt_lat, grid_file,
                                                model_file, epsg_code,
                                                TYPE=var_type,
                                                METHOD='spline',
                                                GRID=model_format)
        # -- calculate complex phase in radians for Euler's
        cph = -1j * ph * np.pi / 180.0
        # -- calculate constituent oscillation
        hc = amp * np.exp(cph)

    else:
        # - model_format = FES
        amp, ph = extract_FES_constants(pt_lon, pt_lat,
                                        model_file, TYPE=model_type,
                                        VERSION=model_version,
                                        METHOD='spline', EXTRAPOLATE=False,
                                        CUTOFF=None,
                                        SCALE=model_scale, GZIP=False)
        # -- available model constituents
        c = ['2n2', 'eps2', 'j1', 'k1', 'k2', 'l2',
             'lambda2', 'm2', 'm3', 'm4', 'm6', 'm8', 'mf', 'mks2', 'mm',
             'mn4', 'ms4', 'msf', 'msqm', 'mtm', 'mu2', 'n2', 'n4', 'nu2',
             'o1', 'p1', 'q1', 'r2', 's1', 's2', 's4', 'sa', 'ssa', 't2']

        # -- calculate complex phase in radians for Euler's
        cph = -1j * ph * np.pi / 180.0
        # -- calculate constituent oscillation
        hc = amp * np.exp(cph)

    # - Other Parameters
    n_sec_x_hour = 60 * 60
    n_sec_x_day = 24 * 60 * 60
    t_date_00 = date1
    t_date_11 = date2

    # - Compute difference in time expressed in hours
    delta_hours = int((t_date_11 - t_date_00).total_seconds() / n_sec_x_hour)

    # - Calculate Number of days relative to Jan 1, 1992 (48622 MJD)
    # - using datetime
    t_jd_ref = datetime.datetime(year=1992, month=1, day=1, hour=0)
    t_est_tide = [t_date_00 + datetime.timedelta(hours=t)
                  for t in range(delta_hours)]

    # - Compute Datetime Values for selected date.
    delta_time = [(t - t_jd_ref).total_seconds() / n_sec_x_day for t in
                  t_est_tide]

    # - Compute Tide Time Series
    tidal_cycle = []
    for h, dt in enumerate(delta_time):
        # -- predict tidal elevations at time and infer minor corrections
        tide_p = predict_tide(dt, hc, c, DELTAT=0,
                              CORRECTIONS=model_format)
        minor_p = infer_minor_corrections(dt, hc, c,
                                          DELTAT=0,
                                          CORRECTIONS=model_format)
        tide_val = tide_p + minor_p
        tidal_cycle.append(tide_val[0])

    return{'tide_ts': tidal_cycle, 'tide_time': t_est_tide}


def main() -> None:
    """
    Main: Compute Hourly Tide Time Series at the selected geographic coordinates
    """
    # - Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Compute Hourly Tide Time Series at the selected
        locations listed inside an ESRI POINT shapefile (.shp).
        """
    )
    # - Positional Arguments
    parser.add_argument('parameters', type=str,
                        help='Processing Parameters File [yml - format].')
    # - Points Coordinates
    parser.add_argument('latitude', type=float,
                        help='Point Latitude [deg].')
    parser.add_argument('longitude', type=float,
                        help='Point Longitude [deg].')
    args = parser.parse_args()

    if not os.path.isfile(args.parameters):
        raise FileNotFoundError('# - Parameters file Not Found.')

    # - Import parameters with PyYaml
    with open(args.parameters, 'r') as stream:
        param_proc = yaml.safe_load(stream)

    # - Create Output Directory
    out_dir = create_dir(param_proc['out_path'], 'tide_ts_pyTMD')

    # - Tide Model
    tide_model = param_proc['model'].upper()
    if tide_model not in ['CATS2008', 'FES2014', 'AOTIM5']:
        raise ValueError(f'# - Unknown Model Selected: {tide_model}')
    else:
        print(f'# - Model Selected: {tide_model}')

    # - Verify Selected Model Data Availability
    tide_model_path = param_proc['path']
    if not os.path.isdir(tide_model_path):
        raise ValueError(f'# - Model Data not available at the provided path:'
                         f' {tide_model_path}')
    else:
        print(f'# - Model Data Path: {tide_model_path}')

    # - Sample Point Coordinates
    pt_lat = args.latitude
    pt_lon = args.longitude
    # - Convert dates to datetime objects
    t_00 = [int(x) for x in param_proc['date1'].split('/')]
    t_11 = [int(x) for x in param_proc['date2'].split('/')]
    t_date_00 = datetime.datetime(year=t_00[2], month=t_00[0], day=t_00[1])
    t_date_11 = datetime.datetime(year=t_11[2], month=t_11[0], day=t_11[1])

    # - Compute Tide Time Series at the selected location
    tide_pt = compute_pt_tide_ts(pt_lon, pt_lat, t_date_00, t_date_11,
                                 tide_model, tide_model_path)
    # - Plot the Computed Daily tidal correction
    tide_ts = tide_pt['tide_ts']
    tide_time = tide_pt['tide_time']

    if param_proc['out_format'].lower() in ['nc', 'nc4', 'netcdf']:
        # - Save Outputs Time Series in NETCDF$ format with xarray
        ds_time_ts = xr.Dataset(data_vars=dict(
            tide_ts=(['time'], tide_ts)),
            coords=dict(time=(['time'], tide_time)),
        )
        # - Dataset Attributes
        ds_time_ts.attrs['model'] = tide_model
        # - Variable Attributes
        ds_time_ts['tide_ts'].attrs['unit'] = 'meters'
        ds_time_ts['tide_ts'].attrs['actual_range']\
            = [np.nanmin(tide_ts), np.nanmax(tide_ts)]
        ds_time_ts['tide_ts'].attrs['_FillValue'] = np.nan
        # - Output File Name
        output_format = 'nc'
        f_name = f'PTide_{tide_model}_Lat{pt_lat}_Lon{pt_lon}_date1_' \
                 f'{t_00[0]:02d}-{t_00[1]:02d}-{t_00[2]}_date2_' \
                 f'{t_11[0]:02d}-{t_11[1]:02d}-{t_11[2]}.{output_format}'
        ds_time_ts.to_netcdf(os.path.join(out_dir, f_name),
                             format="NETCDF4")
    else:
        # - Save tide time series in ascii format.
        # - Output File Name
        output_format = 'txt'
        f_name = f'PTide_{tide_model}_Lat{pt_lat}_Lon{pt_lon}_date1_' \
                 f'{t_00[0]:02d}-{t_00[1]:02d}-{t_00[2]}_date2_' \
                 f'{t_11[0]:02d}-{t_11[1]:02d}-{t_11[2]}.{output_format}'
        # - Save the Computed Daily Correction
        with open(os.path.join(out_dir, f_name),
                  'w', encoding='utf8') as w_fid:
            print('Date'.ljust(25) + 'Tide Height [m]', file=w_fid)
            for cnt, dt in enumerate(list(tide_time)):
                time_str = datetime.datetime.strftime(dt, "%d/%m/%Y %H:%M:%S")
                print(f'{time_str:25}{tide_ts[cnt]}', file=w_fid)


if __name__ == '__main__':
    start_time = datetime.datetime.now()
    main()
    end_time = datetime.datetime.now()
    print(f"# - Computation Time: {end_time - start_time}")
