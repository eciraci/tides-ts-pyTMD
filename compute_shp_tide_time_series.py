#!/usr/bin/env python
u"""
compute_shp_tide_time_series.py
Written by Enrico Ciraci' (05/2022)

Compute Hourly Tide Time Series at the selected locations listed inside
a ESRI shapefile (.shp).

Shapely Geometry Type Accepted:
- LineString.

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
    fiona: Fiona reads and writes geographic data files.
           https://fiona.readthedocs.io
    shapely: Manipulation and analysis of geometric objects in the Cartesian
           plane.
           https://shapely.readthedocs.io/en/stable
    geopandas: GeoPandas is an open source project to add support for
           geographic data to pandas objects.
           https://geopandas.org/en/stable/
    PyYAML: YAML framework for the Python programming language.
           https://pyyaml.org/
UPDATE HISTORY:
06/09/2022 - --stats, -S: option added: Compute Basic Statistics for each of the
          output time series.
"""
# - Python Dependencies
from __future__ import print_function
import os
import argparse
import datetime
import yaml
import numpy as np
import xarray as xr
import fiona
import geopandas as gpd
from multiprocessing import Pool
from functools import partial
# - Utility functions
from utils.create_dir import create_dir
from compute_pt_tide_time_series import compute_pt_tide_ts


def main() -> None:
    """
    Main: Compute Hourly Tide Time Series at the selected locations
          listed inside a ESRI shapefile (.shp).
    """
    # - Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Compute Hourly Tide Time Series at the selected
        locations listed inside an ESRI shapefile (.shp).
        """
    )
    # - Positional Arguments
    parser.add_argument('parameters', type=str,
                        help='Processing Parameters File [yml - format].')
    # - Compute Basic Statistics for each of the output time series
    parser.add_argument('--stats', action='store_true',
                        help='Compute Basic Statistics for each of the output'
                             ' time series')
    args = parser.parse_args()

    if not os.path.isfile(args.parameters):
        raise FileNotFoundError('# - Parameters file Not Found.')

    # - Import parameters with PyYaml
    with open(args.parameters, 'r', encoding='utf8') as stream:
        param_proc = yaml.safe_load(stream)

    # - Tide Model Code
    tide_model = param_proc['model'].upper()

    # - Create Output Directory
    out_dir = create_dir(param_proc['out_path'], 'tide_ts_pyTMD')
    # - Extract mask file name
    mask_name = param_proc['shp'].split('/')[-1][:-4]
    out_dir = create_dir(out_dir, mask_name)

    # - Convert dates to datetime objects
    t_00 = [int(x) for x in param_proc['date1'].split('/')]
    t_11 = [int(x) for x in param_proc['date2'].split('/')]
    t_date_00 = datetime.datetime(year=t_00[2], month=t_00[0], day=t_00[1])
    t_date_11 = datetime.datetime(year=t_11[2], month=t_11[0], day=t_11[1])

    # - Add sub-period output directory
    out_dir = create_dir(out_dir, f'{t_00[1]:02d}-{t_00[0]:02d}-{t_00[2]:04d}_'
                                  f'{t_11[1]:02d}-{t_11[0]:02d}-{t_11[2]:04d}')

    # - Import Sample Locations Coordinates
    s_pt_df = gpd.read_file(param_proc['shp'])
    # - Extract Coordinates Reference System - EPSG Code -> String
    ref_crs = s_pt_df.crs.to_string()

    # - GeoPandas DataFrame number of records
    nr = 0      # - Consider only points included in the first data record
    if s_pt_df.geometry[nr].geom_type == 'LineString':
        # - extract vertexes coordinates
        pt_lon, pt_lat = s_pt_df.geometry[nr].coords.xy
        pt_coords = list(zip(pt_lon, pt_lat))
    else:
        raise TypeError(f'Unsupported geometry type: '
                        f'{s_pt_df.geometry[nr].geom_type} ')

    with Pool(param_proc['nproc']) as p:
        kwargs = {'date1': t_date_00, 'date2': t_date_11,
                  'tide_model': tide_model,
                  'tide_model_path': param_proc['path']}
        map_func = partial(compute_pt_tide_ts, **kwargs)
        pts_tide_ts_list = p.map(map_func, pt_coords)

    # - Save the obtained time series
    for cnt, pt_tide in enumerate(pts_tide_ts_list):
        tide_ts = pt_tide['tide_ts']
        tide_time = pt_tide['tide_time']
        pt_lon, pt_lat = pt_coords[cnt]
        if param_proc['out_format'].lower() in ['nc', 'nc4', 'netcdf']:
            # - Save Outputs Time Series in NETCDF format using xarray
            ds_time_ts = xr.Dataset(data_vars=dict(
                tide_ts=(['time'], tide_ts)),
                coords=dict(time=(['time'], tide_time)),
            )
            # - Dataset Attributes
            ds_time_ts.attrs['model'] = tide_model
            # - Variable Attributes
            ds_time_ts['tide_ts'].attrs['unit'] = 'meters'
            ds_time_ts['tide_ts'].attrs['actual_range'] \
                = [np.nanmin(tide_ts), np.nanmax(tide_ts)]
            ds_time_ts['tide_ts'].attrs['_FillValue'] = np.nan
            # - Output File Name
            output_format = 'nc'
            f_name = f'PT{cnt+1:03d}_{tide_model}' \
                     f'_Lat{pt_lat}_Lon{pt_lon}_date1_' \
                     f'{t_00[0]:02d}-{t_00[1]:02d}-{t_00[2]}_date2_' \
                     f'{t_11[0]:02d}-{t_11[1]:02d}-{t_11[2]}' \
                     f'.{output_format}'
            ds_time_ts.to_netcdf(os.path.join(out_dir, f_name),
                                 format='NETCDF4')
        else:
            # - Save tide time series in ascii format.
            # - Output File Name
            output_format = 'txt'
            f_name = f'PT{cnt+1:03d}_{tide_model}' \
                     f'_Lat{pt_lat}_Lon{pt_lon}_date1_' \
                     f'{t_00[0]:02d}-{t_00[1]:02d}-{t_00[2]}_date2_' \
                     f'{t_11[0]:02d}-{t_11[1]:02d}-{t_11[2]}' \
                     f'.{output_format}'
            # - Save the Computed Time Series
            with open(os.path.join(out_dir, f_name),
                      'w', encoding='utf8') as w_fid:
                print('Date'.ljust(25) + 'Tide Height [m]', file=w_fid)
                for rnt, dt in enumerate(list(tide_time)):
                    time_str = datetime.datetime\
                        .strftime(dt, '%d/%m/%Y %H:%M:%S')
                    print(f'{time_str:25}{tide_ts[rnt]}', file=w_fid)

        if args.stats:
            f_name_st = f'PT{cnt + 1:03d}_{tide_model}' \
                        f'_Lat{pt_lat}_Lon{pt_lon}_date1_' \
                        f'{t_00[0]:02d}-{t_00[1]:02d}-{t_00[2]}_date2_' \
                        f'{t_11[0]:02d}-{t_11[1]:02d}-{t_11[2]}' \
                        f'_STATS.txt'
            # - Save the Computed Time Series
            with open(os.path.join(out_dir, f_name_st),
                      'w', encoding='utf8') as s_fid:
                print(f'PT{cnt + 1:03d}_{tide_model}', file=s_fid)
                print(f'Latitude: {pt_lat} Longitude: {pt_lon}', file=s_fid)
                print(f'Analyzed Period: ', file=s_fid)
                print(f'Date1 {t_00[0]:02d}-{t_00[1]:02d}-{t_00[2]}',
                      file=s_fid)
                print(f'Date2 {t_11[0]:02d}-{t_11[1]:02d}-{t_11[2]}',
                      file=s_fid)
                print('Point Coordinates:', file=s_fid)
                print(f'Latitude: {pt_lat} Longitude: {pt_lon}', file=s_fid)
                print(f'Maximum Annual Value [m]: {np.max(tide_ts)}',
                      file=s_fid)
                print(f'Minimum Annual Value [m]: {np.min(tide_ts)}',
                      file=s_fid)
                print(f'Annual Standard Deviation [m]: {np.std(tide_ts)}',
                      file=s_fid)

    # - Save Sample Point Locations inside a Point Shapefile
    smp_point_mask = os.path.join(out_dir, 'spt_coords_mask.shp')
    # - Define Shapefile Mask Schema
    schema = {
        'geometry': 'Point',
        'properties': [('Name', 'str'), ('id', 'int'),
                       ('Lat', 'float'), ('Lon', 'float')]
    }
    with fiona.open(smp_point_mask, mode='w', driver='ESRI Shapefile',
                    schema=schema, crs=ref_crs) as poly_shp:
        for cnt, pt_tide in enumerate(pts_tide_ts_list):
            pt_lon, pt_lat = pt_coords[cnt]
            # - SP Point
            sp_point = {'type': 'Point', 'coordinates': [pt_lon, pt_lat]}
            # -
            row_dict = {
                # - Geometry [Point]
                'geometry': sp_point,
                # - Properties [based on the schema defined above]
                'properties': {'Name': f'PT{cnt+1:03d}',
                               'id': cnt,
                               'Lat': pt_lat,
                               'Lon': pt_lon
                               },
            }
            poly_shp.write(row_dict)


if __name__ == '__main__':
    start_time = datetime.datetime.now()
    main()
    end_time = datetime.datetime.now()
    print(f'# - Computation Time: {end_time - start_time}')
