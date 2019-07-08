#!/usr/bin/env python3.4
"""Convert LPD output to gridded netCDF4.

Store raw counts, with attribute for conversion
"""
from __future__ import division, print_function
import contextlib
import argparse
import datetime
import warnings
import os.path
import math
import sys

import dateutil.relativedelta
import dateutil.rrule
import dateutil.tz
import numpy as np
import scipy.constants
import scipy.io
import netCDF4
import pyproj

UTC = dateutil.tz.tzutc()
HOURS_PER_DAY = 24
MINUTES_PER_HOUR = 60
SECONDS_PER_MINUTE = 60

SECONDS_PER_HOUR = SECONDS_PER_MINUTE * MINUTES_PER_HOUR
SECONDS_PER_DAY = SECONDS_PER_HOUR * HOURS_PER_DAY

UDUNITS_DATE = "%Y-%m-%d %H:%M:%S%z"
ACDD_DATE = "%Y-%m-%dT%H:%M:%S%z"
CALENDAR = "standard"
RUN_DATE = datetime.datetime.now(tz=UTC)
CLOSE_TO_GROUND = 0.05
"""Distance above ground below which particles are counted as
interacting with it (km).
"""

FLUX_WINDOW = 6
"""Length of flux window in hours.

Note
----
This should divide :const:`HOURS_PER_DAY` evenly
"""
OBS_WINDOW = 1
"""Length of observation window in hours.

Note
----
Must divide :const:`FLUX_WINDOW` evenly.
"""

# conversion from #/box to 1e-6 ppm/(g/grid box/hour)

EARTH_GRAVITY = scipy.constants.g
# From Wolfram|Alpha
AIR_DENSITY = 1.2
"""Average air density at 0C.

This value is standard atmosphere value for 500 m, slightly less than
the average elevation of the CONUS. Sea level is 1.275.
"""
THICKNESS = AIR_DENSITY * EARTH_GRAVITY * (CLOSE_TO_GROUND * 1000)
"""Thickness of bottom layer in Pa.

Derived from hydrostatic relation :math:`dP = - \rho g dz`

The standard atmosphere value (P_500 - P_550) is 6 mb:
this is 5.88 mb.
"""
AIR_MOLAR_MASS = (.78084 * 2 * 14.0067 +
                  .20948 * 2 * 15.9994 +
                  # # this bit varys a lot.
                  # ~0% in deserts to ~5%? in tropics
                  # .01 * (2 * 1.00794 + 15.9994) +
                  .00934 * 39.948 +
                  .000380 * (12.0107 + 2*15.9994))
CO2_MOLAR_MASS = (12.0107 + 2*15.9994)
# source: wolfram alpha
CO2_AIR_MASS_RATIO = AIR_MOLAR_MASS / CO2_MOLAR_MASS
"""Ratio of the molar masses of air and CO2.

Note
----
The denominator needs to change if not working with CO2.
(i.e., with CH4, CO, ...).  Remember also to change the
`long_units` attribute on the influence function to reflect
the new species of interest.
"""
MOLES_TO_PPM = EARTH_GRAVITY * (AIR_MOLAR_MASS / 1000) / THICKNESS * 1e6
GRAMS_TO_PPM = CO2_AIR_MASS_RATIO / (1000.*THICKNESS/EARTH_GRAVITY) * 1e6
"""Conversion from flux units to mixing ratio units

Assumes fluxes are in g/m^2/hr;
I think this is independent of the actual area units,
but I don't know
Converts to mixing ratio tendency in ppmv/hr

Notes
-----
.. math::

    F/M_{CO2}/dz = \\Delta n_{CO2}/dt \\\\
    dz = -dP/\\rho g \\\\
    \\Delta X_{CO2} = \\Delta n_{CO2}/n_{air}
                    = \\Delta n_{CO2}/ (\\rho_{air} / M_{air}) \\\\
    \\Delta X_{CO2} = F dt/(M_{CO2} * -dP/(\\rho_{air} g) * M_{air}/\\rho_{air} \\\\
    \\Delta X_{CO2} / dt = F/M_{CO2} / (-dP/g) * M_{air}

Need to convert F to kg if dP uses Pa
X_{CO2} is here in units of 1; multiply by 1e6 to get ppmv
"""
WRF_EARTH_RADIUS = 6.370e6
"""Radius of earth as used by WRF.

meters

Needed to get correct ellipsoid for projection.

Source: Skamarock et al. 2008 pp. 97
"""


def next_larger_multiple(value, multiple):
    """The next multiple of `multiple` >= `value`.

    Parameters
    ----------
    value: float
    multiple: int

    Returns
    -------
    float
        The smallest integer multiple of `multiple`
        greater than or equal to `value`
    """
    return multiple * math.ceil(value / multiple)


def next_smaller_multiple(value, multiple):
    """The largest multiple of `multiple` <= `value`.

    Parameters
    ----------
    value: float
    multiple: int

    Returns
    -------
    float
        The largest integer multiple of `multiple`
        less than or equal to `value`
    """
    return multiple * math.floor(value / multiple)


def read_configuration(savefile_name):
    """Get the configuration for this run.

    For now, only supports idl .sav files.  Since I can't write this
    from python, I'll need to add support for another format at some
    point. Pickle would be the analogue (version 4 shouldn't have any
    problems), but there are other options (Fortran namelist?)

    Parameters
    ----------
    savefile_name: str
        The name of the file with the information for this run.

    Returns
    -------
    dict
        the configuration information.
    """
    return scipy.io.readsav(savefile_name)["input"]


def read_wrf_grid(wrf_name):
    """Get projection information from the given WRF output.

    Parameters
    ----------
    wrf_name: str
        Name of a WRF output file. Must have variables XLAT and XLONG
        as well as and the global attributes from the projection.

    Returns
    -------
    dict
        the mapping infomation:
        "wrf_lat"=tuple(np.ndarray, dict)
            values and attributes for latitude
        "wrf_lon"=tuple(np.ndarray, dict)
            values and attributes for longitude
        "proj_x_coord"=tuple(np.ndarray, dict)
            values and attributes for x
        "proj_y_coord"=tuple(np.ndarray, dict)
            values and attributes for y
        "coord_sys"=dict
            grid mapping attributes
            "var_name" gives grid_mapping
    """
    ds = netCDF4.Dataset(wrf_name.decode("ascii"), "r")
    lats = ds.variables["XLAT"]
    lons = ds.variables["XLONG"]

    dx = ds.DX
    dy = ds.DY
    map_proj = ds.MAP_PROJ

    # 1 is LPDM code for Lambert Conformal Conic projection
    if map_proj == 1:
        # From wrf-python.projection.py l595
        proj = pyproj.Proj(proj="lcc", lat_1=ds.TRUELAT1, lat_2=ds.TRUELAT2,
                           lon_0=ds.STAND_LON, lat_0=ds.MOAD_CEN_LAT,
                           a=WRF_EARTH_RADIUS, b=WRF_EARTH_RADIUS)
    else:
        warnings.warn("No idea how to parse MAP_PROJ_CHAR=" +
                      ds.getncattr("MAP_PROJ_CHAR"))

    if lats.ndim == 3:
        # WRF adds a time dimension
        # We don't use moving nest
        # Not sure how to test for that.
        proj_start = proj(lons[0,0,0], lats[0,0,0])
    else:
        proj_start = proj(lons[0, 0], lats[0, 0])

    grid_start_x = proj_start[0]
    grid_start_y = proj_start[1]

    if map_proj == 1:
        proj_var = dict(
            grid_mapping_name="lambert_conformal_conic",
            # From Skamarock et al. 2008 p. 97
            earth_radius=6.370e6,
            longitude_of_central_meridian=ds.CEN_LON,
            latitude_of_projection_origin=ds.CEN_LAT,
            central_lon=ds.CEN_LON,
            standard_parallel=sorted((ds.TRUELAT1, ds.TRUELAT2), reverse=True),
            false_easting=-grid_start_x,
            false_northing=-grid_start_y,
            var_name="wrf_proj"
            )
    else:
        # already warned
        pass

    proj_x_coord = (np.arange(
            0, (ds.getncattr("WEST-EAST_GRID_DIMENSION") - 1) * dx, dx),
                    dict(standard_name="projection_x_coordinate",
                         units="m", grid_mapping="wrf_proj", axis="X"))
    proj_y_coord = (np.arange(
            0, (ds.getncattr("SOUTH-NORTH_GRID_DIMENSION") - 1) * dy, dy),
                    dict(standard_name="projection_y_coordinate",
                         units="m", grid_mapping="wrf_proj", axis="Y"))

    lat_coord = (lats, dict(standard_name="latitude", units="degrees_north",
                            long_name="flux_latitude",
                            grid_mapping="wrf_proj",
                            description=lons.getncattr("description"),
                            origin=ds.TITLE))
    lon_coord = (lons, dict(standard_name="longitude", units="degrees_east",
                            long_name="flux_longitude",
                            grid_mapping="wrf_proj",
                            description=lons.getncattr("description"),
                            origin=ds.TITLE))

    # this would destroy lats and lons
    # ds.close()
    return dict(wrf_lat=lat_coord,
                wrf_lon=lon_coord,
                proj_x_coord=proj_x_coord,
                proj_y_coord=proj_y_coord,
                # CRS=iris_crs)
                coord_sys=proj_var,
                dx=dx)


def set_global_attributes(ds):
    """Set global attributes required by conventions.

    Conventions are currently CF and ACDD

    Parameters
    ----------
    ds: netcdf4.Dataset
    """
    ds.title = "LPDM CO2 Concentration Footprints"
    ds.summary = ("Gridded CO2 concentration footprints from the output "
                  "of the Lagrangian Particle Dispersion model "
                  "described in Uliasz 1994.")
    ds.Conventions = "CF-1.6 ACDD-1.3"
    ds.history = ("{date:{acdd_format}} {user:s} "
                  "created by {progname:s}").format(
                      date=RUN_DATE, user=os.environ["USER"],
                      acdd_format=ACDD_DATE,
                      progname=sys.argv[0])
    ds.source = ("Gridded outputs from LPDM v?.?.? "
                 "written by Uliasz et al. and modified by Lauvaux")
    ds.standard_name_vocabulary = "CF Standard Name Table v32"
    ds.date_created = "{date:{acdd_format}}".format(
        date=RUN_DATE, acdd_format=ACDD_DATE)
    ds.creator_name = "Daniel Wesloh, Thomas Lauvaux"
    ds.creator_institution = (
        "The Pennsylvania State University "
        "Department of Meteorology and Atmospheric Science")
    ds.date_modified = "{date:{acdd_format}}".format(
        date=RUN_DATE, acdd_format=ACDD_DATE)
    ds.date_metadata_modified = "{date:{acdd_format}}".format(
        date=RUN_DATE, acdd_format=ACDD_DATE)
    ds.product_version = "Py_v0.0.4"
    ds.references = """Uliasz, M. 1994. Lagrangian particle dispersion modeling in mesoscale applications. Environ Model: Comput Methods and Softw for Simulat Environ Pollut and its Adverse Effects (CMP) 2 : 71-."""

    ds.geospatial_vertical_min = 0
    ds.geospatial_vertical_max = CLOSE_TO_GROUND
    ds.geospatial_vertical_positive = "up"
    ds.geospatial_vertical_units = "km AGL"
    # Kind of a cross between Grid and Trajectory
    # Grid covers the first and last two axes;
    # trajectory covers third-to-last
    ds.cdm_data_type = "Grid"

    ds.institution = ds.creator_institution


def create_grid_mapping(ds, wrf_out):
    """Add grid_mapping information from wrf_out to ds.

    Also adds coordinate variables and metadata
    does not set them

    Parameters
    ----------
    ds: netcdf4.Dataset
        The dataset that needs the grid_mapping
    wrf_out: dict
        A dictionary with the grid_mapping information

    Returns
    -------
    str
        The name of the grid_mapping variable
    """
    coord_sys = wrf_out["coord_sys"]
    grid_mapping = coord_sys["var_name"]
    grid_mapping_var = ds.createVariable(grid_mapping, "i1", ())
    grid_mapping_var.setncatts({key: val
                                for key, val in coord_sys.items()
                                if key != "var_name"})

    ydim_var = ds.variables["dim_y"]
    xdim_var = ds.variables["dim_x"]
    ydim_bounds_var = ds.variables["dim_y_bnds"]
    xdim_bounds_var = ds.variables["dim_y_bnds"]

    ydim_var.setncatts(wrf_out["proj_y_coord"][1])
    xdim_var.setncatts(wrf_out["proj_x_coord"][1])
    ydim_var.setncatts(dict(bounds="dim_y_bnds"))
    xdim_var.setncatts(dict(bounds="dim_x_bnds"))

    ydim_bounds_var.setncatts(
        {att: val
         for att, val in wrf_out["proj_y_coord"][1].items()
         if att not in ("standard_name", "axis")})
    xdim_bounds_var.setncatts(
        {att: val
         for att, val in wrf_out["proj_x_coord"][1].items()
         if att not in ("standard_name", "axis")})

    lat_var = ds.createVariable("latitude", "f4", ("dim_y", "dim_x"))
    lon_var = ds.createVariable("longitude", "f4", ("dim_y", "dim_x"))

    # lat_coord = lat
    # lon_coord = lon
    lat_var.setncatts(wrf_out["wrf_lat"][1])
    lon_var.setncatts(wrf_out["wrf_lon"][1])

    return grid_mapping


def set_coord_values(ds, wrf_out, footprint_nbins):
    """Set the coordinate variables from wrf_out.

    Parameters
    ----------
    ds: netcdf4.Dataset
        output dataset
    wrf_out: dict
        wrf data
    footprint_nbins: int
        number of time bins back
    """
    xdim_var = ds.variables["dim_x"]
    ydim_var = ds.variables["dim_y"]
    xdim_bounds_var = ds.variables["dim_x_bnds"]
    ydim_bounds_var = ds.variables["dim_y_bnds"]
    lon_var = ds.variables["longitude"]
    lat_var = ds.variables["latitude"]

    time_back_var = ds.variables["time_before_observation"]
    time_back_bounds_var = ds.variables["time_before_observation_bnds"]

    height_var = ds.variables["height"]
    height_bounds_var = ds.variables["height_bnds"]

    dx = wrf_out["dx"]

    xdim_data = wrf_out["proj_x_coord"][0]
    ydim_data = wrf_out["proj_y_coord"][0]
    xdim_var[:] = xdim_data[:]
    ydim_var[:] = ydim_data[:]

    xdim_bounds_var[:-1,:] = np.column_stack((xdim_data[:-1], xdim_data[1:]))
    xdim_bounds_var[-1,0] = xdim_data[-1]
    xdim_bounds_var[-1,1] = xdim_data[-1] + dx
    ydim_bounds_var[:-1,:] = np.column_stack((ydim_data[:-1], ydim_data[1:]))
    ydim_bounds_var[-1,0] = ydim_data[-1]
    ydim_bounds_var[-1,1] = ydim_data[-1] + dx

    wrf_lats = wrf_out["wrf_lat"][0][0, :, :]
    wrf_lons = wrf_out["wrf_lon"][0][0, :, :]
    lat_var[:, :] = wrf_lats[:, :]
    lon_var[:, :] = wrf_lons[:, :]

    ds.geospatial_lat_min = wrf_lats.min()
    ds.geospatial_lat_max = wrf_lats.max()
    ds.geospatial_lat_units = "degree_north"
    ds.geospatial_lon_min = wrf_lons.min()
    ds.geospatial_lon_max = wrf_lons.max()
    ds.geospatial_lon_units = "degree_east"

    time_back_vals = np.arange(0, footprint_nbins * FLUX_WINDOW, FLUX_WINDOW)
    time_back_var[:] = time_back_vals
    time_back_bounds_var[:-1,:] = np.column_stack((time_back_vals[:-1],
                                                   time_back_vals[1:]))
    time_back_bounds_var[-1,:] = time_back_vals[-2:] + FLUX_WINDOW

    height_var[...] = 0
    height_bounds_var[:] = (0, CLOSE_TO_GROUND)


def strip_array_wrappers(arry):
    """Turn object array of ndarrays into regular ndarray.

    Parameters
    ----------
    arry: np.ndarray

    Returns
    -------
    np.ndarray
        Not an array of ndarrays.
        Might still be object array
    """
    curr = arry
    if curr.ndim == 0:
        if isinstance(curr[...], np.ndarray):
            return strip_array_wrappers(curr[...])
        return curr

    # there is a possibility for infinite looping
    # e.g. [np.ndarray, str, dict] would stay object array
    # impossible if homogeneous (implied by 1-element wrappers)
    while isinstance(curr[0], np.ndarray):
        if curr.shape[0] == 1:
            curr = curr[0]
        else:
            curr = np.array(tuple(curr))

    return curr


def netcdf_compatible_array(arry):
    """Get an array compatible with netCDF dtypes from arry.

    Return an array whose dtype is not object.
    Assumes object arrays contain a single array.

    Parameters
    ----------
    arry: np.ndarray
        The array processed

    Returns
    -------
    np.ndarray
        The same data with a sensible dtype.
    """
    arry = strip_array_wrappers(arry)

    if arry.ndim > 0:
        for _ in range(3):
            if arry.dtype.char != "O" or arry.ndim == 0:
                break

            if arry.shape[0] == 1:
                arry = np.array(arry[0])
            else:
                arry = np.array(tuple(arry))

    if "S" in arry.dtype.char:
        return np.char.decode(arry, "ascii")
    # TODO: ensure no float16, ...
    return arry


def set_up_file(ds, total_sites, footprint_nbins,
                dimy, dimx, wrf_out, time_unit, site_names):
    """Set up dataset for data.

    Set most metadata.
    All attributes and most values.

    Parameters
    ----------
    ds: netcdf4.Dataset
    total_sites: int
    footprint_nbins: int
    dimy: int
    dimx: int
    wrf_out: dict
    time_unit: str
    site_names: np.ndarray
    """
    if site_names.dtype.str[1] not in "SU":
        site_names = site_names.astype("S")
    if site_names.dtype.str[2:] != "1":
        site_names = netCDF4.stringtochar(site_names)
    name_length = site_names.shape[-1]
    name_str_dim = "string{len:d}".format(len=name_length)

    ds.createDimension("observation_time", 0)
    ds.createDimension("site", total_sites)
    ds.createDimension("time_before_observation", footprint_nbins)
    ds.createDimension("dim_y", dimy)
    ds.createDimension("dim_x", dimx)
    ds.createDimension("bnds2", 2)
    ds.createDimension(name_str_dim, name_length)

    obs_time_var = ds.createVariable("observation_time", "f4",
                                     ("observation_time",))
    obs_time_bounds_var = ds.createVariable("observation_time_bnds", "f4",
                                            ("observation_time", "bnds2"))
    time_back_var = ds.createVariable("time_before_observation", "i2",
                                      ("time_before_observation",))
    time_back_bounds_var = ds.createVariable(
        "time_before_observation_bnds", "i2",
        ("time_before_observation", "bnds2"))
    ds.createVariable("dim_y", "f4", ("dim_y",))
    ds.createVariable("dim_y_bnds", "f4", ("dim_y", "bnds2"))
    ds.createVariable("dim_x", "f4", ("dim_x",))
    ds.createVariable("dim_x_bnds", "f4", ("dim_x", "bnds2"))

    grid_mapping = create_grid_mapping(ds, wrf_out)

    site_name_var = ds.createVariable("site_names", "S1",
                                      ("site", name_str_dim))
    # these are roughly 1MB for a three-week lag
    flux_time_var = ds.createVariable(
        "flux_time", "f4", ("observation_time", "time_before_observation"),
        fill_value=-255,
        zlib=True)
    flux_time_bounds_var = ds.createVariable(
        "flux_time_bnds", "f4",
        ("observation_time", "time_before_observation", "bnds2"),
        fill_value=-255,
        zlib=True)

    height_var = ds.createVariable("height", "f4", ())
    height_bounds_var = ds.createVariable("height_bnds", "f4", ("bnds2",))

    infl_fun_var = ds.createVariable(
        # Empirically, the most particles seen in a grid cell is
        # around 5e3.  35*180 particles/obs_time, 9 files/flux_time on 81km grid
        # We have a factor of six wiggle room with i2
        # u2 may be necessary for 3 hourly 243 km fluxes
        # or more particles/obs_time
        "H", "i2",
        ("observation_time", "site", "time_before_observation",
         "dim_y", "dim_x"),
        zlib=True,
        # This will be written and read by flux time, usually,
        # so that chunksize should be 1
        # not sure if chunk should be total_sites or 1 for site dimension
        # total_size gives a chunk as around 5.3 MiB
        # setting this to 1 may help with file size
        # if some towers were not run all the time
        # NUG has default chunk size of 4 MiB
        #   (roughly a disk read on a high-end system)
        chunksizes=(1, total_sites, 1, dimy, dimx),
        # This requires that every cell be written to.
        # This is my intent, and this (as opposed to fill_value=0)
        # will not have troubles with masking most of the domain.
        # Make sure this isn't what's inflating the size
        fill_value=-1,
        )

    lpdm_opts = ds.createVariable("lpdm_configuration", "i1", ())
    lpdm_opts.setncatts({key: netcdf_compatible_array(config[key]).copy()
                         for key in config.dtype.fields.keys()
                         if key.islower()})

    wrf_opts = ds.createVariable("wrf_configuration", "i1", ())
    with contextlib.closing(netCDF4.Dataset(
            config["wrf_file"][0].decode("ascii"))) as wrf_ds:
        wrf_opts.setncatts({att: wrf_ds.getncattr(att)
                            for att in wrf_ds.ncattrs()})

    ########################################################

    obs_time_var.setncatts(dict(long_name="observation_time",
                                # not entirely sure this applies...
                                standard_name="forecast_reference_time",
                                bounds="observation_time_bnds",
                                units=time_unit,
                                calendar=CALENDAR,
                                coverage_content_type="coordinate",
                                # might be a misapplication of CF 9.5
                                cf_role="timeseries_id"))
    obs_time_bounds_var.setncatts(dict(long_name="observation_time_bounds",
                                       units=time_unit,
                                       calendar=CALENDAR))

    time_back_var.setncatts(dict(long_name="time_before_observation",
                                 standard_name="forecast_period",
                                 units="hours",
                                 bounds="time_before_observation_bnds",
                                 coverage_content_type="coordinate",
                                 ))
    time_back_bounds_var.setncatts(dict(
        description="bounds of time_before_observation",
        units="hours"))

    flux_time_var.setncatts(dict(
            long_name="flux_time",
            standard_name="time",
            bounds="flux_time_bnds",
            units=time_unit,
            calendar=CALENDAR,
            coverage_content_type="coordinate",
            ))
    flux_time_bounds_var.setncatts(dict(
            long_name="flux_time",
            units=time_unit,
            calendar=CALENDAR,
            ))

    infl_fun_var.setncatts(dict(
            long_name="influence_function",
            description=("linearisation of the observation operator "
                         "for carbon dioxide mixing ratios at the "
                         "towers in terms of carbon dioxide mass fluxes"),
            units="ppmv/(mol.m^-2.s^-1)",
            long_units="ppmv/(mol_CO2.m^-2.s^-1)",
            coordinates=("flux_time height latitude longitude "
                         "site_names site_heights site_lats site_lons"),
            # I don't think we can justify more than six or so digits
            # of precision.  The transport is too uncertain.
            # The underlying int type doesn't support more than five.
            # The increased locality should also speed up use.
            scale_factor=np.array(CONVERSION_FACTOR, dtype=np.float32),
            grid_mapping=grid_mapping,
            valid_min=np.array(0, dtype=infl_fun_var.dtype),
            # description of coordinate relationships
            cell_methods=(
                # not entirely sure if space and obs time should be in
                # same sum. The two times are another possible
                # combination.
                "height: dim_y: dim_x: sum "
                "observation_time: sum "
                "(interval: {lpdm_timestep:f} seconds) "
                "site: point "
                # this sum is done later than the others
                "flux_time: sum "
                "(interval: {minutes_per_file:d} minutes)"
                "").format(minutes_per_file=(MINUTES_PER_HOUR //
                                             int(config["num_file_per_h"])),
                           lpdm_timestep=float(config["lpdm_timestep"])),
            # What type of thing this is:
            coverage_content_type="modelResult",
            ))
    # I want to store the counts directly
    infl_fun_var.set_auto_maskandscale(False)

    site_lats_var = ds.createVariable("site_lats", "f4", ("site",))
    site_lons_var = ds.createVariable("site_lons", "f4", ("site",))
    site_heights_var = ds.createVariable("site_heights", "f4", ("site",))
    site_lats_var.setncatts(dict(
            units="degrees_north", standard_name="latitude",
            long_name="site_latitude",
            coverage_content_type="coordinate",
            description="latitude of the observation tower site",
            origin="Set in LPD run script"))
    site_lons_var.setncatts(dict(
            units="degrees_east", standard_name="longitude",
            long_name="site_longitude",
            coverage_content_type="coordinate",
            description="longitude of the observation tower site",
            origin="Set in LPD run script"))
    site_name_var.setncatts(dict(
            long_name="name_of_observation_site",
            # most likely an abuse of CF section 9.5
            # cf_role="trajectory_id"
            coverage_content_type="referenceInformation",
            ))
    site_heights_var.setncatts(dict(
            standard_name="height",
            long_name="site_heights",
            description="height of the observation tower intake",
            origin="Set in LPD run script",
            coverage_content_type="coordinate",
            positive="up",
            units="m"))

    height_var.setncatts(dict(
            standard_name="height",
            long_name="flux_influence_height",
            description=("How low the particles have to be "
                         "to be \"influenced\" by the ground"),
            origin="Constant CLOSE_TO_GROUND in carsurf_loop.py",
            coverage_content_type="referenceInformation",
            positive="up",
            units="km", bounds="height_bnds"))
    height_bounds_var.setncatts(dict(
            long_name="height_bounds",
            units="km"))

    # pretty sure this fails somewhat badly at encapsulization
    set_coord_values(ds, wrf_out, footprint_nbins)
    site_name_var[:] = site_names

    return infl_fun_var


def carsurf_loop(config):
    """The main workhorse of the script.

    Loop through the flux times,
    get the particle counts by obs time, site, y, and x
    for those within :const:`CLOSE_TO_GROUND` of the ground
    add those up,
    then write to netCDF4.
    """
    site_names = config["site_names"][0]
    print(site_names)
    total_sites = len(site_names)

    dx = config["dx"].copy()

    # how many hours back should footprints be calculated?
    # roughly, how long should time_before_observation be?
    length = config["length"].copy()
    # I think this is the right way to adjust it
    # We have particle releases at the beginning and end of the `OBS_WINDOW`
    # both of these will need `length` bins back to put observations in
    footprint_nbins = math.ceil((length + OBS_WINDOW) / FLUX_WINDOW)
    # how many days before the first day of the month the simulation goes
    # how far back did LPD calculate trajectories?
    lag = int(config["lag"])

    site_alt = config["alt"][0].copy()
    site_lon = config["lon"][0].copy()
    site_lat = config["lat"][0].copy()

    dimx = int(config["dimx"])
    dimy = int(config["dimy"])

    days_tot = config["num_days"]

    out_dir = config["outdir"][0].decode("ascii")

    year = int(config["year"][0])
    month = int(config["month"][0])

    simulation_earliest_obs = datetime.datetime(year, month, 1)
    # technically the start of the first observation of the next month
    simulation_latest_obs = (simulation_earliest_obs +
                             dateutil.relativedelta.relativedelta(months=+1))
    simulation_zero = (simulation_earliest_obs -
                       datetime.timedelta(days=lag))
    # obs_time_bounds = dateutil.rrule.rrule(
    #     dateutil.rrule.HOURLY, dtstart=simulation_earliest_obs,
    #     interval=OBS_WINDOW, until=simulation_latest_obs,
    #     cache=True)
    n_obs_bins = ((simulation_latest_obs - simulation_earliest_obs) //
                  datetime.timedelta(hours=OBS_WINDOW))

    print("Simulation zero:      ", simulation_zero)
    print("Earliest release time:", simulation_earliest_obs)
    print("Last release time:    ", simulation_latest_obs)

    def obs_var_to_index(sec_since_start):
        """Get the index for the bin.

        Parameters
        ----------
        bin: int

        Returns
        -------
        int
            The index in the NetCDF file created
            0 is the beginning of the simulation,
            at the end of the time window.
        """
        sec_since_first_obs = (sec_since_start -
                               (simulation_earliest_obs -
                                simulation_zero).total_seconds())
        bin_num = int(sec_since_first_obs // (SECONDS_PER_HOUR * OBS_WINDOW))

        # alternate: netCDF4.numtodate(sec_since_start, lpdm_obs_time_unit)
        #            - simulation_unit
        #            // datetime.timedelta(hours=OBS_WINDOW)

        # use time at the end of the window, not the start
        return n_obs_bins - bin_num

    print("Bin index for last release: ",
          obs_var_to_index((simulation_latest_obs -
                            simulation_zero).total_seconds()))
    print("Bin index for first release:",
          obs_var_to_index((simulation_earliest_obs -
                            simulation_zero).total_seconds()))

    # int is more precise than float for this range (up to 4 billion)
    # as we can only have a million particles at a time
    # (for now, run_lprm maxnp)
    # this should also be faster
    # final = np.zeros((total_sites, length, dimy, dimx),
    #                  dtype=np.int32)

    # list of cubes with influence function
    # final_list = collections.deque((), config["lpdm_terase"]//3600)
    # list of release times corresponding to those cubes
    # release_times = collections.deque((), config["lpdm_terase"]//3600)
    # file_name_list = collections.deque((), config["lpdm_terase"]//3600)

    wrf_out = read_wrf_grid(config["wrf_file"][0])

    # LPDM works in minutes for the most part
    time_unit = "minutes since {start:{date_fmt:s}}".format(
        start=simulation_zero, date_fmt=UDUNITS_DATE)

    print("About to create file")
    ds = netCDF4.Dataset(
        os.path.join(
            out_dir,
            "LPDM_{year:04d}_{month:02d}_{flux_window:02d}"
            "hrly_{dx:03d}km_molar_footprints.nc4".format(
                year=year, month=month, flux_window=FLUX_WINDOW,
                dx=int(dx))),
        "w", format="NETCDF4")
    set_global_attributes(ds)

    ds.time_coverage_start = simulation_earliest_obs.strftime(ACDD_DATE)
    ds.time_coverage_end = simulation_latest_obs.strftime(ACDD_DATE)
    ds.time_coverage_duration = "P0000-01-00T00:00:00"
    ds.time_coverage_resolution = "P0000-00-00T{obs_window:02d}:00:00".format(
        obs_window=OBS_WINDOW)

    infl_fun_var = set_up_file(
        ds, total_sites, footprint_nbins,
        dimy, dimx, wrf_out, time_unit, site_names)

    ds.variables["site_lats"][:] = site_lat
    ds.variables["site_lons"][:] = site_lon
    ds.variables["site_heights"][:] = site_alt
    print("Created file")
    # ds.variables["site_names"][:] = np.char.ljust(
    #     site_names, int(site_names.dtype.str[2:]), " ")

    # loop over input files
    # loop goes backward in time from first file output to last
    for step, current_time in zip(
            range(int(HOURS_PER_DAY * days_tot), 0, -FLUX_WINDOW),
            reversed(tuple(dateutil.rrule.rrule(
                    dateutil.rrule.HOURLY,
                    simulation_zero,
                    FLUX_WINDOW,
                    until=simulation_latest_obs)))):
        print("Day: ", step // HOURS_PER_DAY - 1, step/HOURS_PER_DAY,
              "\tHour: ", step % HOURS_PER_DAY)

        # which file to open first (r_{first_file:d}m.dat)
        first_file = step * MINUTES_PER_HOUR

        # end of the period for flux integration
        # earliest file to open? (minutes)
        # now unused.
        # end_flights = first_file - length * MINUTES_PER_HOUR

        print("Current output time:", current_time)

        # set up the cube to receive the data
        # current_time = simulation_zero + datetime.timedelta(minutes=flights)

        # LPDM output codes release time in seconds
        # in a given file, we will have observations from the current time
        # forward for lag days (fluxes influence future obs)

        # last release we care about
        # oldest particles in the first file this iteration
        last_obs = next_larger_multiple(
                (min(current_time + datetime.timedelta(hours=float(length)),
                     simulation_latest_obs) -
                 simulation_zero).total_seconds(),
                OBS_WINDOW * SECONDS_PER_HOUR)
        print("Last release in this iteration:",
              simulation_zero + datetime.timedelta(seconds=last_obs))
        # first release we care about
        # newest particles in the last file this iteration
        first_obs = next_smaller_multiple(
                (max(current_time - datetime.timedelta(hours=FLUX_WINDOW),
                     simulation_earliest_obs) -
                 simulation_zero).total_seconds(),
                OBS_WINDOW * SECONDS_PER_HOUR)
        print("First release in this iteration:",
              simulation_zero + datetime.timedelta(seconds=first_obs))
        print("Last release should be no later than:",
              simulation_zero + datetime.timedelta(seconds=first_obs) +
              datetime.timedelta(hours=FLUX_WINDOW))

        n_obs_bins_here = (last_obs - first_obs) // SECONDS_PER_HOUR // OBS_WINDOW

        # new_cube = create_vars...
        # final_list.append(new_cube)
        # file_name_list.append(
        #     "INFUN_{date:M%m_D%d_H%H}.nc4".format(date=current_time))
        # release_times.append(current_time)

        # go through the files for the hour
        # increase the dtype if LPDM maxnp * n_files_per_hour
        # goes above about 3 billion
        # length should be footprint_nbins
        flux_window_data = np.zeros((dimx, dimy, total_sites, n_obs_bins_here),
                                    dtype=np.int16)
        file_per_hour = int(config["num_file_per_h"])
        minutes_per_file = MINUTES_PER_HOUR // file_per_hour

        flux_time_var = ds.variables["flux_time"]
        flux_time_bounds_var = ds.variables["flux_time_bnds"]

        obs_time_var = ds.variables["observation_time"]
        obs_time_bounds_var = ds.variables["observation_time_bnds"]

        print("Reading LPD output")

        # loop over flux files in this window
        for i in range(FLUX_WINDOW):
            for minute in range(MINUTES_PER_HOUR, 0, -minutes_per_file):
                # get the data from the file
                # Probably in C
                data = np.genfromtxt(
                    os.path.join(
                        config["indir"][0].decode("ascii"),
                        "r_{fli:d}m.dat".format(
                            fli=(first_file - i * MINUTES_PER_HOUR - minute))),
                        # number of lines determined from file
                        skip_header=1,
                        # particle id not needed
                        usecols=(1, 2, 3, 4, 5),
                    )

                # given as x, y, z, site, obs_time?
                # obs_time in seconds, apparently
                mins = (      0,       0,               0,
                              1,     first_obs)
                maxs = (float(dimx*dx), float(dimy*dx), CLOSE_TO_GROUND,
                        total_sites, last_obs)

                # probably in C
                binned_data, bin_desc = np.histogramdd(
                    data, bins=(dimx, dimy, 1, total_sites, n_obs_bins_here),
                    # also kind of cheating
                    range=np.column_stack((mins, maxs))
                )
                del data

                # drop the z dimension from the counts
                # flux_window_data += np.asanyarray(binned_data[:,:,0,:,:],
                #                                   dtype=np.int32)
                # binned_data is a float array,
                # so need unsafe casting to bring back integer counts
                # C
                np.add(flux_window_data, binned_data[:,:,0,:,:],
                       out=flux_window_data, casting="unsafe")
                del binned_data, bin_desc

        print("Read LPD output; writing data")

        # find the indicies where the data should go
        # no data for particles released before first_obs yet
        # problem does not seem to be here, given range semantics
        obs_start = obs_var_to_index(first_obs)
        obs_end = obs_var_to_index(last_obs)
        # print(obs_end, obs_start)

        # simplify the logic and write all times
        # it's a 1-D coord with bounds
        all_dates = tuple(dateutil.rrule.rrule(
            dateutil.rrule.HOURLY,
            simulation_zero + datetime.timedelta(seconds=first_obs),
            OBS_WINDOW,
            until=simulation_zero + datetime.timedelta(seconds=last_obs)))
        # print(all_dates[0], all_dates[-1])
        # print(simulation_zero + datetime.timedelta(seconds=first_obs))
        # print(simulation_zero + datetime.timedelta(seconds=last_obs))

        # observation_time is monotone decreasing by design
        # so the index for the chronologically last time will be
        # lower than that of the chronologically earlier time
        for obs_ind, obs_time_val in zip(
              range(obs_end, obs_start),
              reversed(all_dates)):
            # print("time is", obs_time_val, "mapping to index",
            #       obs_var_to_index((obs_time_val -
            #                         simulation_zero).total_seconds()),
            #       "\nIndex being used:", obs_ind)
            obs_time_var[obs_ind] = netCDF4.date2num(
                obs_time_val,
                time_unit, CALENDAR)
            obs_time_bounds_var[obs_ind, :] = netCDF4.date2num(
                (obs_time_val,
                 obs_time_val - datetime.timedelta(hours=OBS_WINDOW)),
                time_unit, CALENDAR)
        print("Wrote obs times")

        # get loop invariants
        curr_flux_time = netCDF4.date2num(current_time, time_unit, CALENDAR)
        curr_flux_bounds = netCDF4.date2num((
                current_time,
                current_time - datetime.timedelta(hours=FLUX_WINDOW)),
                                            time_unit, CALENDAR)
        print(curr_flux_time, time_unit, "corresponds to", current_time)

        # first obs is at simulation_earliest_obs
        # we are looking at times from current_time
        # to current_time - FLUX_WINDOW
        # first obs is in this window if obs_start == n_obs_bins
        # if obs_start - obs_end is less than footprint_nbins,
        #   need to start writing at time_back=difference
        write_offset = footprint_nbins - ((obs_start - obs_end) *
                                          OBS_WINDOW // FLUX_WINDOW)
        if obs_end == 0:
            # the other time this can occur (the beginning)
            write_offset = 0
        print("Writing data with an offset of", write_offset)

        # now add the data to the file
        # reversing a range is rather annoying.
        for obs_bin_num, obs_ind in enumerate(
              reversed(range(obs_end, obs_start))):
            # np.transpose reverses all dimensions if no spec given
            data_part = np.transpose(flux_window_data[:,:,:,obs_bin_num])

            # final_list[travel_time][:,-travel_time,:,:] = (
            #     CONVERSION_FACTOR * data_part)
            # dataset = netCDF4.Dataset(file_name_list[travel_time], "a")
            # infl_fun = dataset.variables["H"]

            print(infl_fun_var.shape, obs_start, obs_ind, obs_end,
                  data_part.shape)

            # This should support OBS_WINDOW != FLUX_WINDOW, in as
            # much generality as necessary.
            print(obs_start - obs_end, n_obs_bins, obs_bin_num)

            back_bin_num = obs_bin_num * OBS_WINDOW // FLUX_WINDOW
            infl_fun_var[obs_ind, :, back_bin_num+write_offset, :, :] = data_part
            flux_time_var[obs_ind, back_bin_num+write_offset] = curr_flux_time
            flux_time_bounds_var[
                obs_ind, back_bin_num+write_offset, :] = curr_flux_bounds
            del data_part
        ds.sync()
        del flux_window_data
        del curr_flux_time, curr_flux_bounds
        print("Data written")

        # if len(final_list) == final_list.maxlen:
        #     # no more data to be added to the cube
        #     # time to save it and free the memory
        #     finished_cube = final_list.popleft()
        #     release_time = release_times.popleft()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Turn LPDM text outputs into netCDF footprint files"
        )
    parser.add_argument("savefile_name", type=str,
                        help="The file with the run configuration")
    args = parser.parse_args()
    config = read_configuration(args.savefile_name)

    CONVERSION_FACTOR = (
        MOLES_TO_PPM /
        # Normalize by number of particles
        # LPDM releases rel_rate particles / LPDM time step
        # No idea if this needs to be OBS_WINDOW (particles released),
        # FLUX_WINDOW (particles read in), or what
        # TODO: figure that out.
        (
            # How long the particle spent in the box for each count
            (SECONDS_PER_HOUR / config["num_file_per_h"]) *
            # How many particles were released during observation
            OBS_WINDOW *
            config["rel_rate"] *
            (SECONDS_PER_HOUR / config["lpdm_timestep"]))
    )
    carsurf_loop(config)
