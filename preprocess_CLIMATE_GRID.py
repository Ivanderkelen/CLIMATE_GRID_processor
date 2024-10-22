#!/usr/bin/env python3

"""
PREPROCESS CLIMATE_GRID

Script to download and grid CLIMATE_GRID data into netcdf files
! This needs to be executed on kili

Steps:
1. Connect to RMI Oracle DB and download raw data into .csv files (only works on kili)
2. Grid the raw .csv files into netcdf files for further use.

Author: I. Vanderkelen, June 2024
"""

# Load necessary modules
import warnings

warnings.filterwarnings("ignore")

import os
import pandas as pd
import numpy as np
import getpass
import oracledb
import xarray as xr
from datetime import date
import cartopy.crs as ccrs

# import rioxarray


# User option: download raw oracle data. If False, load intermediate .csv files
download_from_oracle = False  # if True, has to be executed on kili

# User settings
data_dir = "/mnt/HDS_CLIMATE/CLIMATE/CLIMATE_GRID/"
dataset = "CLIMATE_GRID"
variables = [
    "EVAPOTRANS_REF",
    "SUN_INT",
    "SUN_DURATION",
    "PRECIP_DURATION",
    "WIND_PEAK_SPEED",
    "PRECIP_1H_MAX",
    "EVAPOTRANS_REF",
    "TEMP_MAX",
    "HUMIDITY_RELATIVE",
    "TEMP_AVG",
    "WIND_SPEED",
    "PRESSURE",
    "SHORT_WAVE_FROM_SKY",
    "SUN_INT_HORIZ",
    "PRECIP_QUANTITY",
    "TEMP_MIN",
]

# Initial and end year for the data request
init_yr = 1950
end_yr = 2023

if download_from_oracle:
    # Oracle database connection - user interference only once
    username = getpass.getpass("Enter oracle username: ")
    password = getpass.getpass("Enter password: ")

# Process data for each variable
for variable in variables:
    # Define filenames for intermediate files
    filename_csv = f"climate_atlas_{variable}_{dataset}_{init_yr}_{end_yr}.csv"
    filename_municipalities_csv = (
        f"climate_atlas_{variable}_{dataset}_municipalities_{init_yr}_{end_yr}.csv"
    )

    if download_from_oracle and not os.path.isfile(data_dir + filename_csv):

        host = "delphi.oma.be"
        service_name = "rmidbs1.oma.be"
        port = 1521

        params = oracledb.ConnectParams(host=host, port=port, service_name=service_name)
        connection = oracledb.connect(user=username, password=password, params=params)

        init_date = f"{init_yr}0101"
        end_date = f"{end_yr + 1}0101"

        # Get info from grid
        df_grid = pd.read_sql(
            "SELECT PIXEL_ID, PIXEL_LON_CENTER, PIXEL_LAT_CENTER FROM CLIMATE_GRID_PIXEL ORDER BY PIXEL_ID",
            connection,
        )
        pixel_ids = df_grid["PIXEL_ID"].tolist()

        # Get info from municipalities
        df_municipality_ids = pd.read_sql(
            "SELECT NAME, CODE_INS from MUNICIPALITY", connection
        )
        municipalities_codes = df_municipality_ids["CODE_INS"].tolist()

        print(f"Oracle climate grid query for variable: {variable}")

        req_str = f"""
            SELECT CODE_PIXEL_ID_CODE_INS, DATE_END, {variable}
            FROM {dataset}
            WHERE DATE_END = DATE_BEGIN
            AND DATE_BEGIN >= TO_DATE('{init_date}', 'YYYYMMDD')
            AND DATE_END < TO_DATE('{end_date}', 'YYYYMMDD')
            ORDER BY DATE_END, CODE_PIXEL_ID_CODE_INS
        """

        df_data = pd.read_sql(req_str, connection)
        connection.close()

        df_data.columns = ["location", "time", "value"]
        df_data["time"] = pd.to_datetime(df_data["time"]).dt.date

        # Separate grid cells and municipalities
        df_gridcells = df_data[df_data["location"].isin(pixel_ids)]
        df_municipality = df_data[df_data["location"].isin(municipalities_codes)]

        # Save intermediate data
        df_gridcells.to_csv(data_dir + filename_csv)
        df_municipality.to_csv(data_dir + filename_municipalities_csv)

    else:
        print(f"Loading csv file: {data_dir + filename_csv}")
        df_gridcells = pd.read_csv(data_dir + filename_csv)

    # Prepare for gridding into netcdf
    df_pivotted = df_gridcells.pivot_table(
        index="location", columns="time", values="value", fill_value=np.nan
    )
    data = df_pivotted.values
    dates = pd.to_datetime(df_pivotted.columns)

    # Load climate grid metadata on variables and units
    meta = pd.read_csv("CLIMATE_GRID_meta.csv", delimiter=";")

    # OLD, wrong PROJ string corresponding to LAMBERT1972
    # proj_string = "+proj=lcc +lat_1=49.83333388888889 +lat_2=51.16666722222222 +lat_0=90 +lon_0=4.367486666666666 +x_0=150000.013 +y_0=5400088.438 +ellps=intl +units=m +no_defs"

    # LAMBERT2008 projection
    proj_string = "+proj=lcc +lat_2=50.569898649999999 +lat_1=50.569898649999999 +lon_0=4.553615160000000 +units=m +no_defs +a=6371229.0 +b=6371229.0"

    # Load the pixel lat and lon variable and use this to transpose to own defined grid
    df_coords_points = pd.read_csv(
        "grid_5kmx5km.csv", header=1, delimiter=" "
    )  # Lat/lon and Lambert coordinates for all pixels in CLIMATE_DATA

    # Load the full grid, created based on the proj_string and following bounding points
    df_full_grid = pd.read_csv(
        "lambert_coordinates_full_climate_grid.csv"
    )  # Made using R script by Michel Journee
    lambert_x_grid_raw = df_full_grid["x1"].unique()
    lambert_y_grid_raw = df_full_grid["x2"].unique()

    # Cut this grid to bounding box including grid cells from CLIMATE_GRID
    lambert_x_grid = lambert_x_grid_raw[
        (lambert_x_grid_raw >= df_coords_points["LAMBERT_X"].min())
        & (lambert_x_grid_raw <= df_coords_points["LAMBERT_X"].max())
    ]
    lambert_y_grid = lambert_y_grid_raw[
        (lambert_y_grid_raw >= df_coords_points["LAMBERT_Y"].min())
        & (lambert_y_grid_raw <= df_coords_points["LAMBERT_Y"].max())
    ]

    # Find the nearest index in the lons and lats grids and add this to coordinates dataframe
    def find_nearest(array, value):
        idx = (np.abs(array - value)).argmin()
        return idx

    df_coords_points["LAMBERT_X_INDEX"] = df_coords_points["LAMBERT_X"].apply(
        lambda x: find_nearest(lambert_x_grid, x)
    )
    df_coords_points["LAMBERT_Y_INDEX"] = df_coords_points["LAMBERT_Y"].apply(
        lambda x: find_nearest(lambert_y_grid, x)
    )

    print(f"Converting to netcdf: {variable}")

    # Create empty arrays to fill with gridded data and lat/lon
    grid_data = np.full((len(dates), len(lambert_y_grid), len(lambert_x_grid)), np.nan)
    lat_2d = np.full((len(lambert_y_grid), len(lambert_x_grid)), np.nan)
    lon_2d = np.full((len(lambert_y_grid), len(lambert_x_grid)), np.nan)

    # Fill the grid data array
    for _, row in df_coords_points.iterrows():
        lambert_x_idx = int(row["LAMBERT_Y_INDEX"])
        lambert_y_idx = int(row["LAMBERT_X_INDEX"])
        pixel_id = int(row["PIXEL_ID"])

        grid_data[:, lambert_x_idx, lambert_y_idx] = data[pixel_id - 1, :]
        lat_2d[lambert_x_idx, lambert_y_idx] = df_coords_points.loc[
            df_coords_points["PIXEL_ID"] == pixel_id, "LAT"
        ].values[0]
        lon_2d[lambert_x_idx, lambert_y_idx] = df_coords_points.loc[
            df_coords_points["PIXEL_ID"] == pixel_id, "LON"
        ].values[0]

    # Get metadata from meta dataframe
    unit = meta.loc[meta["variable"] == variable, "unit"].values[0]
    long_name = meta.loc[meta["variable"] == variable, "long_name"].values[0]
    description = meta.loc[meta["variable"] == variable, "description"].values[0]

    # this is necessary for remapping by cdo
    coordinates = "lat lon"

    # Create data array
    da = xr.DataArray(
        data=grid_data,
        dims=["time", "y", "x"],
        coords=dict(
            y=lambert_y_grid,
            x=lambert_x_grid,
            time=dates,
        ),
        attrs=dict(
            long_name=long_name,
            description=description,
            units=unit,
            coordinates=coordinates,
        ),
    )

    da["x"].attrs = {
        "units": "E[east]: Easting (meters)",
        "long_name": "x coordinate Lambert Conic Conformal (2SP)",
    }
    da["y"].attrs = {
        "units": "N[north]: Northing (meters)",
        "long_name": "y coordinate Lambert Conic Conformal (2SP)",
    }

    # Convert to dataset and give dataset attributes
    ds = da.to_dataset(name=variable)

    # Add 2D lat and lon in lambert coordinates
    ds["lat"] = xr.DataArray(
        data=lat_2d,
        dims=["y", "x"],
        coords=dict(
            y=lambert_y_grid,
            x=lambert_x_grid,
        ),
        attrs=dict(
            long_name="latitude",
            description="WGS84 latitude, from values of CLIMATE_GRID, provided per grid point",
            units="degrees_north",
            bounds="lat_bounds",
        ),
    )

    # Interpolate to also have lat values outside of Belgium
    ds["lat"] = ds["lat"].interpolate_na(
        dim="x", method="linear", fill_value="extrapolate"
    )

    ds["lon"] = xr.DataArray(
        data=lon_2d,
        dims=["y", "x"],
        coords=dict(
            y=lambert_y_grid,
            x=lambert_x_grid,
        ),
        attrs=dict(
            long_name="longitude",
            description="WGS84 longitude, from values of CLIMATE_GRID, provided per grid point",
            units="degrees_east",
            bounds="lon_bounds",
        ),
    )

    # Interpolate to also have lon values outside of Belgium
    ds["lon"] = ds["lon"].interpolate_na(
        dim="x", method="linear", fill_value="extrapolate"
    )

    # Pass CRS using rioxarray - don't do this because it inhibits regridding using CDO.
    # ds.rio.write_crs(ccrs.Projection(proj_string), inplace=True)

    # add lat and lon bounds (necessary for conservative remapping)

    def estimate_corners(lat, lon):
        """
        Estimate the corners of each grid cell for 2D curvilinear grids.
        The corners are estimated by averaging the midpoints of neighboring cells.
        Args:
            lat (np.array): 2D array of latitudes for cell midpoints.
            lon (np.array): 2D array of longitudes for cell midpoints.
        Returns:
            lat_corners, lon_corners: 2D arrays of shape (nlat, nlon, 4) containing the estimated
                                    latitude and longitude of each corner for each grid cell.
        """
        nlat, nlon = lat.shape
        lat_corners = np.zeros((nlat, nlon, 4))  # Initialize arrays for corners
        lon_corners = np.zeros((nlat, nlon, 4))

        # Internal cells
        for i in range(0, nlat):
            for j in range(0, nlon):
                im1 = np.max([0, i - 1])
                ip1 = np.min([nlat - 1, i + 1])
                jm1 = np.max([0, j - 1])
                jp1 = np.min([nlon - 1, j + 1])

                # Average to get the SW corner of the (i, j) cell
                lat_corners[i, j, 0] = (
                    lat[i, j] + lat[i, jm1] + lat[im1, j] + lat[im1, jm1]
                ) / 4
                lon_corners[i, j, 0] = (
                    lon[i, j] + lon[i, jm1] + lon[im1, j] + lon[im1, jm1]
                ) / 4

                # Average to get the NW corner of the (i, j) cell
                lat_corners[i, j, 1] = (
                    lat[i, j] + lat[i, jp1] + lat[im1, j] + lat[im1, jp1]
                ) / 4
                lon_corners[i, j, 1] = (
                    lon[i, j] + lon[i, jp1] + lon[im1, j] + lon[im1, jp1]
                ) / 4

                # Average to get the NE corner of the (i, j) cell
                lat_corners[i, j, 2] = (
                    lat[i, j] + lat[i, jp1] + lat[ip1, j] + lat[ip1, jp1]
                ) / 4
                lon_corners[i, j, 2] = (
                    lon[i, j] + lon[i, jp1] + lon[ip1, j] + lon[ip1, jp1]
                ) / 4

                # Average to get the SE corner of the (i, j) cell
                lat_corners[i, j, 3] = (
                    lat[i, j] + lat[i, jm1] + lat[ip1, j] + lat[ip1, jm1]
                ) / 4
                lon_corners[i, j, 3] = (
                    lon[i, j] + lon[i, jm1] + lon[ip1, j] + lon[ip1, jm1]
                ) / 4

        return lat_corners, lon_corners

    lat_bounds, lon_bounds = estimate_corners(ds["lat"].values, ds["lon"].values)

    ds["lat_bounds"] = xr.DataArray(
        data=lat_bounds,
        dims=["y", "x", "nv"],
        coords=dict(y=lambert_y_grid, x=lambert_x_grid),
    )

    ds["lon_bounds"] = xr.DataArray(
        data=lon_bounds,
        dims=["y", "x", "nv"],
        coords=dict(y=lambert_y_grid, x=lambert_x_grid),
    )

    # Add global attributes
    ds.attrs = {
        "creation_date": date.today().strftime("%d-%m-%Y"),
        "creators": "Ghilain N., Van Schaeybroeck B., Vanderkelen I.",
        "contact": "inne.vanderkelen@meteo.be",
        "version": "1.1",
        "affiliation": "Royal Meteorological Institute of Belgium",
        "projection": proj_string,
    }

    filename_out = (
        f"{variable}_CLIMATE_GRID_{dates.year.min()}_{dates.year.max()}_daily.nc"
    )

    # Export to netcdf
    ds.to_netcdf(data_dir + filename_out, encoding={"time": {"dtype": "int32"}})
    print(f"Saved as: {data_dir + filename_out}")
