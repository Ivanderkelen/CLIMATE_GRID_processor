#!/usr/bin/env python3

"""
EXTRACT CLIMATE_GRID ALTITUDE

Script to grid the pixel altitude (PIXEL_ALT) of CLIMATE_GRID into a static
netcdf file, on the same LAMBERT2008 grid and with the same lat/lon, bounds and
fill value as the files written by preprocess_CLIMATE_GRID.py.
The output can be remapped with remap_CLIMATE_GRID_to_latlon.py
(add "ALTITUDE" to its variable list).

Input: ASCII dump of the CLIMATE_GRID_PIXEL table from the RMI Oracle database
(tab-separated, '#' header lines), containing the columns PIXEL_ID and PIXEL_ALT.

Run from the repository root (grid csv files are read with relative paths).
"""

import os
import numpy as np
import pandas as pd
import xarray as xr
from datetime import date

# User settings
data_dir = "grid_files/"  # holds the pixel table dump; the netcdf is written here too
pixel_table_file = data_dir + "CLIMAT_CLIMATE_GRID_PIXEL_13cols_FULL_DUMP_FIRSTALL--20260623_225142--20260623_225142.txt"
variable = "ALTITUDE"
filename_out = f"{variable}_CLIMATE_GRID_static.nc"

# LAMBERT2008 projection (same as preprocess_CLIMATE_GRID.py)
proj_string = "+proj=lcc +lat_2=50.569898649999999 +lat_1=50.569898649999999 +lon_0=4.553615160000000 +units=m +no_defs +a=6371229.0 +b=6371229.0"

# ---------------------------------------
# Load altitude per pixel
df_pixels = pd.read_csv(pixel_table_file, sep="\t", comment="#")
# drop the '-----' separator line below the column names
df_pixels = df_pixels[pd.to_numeric(df_pixels["PIXEL_ID"], errors="coerce").notna()]
df_pixels = df_pixels[["PIXEL_ID", "PIXEL_ALT"]].astype({"PIXEL_ID": int, "PIXEL_ALT": float})

# Pixel lat/lon and Lambert coordinates: same source as preprocess_CLIMATE_GRID.py,
# so that the output grid (incl. lat/lon) is identical to the meteo files.
df_coords_points = pd.read_csv("grid_5kmx5km.csv", header=1, delimiter=" ")
df_coords_points = df_coords_points.merge(
    df_pixels, on="PIXEL_ID", how="left", validate="one_to_one"
)
assert (
    df_coords_points["PIXEL_ALT"].notna().all()
), "Not all CLIMATE_GRID pixels have an altitude in the pixel table"

# ---------------------------------------
# Build the Lambert grid (same procedure as preprocess_CLIMATE_GRID.py)
df_full_grid = pd.read_csv(
    "lambert_coordinates_full_climate_grid.csv"
)  # Made using R script by Michel Journee
lambert_x_grid_raw = np.sort(df_full_grid["x1"].unique())
lambert_y_grid_raw = np.sort(df_full_grid["x2"].unique())

# Grid spacing (5000 m) and a half-cell tolerance
dx = np.median(np.diff(lambert_x_grid_raw))
dy = np.median(np.diff(lambert_y_grid_raw))
tol_x, tol_y = dx / 2, dy / 2

# Cut this grid to bounding box including grid cells from CLIMATE_GRID
lambert_x_grid = lambert_x_grid_raw[
    (lambert_x_grid_raw >= df_coords_points["LAMBERT_X"].min() - tol_x)
    & (lambert_x_grid_raw <= df_coords_points["LAMBERT_X"].max() + tol_x)
]
lambert_y_grid = lambert_y_grid_raw[
    (lambert_y_grid_raw >= df_coords_points["LAMBERT_Y"].min() - tol_y)
    & (lambert_y_grid_raw <= df_coords_points["LAMBERT_Y"].max() + tol_y)
]

# Nearest grid index for each pixel
xi = np.abs(
    lambert_x_grid[None, :] - df_coords_points["LAMBERT_X"].values[:, None]
).argmin(axis=1)
yi = np.abs(
    lambert_y_grid[None, :] - df_coords_points["LAMBERT_Y"].values[:, None]
).argmin(axis=1)

# Safety checks: fail loudly instead of silently corrupting the grid
assert (
    np.abs(lambert_x_grid[xi] - df_coords_points["LAMBERT_X"]) < tol_x
).all(), "X snapping too far"
assert (
    np.abs(lambert_y_grid[yi] - df_coords_points["LAMBERT_Y"]) < tol_y
).all(), "Y snapping too far"
assert len(set(zip(xi, yi))) == len(df_coords_points), "Pixel collision"

# ---------------------------------------
# Fill 2D arrays
altitude = np.full((len(lambert_y_grid), len(lambert_x_grid)), np.nan)
lat_2d = np.full_like(altitude, np.nan)
lon_2d = np.full_like(altitude, np.nan)

altitude[yi, xi] = df_coords_points["PIXEL_ALT"].values
lat_2d[yi, xi] = df_coords_points["LAT"].values
lon_2d[yi, xi] = df_coords_points["LON"].values

coords = dict(y=lambert_y_grid, x=lambert_x_grid)

ds = xr.Dataset(coords=coords)
ds[variable] = xr.DataArray(
    data=altitude,
    dims=["y", "x"],
    coords=coords,
    attrs=dict(
        long_name="altitude",
        standard_name="surface_altitude",
        description="Mean altitude of the CLIMATE_GRID pixel (PIXEL_ALT in table CLIMATE_GRID_PIXEL)",
        units="m",
        coordinates="lat lon",  # necessary for remapping by cdo
    ),
)

ds["x"].attrs = {
    "units": "E[east]: Easting (meters)",
    "long_name": "x coordinate Lambert Conic Conformal (2SP)",
}
ds["y"].attrs = {
    "units": "N[north]: Northing (meters)",
    "long_name": "y coordinate Lambert Conic Conformal (2SP)",
}

# Add 2D lat and lon, interpolated/extrapolated to also have values outside of Belgium
ds["lat"] = xr.DataArray(
    data=lat_2d,
    dims=["y", "x"],
    coords=coords,
    attrs=dict(
        long_name="latitude",
        description="WGS84 latitude, from values of CLIMATE_GRID, provided per grid point",
        units="degrees_north",
        bounds="lat_bounds",
    ),
).interpolate_na(dim="x", method="linear", fill_value="extrapolate")

ds["lon"] = xr.DataArray(
    data=lon_2d,
    dims=["y", "x"],
    coords=coords,
    attrs=dict(
        long_name="longitude",
        description="WGS84 longitude, from values of CLIMATE_GRID, provided per grid point",
        units="degrees_east",
        bounds="lon_bounds",
    ),
).interpolate_na(dim="x", method="linear", fill_value="extrapolate")


# add lat and lon bounds (necessary for conservative remapping)
def estimate_corners(lat, lon):
    """Cell corners for a 2-D curvilinear grid, CF order SW, SE, NE, NW.

    Each corner is the centroid of the four cell centres meeting at it.
    Boundary centres are extended by linear extrapolation (one ghost ring),
    so that edge cells get their full width/height instead of collapsing
    onto their own centres.
    """

    def _pad(a):
        ny, nx = a.shape
        b = np.empty((ny + 2, nx + 2), dtype=float)
        b[1:-1, 1:-1] = a
        b[0, 1:-1] = 2.0 * a[0, :] - a[1, :]  # ghost row south
        b[-1, 1:-1] = 2.0 * a[-1, :] - a[-2, :]  # ghost row north
        b[:, 0] = 2.0 * b[:, 1] - b[:, 2]  # ghost col west
        b[:, -1] = 2.0 * b[:, -2] - b[:, -3]  # ghost col east
        return b

    def _corners(a):
        b = _pad(a)
        c = 0.25 * (b[:-1, :-1] + b[1:, :-1] + b[:-1, 1:] + b[1:, 1:])
        out = np.empty(a.shape + (4,), dtype=float)
        out[..., 0] = c[:-1, :-1]  # SW
        out[..., 1] = c[:-1, 1:]  # SE
        out[..., 2] = c[1:, 1:]  # NE
        out[..., 3] = c[1:, :-1]  # NW
        return out

    return _corners(np.asarray(lat, float)), _corners(np.asarray(lon, float))


lat_bounds, lon_bounds = estimate_corners(ds["lat"].values, ds["lon"].values)
ds["lat_bounds"] = xr.DataArray(data=lat_bounds, dims=["y", "x", "nv"], coords=coords)
ds["lon_bounds"] = xr.DataArray(data=lon_bounds, dims=["y", "x", "nv"], coords=coords)

# Global attributes
ds.attrs = {
    "creation_date": date.today().strftime("%d-%m-%Y"),
    "creators": "Ghilain N., Van Schaeybroeck B., Vanderkelen I.",
    "contact": "inne.vanderkelen@meteo.be",
    "version": "1.1",
    "affiliation": "Royal Meteorological Institute of Belgium",
    "projection": proj_string,
    "source": os.path.basename(pixel_table_file),
}

# Finite fill value instead of NaN, so CDO remapping does not propagate NaNs
# (see preprocess_CLIMATE_GRID.py)
fillval = -9e33
encoding = {
    name: {"_FillValue": fillval}
    for name in [variable, "lat", "lon", "lat_bounds", "lon_bounds"]
}

if os.path.exists(data_dir + filename_out):
    os.remove(data_dir + filename_out)
ds.to_netcdf(data_dir + filename_out, encoding=encoding)
print(f"Saved as: {data_dir + filename_out}")
