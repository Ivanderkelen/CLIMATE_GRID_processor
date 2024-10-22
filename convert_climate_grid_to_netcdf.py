#!/usr/bin/env python3

"""
PREPROCESS CLIMATE_GRID

Script to grid CLIMATE_GRID data into netcdf files

Author: I. Vanderkelen, October 2024
"""

# Load necessary modules
import os
import pandas as pd
import numpy as np
import xarray as xr
from datetime import date
import cartopy.crs as ccrs


##########################
### 0. Initialisation

creator = "Inne Vanderkelen"
contact = "inne.vanderkelen@meteo.be"
version = "1.1"

# User settings
path = "./climateGrid/"
outpath = "./climateGrid_netcdf/"
file_metadata = 'climategrid_pixel_metadata.csv'



##########################
### 1. Reading in data

df_pixel_metadata = pd.read_csv(os.path.join(path, file_metadata), delimiter=";")

# read all files into one dataframe
files = [file for file in os.listdir(path) if (file.endswith('.csv') and not file_metadata in file)]

dfs = []
for file in files:
    df = pd.read_csv(os.path.join(path, file), delimiter=";")
    dfs.append(df)

df = pd.concat(dfs, ignore_index=True)

df_latlon = pd.merge(df, df_pixel_metadata, left_on="pixel_id", right_on="PIXEL_ID")



##########################
### 2. Read in variable meta data

# Initialize an empty dictionary to store the variable metadata
vars_dict = {}

# Open and read the file
with open(os.path.join(path,'climategrid_parameters_description.txt'), 'r') as file:
    for line in file:
        # Split the line into three parts: key, unit, and description
        parts = line.split('\t')
        
        # The first part is the key, the second part is the unit, and the rest is the description
        variable = parts[0].strip()
        unit = parts[2].strip()
        description = parts[3].strip()
        
        # Combine the unit and description into a tuple or a string if you prefer
        vars_dict[variable] = {'unit' : unit, 'long_name' : description }
        
# manual correction
vars_dict['SUN_DURATION']['unit'] = 'hours/day'
vars_dict['SHORT_WAVE_FROM_SKY']['unit'] = 'kWh/m2/day'



##########################
### 4. Convert variable per variable to netcdf


# loop over the different variables 
for variable in vars_dict.keys(): 
    print(f"Converting to netcdf: {variable}")

    # subselect the variable, lat, lon and time
    df_variable = df_latlon[['day', variable.lower(), 'PIXEL_LAT_CENTER', 'PIXEL_LON_CENTER']]

    # only keep the days for which there are obs in the variable 
    df_variable = df_variable.dropna(subset=[variable.lower()])

    # create coordinate arrays
    lats = df_variable['PIXEL_LAT_CENTER'].sort_values().unique()
    lons = df_variable['PIXEL_LON_CENTER'].sort_values().unique()


    # to determine the dates, only include 

    dates = df_variable['day'].unique()
    years = np.unique([date.split('/')[0] for date in dates])


    # create empty list for data arrays
    all_dataarrays = []

    # process year per year
    for year in [years[0]]:
        print("processing year " + year)
        # Process chunk of data and save intermediate results to file
        dates_by_year = pd.to_datetime([date for date in dates if date.startswith(year)])

        # select the corresponding year in the dataframe 
        df_variable_year = df_variable[df_variable['day'].str.startswith(year+'/')]

        grid_data = np.full((len(dates_by_year), len(lats), len(lons)), np.nan)

        # Fill the grid data array
        for _, row in df_variable_year.iterrows():
            time_idx = np.where(dates_by_year == pd.to_datetime(row['day']))[0][0]
            lat_idx = np.where(lats == row['PIXEL_LAT_CENTER'])[0][0]
            lon_idx = np.where(lons == row['PIXEL_LON_CENTER'])[0][0]
            
            grid_data[time_idx, lat_idx, lon_idx] = row[variable.lower()]

        # Get metadata from variable dictionary
        unit = vars_dict[variable]["unit"]
        long_name = vars_dict[variable]["long_name"]

        # this is necessary for remapping by cdo

        # Create data array
        da_year = xr.DataArray(
            data=grid_data,
            dims=["time", "lat", "lon"],
            coords=dict(
                lat=lats,
                lon=lons,
                time=dates_by_year,
            ),
            attrs=dict(
                long_name=long_name,
                units=unit,
            ),
        )

        # add data array for year to the list
        all_dataarrays.append(da_year)


    # merge data arrays for every year along the time dimension
    da = xr.concat(all_dataarrays, dim='time')


    da["lat"].attrs = {
        "units": "degrees_north",
        "long_name": "WGS84 latitude, from values of CLIMATE_GRID, provided per grid point",
    }
    da["lon"].attrs = {
        "units": "degrees_east",
        "long_name": "WGS84 longitude, from values of CLIMATE_GRID, provided per grid point",
    }


    # Convert to dataset and give dataset attributes
    ds = da.to_dataset(name=variable)



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
        dims=["lat", "lon", "nv"],
        coords=dict(lat=lats, lon=lons),
    )


    ds["lon_bounds"] = xr.DataArray(
        data=lon_bounds,
        dims=["lat", "lon", "nv"],
        coords=dict(lat=lats, lon=lons),
    )


    # Add global attributes
    ds.attrs = {
        "creation_date": date.today().strftime("%d-%m-%Y"),
        "creators": creator,
        "contact": contact,
        "version": version,
        "affiliation": "Royal Meteorological Institute of Belgium"
    }

    filename_out = (
        f"{variable}_CLIMATE_GRID_{dates.year.min()}_{dates.year.max()}_daily.nc"
    )


    # if output directory doesn't exist yet, create it. 
    if not os.path.exists(outpath):
        os.makedirs(outpath)


    # Export to netcdf
    ds.to_netcdf(outpath + filename_out, encoding={"time": {"dtype": "int32"}})
    print(f"Saved as: {outpath + filename_out}")

