# Scripts to Preprocess CLIMATE_GRID Dataset

This repository contains a collection of Python scripts to:
1. Download the CLIMATE_GRID dataset from the RMI Oracle database as vector data and store it as CSV.
2. Convert the dataset to NetCDF format with the necessary attributes.
3. Remap the dataset to a self-defined lat-lon grid.

(c) Inne Vanderkelen, Bert van Schaeybroeck, Nicolas Ghilain
June 2024

## Overview of Files

### Scripts
- **`preprocess_CLIMATE_GRID.py`**:
  - Connects to the RMI Oracle database and downloads the CLIMATE_GRID data.
  - Converts and saves the data as a NetCDF file.
- **`extract_altitude_CLIMATE_GRID.py`**:
  - Reads the pixel altitude (`PIXEL_ALT`) from an ASCII dump of the `CLIMATE_GRID_PIXEL` table.
  - Saves it as a static NetCDF (`ALTITUDE_CLIMATE_GRID_static.nc`) on the same grid as the output of `preprocess_CLIMATE_GRID.py`.
- **`remap_CLIMATE_GRID_to_latlon.py`**:
  - Remaps the projected (raw) CLIMATE_GRID NetCDF data to a lat/lon grid.
  - The target grid can be user-defined within the script.

### Data Files
- **`grid_5kmx5km.csv`**:
  - Contains pixel lat and lon variables to transpose to a user-defined grid.
  - Used by `preprocess_CLIMATE_GRID.py`.
- **`CLIMATE_GRID_meta.csv`**:
  - Contains the variable names and units for the CLIMATE_GRID dataset.
  - Used to write the NetCDF file by `preprocess_CLIMATE_GRID.py`.
- **`lambert_coordinates_full_climate_grid.csv`**:
  - Lists pixel IDs and corresponding x, y coordinates in LAMBERT2008.
  - Used by `preprocess_CLIMATE_GRID.py`.

## Usage

0. **Requirements**

Installed version of python (eg through a conda installation, possible with [miniconda](https://docs.anaconda.com/miniconda/) or [Anaconda Navigator](https://www.anaconda.com/download) )

The environment with packages used by the script is included in the [environment.yml](./environent.yml) file and can be installed using conda as follows: 

```
# Create a conda environment based on environment.yml
conda env create -f environment.yml

# Activate the environment
conda activate env_climategrid
```

1. **Download and Preprocess the Data:**
   ```bash
   python preprocess_CLIMATE_GRID.py
    ```

2. **(Optional) Extract the altitude:**
Set the path to the `CLIMATE_GRID_PIXEL` table dump in [extract_altitude_CLIMATE_GRID.py](./extract_altitude_CLIMATE_GRID.py), then run:
   ```bash
   python extract_altitude_CLIMATE_GRID.py
    ```
   To remap it as well, add `"ALTITUDE"` to the variable list of the remap script.

3. **Remapping the netcdf files**
Open [remap_CLIMATE_GRID_to_latlon.py](./remap_CLIMATE_GRID_to_latlon.py) and do user adjustments (grid, metadata etc). Then execute the script: 

   ```bash
   python remap_CLIMATE_GRID_to_latlon.py
    ```