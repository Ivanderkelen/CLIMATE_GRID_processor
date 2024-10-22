# CLIMATE GRID preprocessor

Convert the downloaded CLIMATE_GRID dataset from RMI into netcdf data format
The Climate grid dataset can be downloaded here: https://opendata.meteo.be/download (free for academics). 

This repository contains a python script to convert the dataset from .csv to NetCDF format with the necessary attributes.

(c) Inne Vanderkelen, Bert van Schaeybroeck, Nicolas Ghilain
October 2024

## 0. Requirements

Installed version of python (eg through a conda installation, possible with [miniconda](https://docs.anaconda.com/miniconda/) or [Anaconda Navigator](https://www.anaconda.com/download) )

The environment with packages used by the script is included in the [environment.yml](./environent.yml) file and can be installed using conda as follows: 

```
# Create a conda environment based on environment.yml
conda env create -f environment.yml

# Activate the environment
conda activate env_climategrid
```

## 1. Running the script

Run the script **[convert_climate_grid_to_netcdf.py](./convert_climate_grid_to_netcdf.py)** as follows: 


```
python convert_climate_grid_to_netcdf.py
```

