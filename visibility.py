"""

script to extract era5 temperature and water content data
and calculate visibility

functions in file:

main(): runs other functions
process_data(): runs get stats loop
get_stats(): gets data from era5 cube
conver_wc_units(): converts water content units
calculate_visibility(): calculate visibility approximation

"""

import glob
import math
import warnings
from pathlib import Path

import iris
import numpy as np
import pandas as pd

## import functions module
import functions as func
from constants import SHP_FILE, SITES, MONTHS, BAS_PATH

warnings.filterwarnings("ignore")

## file path
TEMP_F_PATH = "/data/scratch/lewis.davies/bas/temperature"
VIS_F_PATH = "/data/scratch/lewis.davies/bas/visibility"

LWC_F_NAME = "liq_water_*"
IWC_F_NAME = "ice_water_*"
TEMP_F_NAME = "*.grib"

TEMP_FILES = sorted(glob.glob(f"{TEMP_F_PATH}/{TEMP_F_NAME}"))
LWC_FILES = sorted(glob.glob(f"{VIS_F_PATH}/{LWC_F_NAME}"))
IWC_FILES = sorted(glob.glob(f"{VIS_F_PATH}/{IWC_F_NAME}"))


## run other functions
def main():
    """ 
    Runs all functions within file.
    """
    
    process_data()
    func.make_heatmap_plots('Visibility')

    return


## runs get stats functions
def process_data():
    """
    Processes era5 files into stats csv files.
    """
    
    ## loops through sites
    for index, site in enumerate(SITES):

        ## path for csv files
        csv_file = Path(f"csv_ouputs/Visibility_{site}_stats.csv")

        ## checks if csv file exists
        if csv_file.is_file():

            print(f"Visibility_{site}_stats.csv exists")

        # ...if it doesn't - produces data
        else:

            ## get stats from cubes
            get_stats(index, site, SHP_FILE)


## extract data from era5 files
def get_stats(index, site, shape_file):
    """ 
    Creates dataframe of hourly era5 weather data, separated into grid
    point.
    
    Args:
        
        index (int): index of shape file in list
        site (str): Antarctica site.
        shape_file (path): path to shape file
    """
    
    big_df = pd.DataFrame()

    ## loop through lwc, iwc and temp files together
    for m_name in MONTHS.keys():

        ## load cubes - [0] for files that contain duplicate cubes
        temp_cube = iris.load_cube(f'{TEMP_F_PATH}/temp_{m_name}.grib')[0]
        lwc_cube = iris.load(f'{VIS_F_PATH}/liq_water_{m_name}.grib')[0]
        iwc_cube = iris.load(f'{VIS_F_PATH}/ice_water_{m_name}.grib')[0]

        ## calculate visibility using cubes
        vis_cube = calculate_visibility(temp_cube, lwc_cube, iwc_cube)

        ## mask cube for correct site using shape_file
        site_cube = func.mask_cube(index, vis_cube, shape_file)

        # Extract coordinates
        time_points = site_cube.coord("time").points
        time = site_cube.coord("time").units.num2date(time_points)
        lat = site_cube.coord("latitude").points
        lon = site_cube.coord("longitude").points

        # Create full coordinate meshgrid
        time_grid, lat_grid, lon_grid = np.meshgrid(time, lat, 
                                                    lon, indexing="ij")

        # Put data into dataframe
        df = pd.DataFrame(
            {
                "Date and Time": time_grid.flatten(),
                "Latitude": lat_grid.flatten(),
                "Longitude": lon_grid.flatten(),
                "Visibility": site_cube.data.flatten(),
            }
        )

        # Drop masked values
        df = df.dropna(subset=["Visibility"])

        # Add in month, day and hour columns
        df["Year"] = [d_time.year for d_time in df["Date and Time"]]
        df["Month"] = [d_time.month for d_time in df["Date and Time"]]
        df["Day"] = [d_time.day for d_time in df["Date and Time"]]
        df["Hour"] = [d_time.hour for d_time in df["Date and Time"]]

        # Concatenate to big dataframe
        big_df = pd.concat([big_df, df], ignore_index=True)

    ## sort values by date
    big_df = big_df.sort_values(by="Date and Time")

    ## save dataframe as csv
    big_df.to_csv(
        f"{BAS_PATH}/csv_ouputs/Visibility_{site}_stats.csv"
    )

    return


## convert units
def convert_wc_units(temp_cube, wc):
    """ 
    Converts ice and liquid water content data using temperature data.
    
    Args:
    
        temp_cube (iris cube): iris cube of temperature data.
        wc (iris cube): iris cube of water content data.
        
    Returns:
    
        converted water content data.
    """
    
    ## constant
    R = 287.05

    air_density = 100000 / (R * temp_cube.data)  ## convert pressure units

    converted_data = (wc.data * air_density) * 1000  ## converts kg/kg to gm-3

    ## return converted data
    return converted_data


## calculate visbility
def calculate_visibility(temp_cube, lwc_cube, iwc_cube):
    """ 
    Calculates visibility from temperature & liquid and ice water 
    content.
    
    Args:
    
        temp_cube (iris cube): temperature iris cube.
        lwc_cube (iris cube): liquid water content iris cube.
        iwc_cube (iris cube): ice water content iris cube.
        
    Returns:

        iris cube: visibility iris cube.
    """

    ## liquid water data
    lwc_data = convert_wc_units(temp_cube, lwc_cube)
    
    ## ice water data
    iwc_data = convert_wc_units(temp_cube, iwc_cube)

    ## constants
    N = 50
    e = 0.02

    ## calculate liq, ice, rain and snow coefficients
    liq_ext_coeff = -(math.log(e) / 1.002e3) * ((lwc_data * N) ** 0.6473)
    ice_ext_coeff = (163.0e-3) * iwc_data
    rain_ext_coeff = (5.0e-3) * (lwc_data**0.75)
    snow_ext_coeff = (4.0e-3) * (iwc_data**0.78)

    ## cloud coefficient
    cld_ext_coeff = liq_ext_coeff + ice_ext_coeff + snow_ext_coeff + rain_ext_coeff
    
    ## air coefficient
    air_ext_coeff = (math.log(e)) / (10**5)

    ## total coefficient
    total_ext_coeff = -air_ext_coeff + cld_ext_coeff

    ## equation to calculate visibility from coefficients
    visibility_data = -(math.log(e)) / total_ext_coeff

    ## copy data to cube
    vis_cube = lwc_cube.copy(visibility_data)

    ## change standard name
    vis_cube.standard_name = "visibility_in_air"

    ## assign units
    vis_cube.units = "meters"

    return vis_cube

if __name__ == "__main__":
    main()
print("finished")
