"""

script to extract data from era5 u and v wind files and
calculate wind speed

functions in file:

main(): runs other functions
process_data(): runs get stats loop
get_stats(): gets data from era5 cube
calculate_wind_speed(): calculates wind speed

"""
import glob
import warnings
from pathlib import Path

import iris
import numpy as np
import pandas as pd

## import functions module
import functions as func
from constants import BAS_PATH, SEASONS, SHP_FILE, SITES

warnings.filterwarnings("ignore")


## file path
DATA_DIR = "/data/scratch/lewis.davies/bas/wind"

## file name structure
U_FILE_NAME = "uwind_*.grib"
V_FILE_NAME = "vwind_*.grib"

## runs other functions
def main():
    """ 
    Runs all functions within file.
    """
    
    ## process data from era5 files
    for season in SEASONS:
        
        process_data(season)
        
    #for combo in COMBOS:
        
        #func.make_heatmap_plots(combo, 'Wind_Speed')

    ## make bar plots
    func.make_bar_plots('Wind_Speed', 0)
    func.make_bar_plots('Wind_Speed', 1)
    func.make_bar_plots('Wind_Speed', 2)
    
    ## make seasonal plots
    func.seasonal_plots('Wind_Speed', 'summer')
    func.seasonal_plots('Wind_Speed', 'winter')


## run get stats functions
def process_data(season):
    """
    Processes era5 files into stats csv files.
    
    Args:
    
    season: Season to analyse
    """
    
    ## loops through sites
    for index, site in enumerate(SITES):

        ## path for csv files
        csv_file = Path(f"csv_ouputs/Wind_Speed_{site}_stats_{season}.csv")

        ## checks if csv file exists
        if csv_file.is_file():

            print(f"Wind_Speed_{site}_stats_{season}.csv exists")

        # ...if it doesn't - produces data
        else:

            ## get stats from cubes
            get_stats(index, site, SHP_FILE, season)


## get stats from era5 data
def get_stats(index, site, shape_file, season):
    """ 
    Creates dataframe of hourly era5 weather data, separated into grid
    point.
    
    Args:
        
        index (int): index of shape file in list
        site (str): Antarctica site.
        shape_file (path): path to shape file
    """
    
    if season == 'flying':
        
        ## files
        U_WIND_FILES = sorted(glob.glob(f"{DATA_DIR}/{U_FILE_NAME}"))
        V_WIND_FILES = sorted(glob.glob(f"{DATA_DIR}/{V_FILE_NAME}"))
        
    else:
        
        U_WIND_FILES = sorted(glob.glob(f"{DATA_DIR}/{season}/{U_FILE_NAME}"))
        V_WIND_FILES = sorted(glob.glob(f"{DATA_DIR}/{season}/{V_FILE_NAME}"))
    
    big_df = pd.DataFrame()

    ## loop through u and v wind files together
    for u_file, v_file in zip(U_WIND_FILES, V_WIND_FILES):

        ## load u and v files
        u_cube = iris.load_cube(u_file)
        v_cube = iris.load_cube(v_file)

        ## calculate new wind speed cube from u and v
        ws_cube = calculate_wind_speed(u_cube, v_cube)

        ws_cube.convert_units('knots')
        
        ws_cube = func.convert_timezone(ws_cube)
        
        ## mask wind speed cube based on site
        ws_site_cube = func.mask_cube(index, ws_cube, shape_file)

        # Extract coordinates
        time_points = ws_site_cube.coord("time").points
        time = ws_site_cube.coord("time").units.num2date(time_points)
        lat = ws_site_cube.coord("latitude").points
        lon = ws_site_cube.coord("longitude").points

        # Create full coordinate meshgrid
        time_grid, lat_grid, lon_grid = np.meshgrid(time, lat, lon, indexing="ij")

        # Put data into dataframe
        df = pd.DataFrame(
            {
                "Date and Time": time_grid.flatten(),
                "Latitude": lat_grid.flatten(),
                "Longitude": lon_grid.flatten(),
                "Wind_Speed": ws_site_cube.data.flatten(),
            }
        )

        # Drop masked values
        df = df.dropna(subset=["Wind_Speed"])

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
        f"{BAS_PATH}/csv_ouputs/Wind_Speed_{site}_stats_{season}.csv"
    )


## calculate wind speed from u and v cubes
def calculate_wind_speed(u_cube, v_cube):
    """ 
    Calculates wind speed from u and v wind component data.
    
    Args:

        u_cube (iris cube): u component of wind iris cube.
        v_cube (iris cube): v component of wind iris cube.
        
    Returns:
    
        iris cube: wind speed iris cube.
    """
    # multiply u and v components
    wind_speed_data = np.sqrt(u_cube.data**2 + v_cube.data**2)

    # copy to u_cube to retain other coordinates
    ws_cube = u_cube.copy(wind_speed_data)

    # change cube name
    ws_cube.standard_name = "wind_speed"

    return ws_cube


if __name__ == "__main__":
    main()
print("finished")
