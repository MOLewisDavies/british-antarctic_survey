"""

script to analyse precipitaion data from era5

threshold for precipiation used = 0.1mm/hr

functions in file:

main(): runs other functions
process_data(): process era 5 files
get_stats(): gets data from era5 cube
run_all_points_heatmap_func(): creates heatmap for all era5 grid points
make_box_whisker_plots(): makes box plots
make_box_whisker_east_west(): make box plots for east and west of sites

"""

import glob
import warnings
from pathlib import Path

import iris
import numpy as np
import pandas as pd
import seaborn as sns

## import fuctions and constants
import functions as func
from constants import BAS_PATH, SHP_FILE, SITES

warnings.filterwarnings("ignore")

# Seaborn settings - makes it look nice
sns.set_style("darkgrid")

## precipitation files
PRECIP_F_PATH = "/data/scratch/lewis.davies/bas/precip"
PRECIP_F_NAME = "*.grib"
PRECIP_FILES = sorted(glob.glob(f"{PRECIP_F_PATH}/{PRECIP_F_NAME}"))


## runs other functions
def main():

    """ 
    Runs all functions within file.
    """
    
    ## process ERA5 data files
    process_data()
    
    #for combo in COMBOS:
        
        #func.make_heatmap_plots(combo, 'Precip')  
    
    ## make bar plots
    func.make_bar_plots('Precip', 0)
    func.make_bar_plots('Precip', 1)

    return


## run get stats functions
def process_data():
    """
    Processes era5 files into stats csv files.
    """

    ## loops through sites
    for index, site in enumerate(SITES):

        ## path for csv files
        csv_file = Path(f"csv_ouputs/Precip_{site}_stats.csv")

        ## checks if csv file exists
        if csv_file.is_file():

            print(f"Precip_{site}_stats.csv exists")

        # ...if it doesn't - produces data
        else:

            ## get stats from cubes
            get_stats(index, site, SHP_FILE)


## get data from era 5
def get_stats(index, site, shape_file):
    """ 
    Creates dataframe of hourly era5 weather data, separated into grid
    point.
    
    Args:
        
        index (int): index of shape file in list
        site (str): Antarctica site.
        shape_file (path): path to shape file
    """
    
    # Empty dataframe to concatenate to
    big_df = pd.DataFrame()

    ## loop through precipitation files
    for file in PRECIP_FILES:

        ## load file
        precip_cube = iris.load_cube(file)

        precip_cube = func.convert_timezone(precip_cube)
        
        precip_cube.units = 'm'
        
        precip_cube.convert_units('mm')
        
        ## constrain cube to site location
        site_cube = func.mask_cube(index, precip_cube, shape_file)

        # Extract coordinates
        time_points = site_cube.coord("time").points
        time = site_cube.coord("time").units.num2date(time_points)
        lat = site_cube.coord("latitude").points
        lon = site_cube.coord("longitude").points

        # Create full coordinate meshgrid
        time_grid, lat_grid, lon_grid = np.meshgrid(time, lat, lon, indexing="ij")

        # Put data into dataframe
        df = pd.DataFrame(
            {
                "Date and Time": time_grid.flatten(),
                "Latitude": lat_grid.flatten(),
                "Longitude": lon_grid.flatten(),
                "Precip": site_cube.data.flatten(),
            }
        )

        # Drop masked values
        df = df.dropna(subset=["Precip"])

        # Add in month, day and hour columns
        df["Year"] = [d_time.year for d_time in df["Date and Time"]]
        df["Month"] = [d_time.month for d_time in df["Date and Time"]]
        df["Day"] = [d_time.day for d_time in df["Date and Time"]]
        df["Hour"] = [d_time.hour for d_time in df["Date and Time"]]

        # Concatenate to big dataframe
        big_df = pd.concat([big_df, df], ignore_index=True)

    ## sort values by date
    big_df = big_df.sort_values(by="Date and Time")

    # print(big_df)

    ## save dataframe as csv
    big_df.to_csv(
        f"{BAS_PATH}/csv_ouputs/Precip_{site}_stats.csv"
    )


if __name__ == "__main__":
    main()
    
print("finished")
