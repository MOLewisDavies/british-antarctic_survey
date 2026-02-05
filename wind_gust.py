"""

script to extract wind gust data from era5 files

functions in file:

main(): runs other functions
process_data(): runs get stats loop
get_stats(): gets data from era5 cube


 """
 
import glob
import warnings
from pathlib import Path

import iris
import numpy as np
import pandas as pd

## import functions module
import functions as func
from constants import BAS_PATH, SHP_FILE, SITES

warnings.filterwarnings("ignore")


## wind gust era5 files
GUST_F_PATH = "/data/scratch/lewis.davies/bas/wind"
GUST_F_NAME = "gust*.grib"
GUST_FILES = sorted(glob.glob(f"{GUST_F_PATH}/{GUST_F_NAME}"))


## runs other functions
def main():
    """ 
    Runs all functions within file.
    """
    
    ## process ERA5 data files
    process_data()
    
    #for combo in COMBOS:
        
        #func.make_heatmap_plots(combo, 'Gust')
        
    #func.make_box_plot('Gust', 'EVERY', 'ALL')
    #func.make_box_plot('Gust', 'EVERY', 'TOP THREE')

    ## make bar plots
    func.make_bar_plots('Gust', 0)
    func.make_bar_plots('Gust', 1)
    
 
## runs get stats functions
def process_data():
    """
    Processes era5 files into stats csv files.
    """

    ## loops through sites
    for index, site in enumerate(SITES):

        ## path for csv files
        csv_file = Path(f"csv_ouputs/Gust_{site}_stats.csv")

        ## checks if csv file exists
        if csv_file.is_file():

            print(f"Gust_{site}_stats.csv exists")

        # ...if it doesn't - produces data
        else:

            ## get stats from cubes
            get_stats(index, site, SHP_FILE)


## get stats from cubes
def get_stats(index, site, shape_file):
    """ 
    Creates dataframe of hourly era5 weather data, separated into grid
    point.
    
    Args:
        
        index (int): index of shape file in list
        site (str): Antarctica site.
        shape_file (path): path to shape file
    """
    
    ## empty dataframe
    big_df = pd.DataFrame()

    for file in GUST_FILES:

        gust_cube = iris.load(file)[0]
        
        gust_cube.units = 'm s-1'
        
        gust_cube.convert_units('knots')
        
        gust_cube = func.convert_timezone(gust_cube)
        site_cube = func.mask_cube(index, gust_cube, shape_file)

        ## Extract coordinates
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
                "Gust": site_cube.data.flatten(),
            }
        )

        # Drop masked values
        df = df.dropna(subset=["Gust"])

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
        f"{BAS_PATH}/csv_ouputs/Gust_{site}_stats.csv"
    )

    
if __name__ == "__main__":
    main()
print("finished")
