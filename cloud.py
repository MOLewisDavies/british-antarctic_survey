"""

script to extract era5 cloud cover data

functions in file:

main(): runs other functions
process_data(): runs get stats loop
get_stats(): gets data from era5 cube
get_cloud_base(): calculates lowest cloud base


"""

import glob
import warnings
from pathlib import Path

import iris
import numpy as np
import pandas as pd
import seaborn as sns

## import functions module
import functions as func
from constants import SHP_FILE, SITES, BAS_PATH

warnings.filterwarnings("ignore")

# Seaborn settings - makes it look nice
sns.set_style("darkgrid")

## file path
CLOUD_F_PATH = "/data/scratch/lewis.davies/bas/cloud"
CLOUD_F_NAME = "*.grib"
CLOUD_FILES = sorted(glob.glob(f"{CLOUD_F_PATH}/{CLOUD_F_NAME}"))


## runs other functions
def main():
    """ 
    Runs all functions within file.
    """
    
    process_data()
    func.make_heatmap_plots('cloud_base')


## processes files into csv files
def process_data():
    """
    Processes era5 files into stats csv files.
    """
    
    ## loops through sites
    for index, site in enumerate(SITES):

        ## path for csv files
        csv_file = Path(f"csv_ouputs/cloud_base_{site}_stats.csv")

        ## checks if csv file exists
        if csv_file.is_file():

            print(f"csv_ouputs/cloud_base_{site}_stats.csv", 'exists')
            
            continue

        ## ...if it doesn't - produces data
        else:

            ## get data from era5
            get_stats(index, site, SHP_FILE)

    return


##  get data from era 5
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

    ## loops through cloud files
    for file in CLOUD_FILES:

        ## load cube
        cloud_cube = iris.load(file)[0]

        ## mask cube to polygon based on site
        site_cube = func.mask_cube(index, cloud_cube, shape_file)

        # Extract coordinates
        time_points = site_cube.coord("time").points
        time = site_cube.coord("time").units.num2date(time_points)
        lat = site_cube.coord("latitude").points
        lon = site_cube.coord("longitude").points

        # Create full coordinate meshgrid
        time_grid, lat_grid, lon_grid = np.meshgrid(time, lat, lon, indexing="ij")

        ## spressure level constraints
        con_950 = iris.Constraint(pressure=950)
        con_975 = iris.Constraint(pressure=975)
        con_1000 = iris.Constraint(pressure=1000)

        ## Separate cubes by pressure levels
        cube_950 = site_cube.extract(con_950)
        cube_975 = site_cube.extract(con_975)
        cube_1000 = site_cube.extract(con_1000)

        ## get cube data
        cube_950_data = cube_950.data.flatten()
        cube_975_data = cube_975.data.flatten()
        cube_1000_data = cube_1000.data.flatten()

        # Put data into dataframe
        df = pd.DataFrame(
            {
                "Date and Time": time_grid.flatten(),
                "Latitude": lat_grid.flatten(),
                "Longitude": lon_grid.flatten(),
                "Cloud_950": cube_950_data,
                "Cloud_975": cube_975_data,
                "Cloud_1000": cube_1000_data,
            }
        )

        # Drop masked values
        df = df.dropna(subset=["Cloud_950"])
        df = df.dropna(subset=["Cloud_975"])
        df = df.dropna(subset=["Cloud_1000"])

        # Add in month, day and hour columns
        df["Year"] = [d_time.year for d_time in df["Date and Time"]]
        df["Month"] = [d_time.month for d_time in df["Date and Time"]]
        df["Day"] = [d_time.day for d_time in df["Date and Time"]]
        df["Hour"] = [d_time.hour for d_time in df["Date and Time"]]

        # Concatenate to big dataframe
        big_df = pd.concat([big_df, df], ignore_index=True)

    heights = get_cloud_base(big_df)

    print(heights)
    print(big_df)

    big_df["cloud_base"] = heights

    ## sort values by date
    big_df = big_df.sort_values(by="Date and Time")

    ## save dataframe as csv
    big_df.to_csv(
        f"{BAS_PATH}/csv_ouputs/cloud_base_{site}_stats.csv"
    )


## finds lowest cloud base height
def get_cloud_base(big_df):
    """ 
    Calculates cloud base height based on okta cloud fraction.
    
    Args:
    
        big_df (DataFrame): cloud fractions on pressure levels.
        
    Returns:
    
        list: cloud base height in feet.
    """
    
    
    heights = []

    for index, row in big_df.iterrows():

        clouds = {
            row["Cloud_1000"]: [1000, 350],
            row["Cloud_975"]: [975, 1100],
            row["Cloud_950"]: [950, 1850],
        }

        ## loop through clouds dict
        for key, height in clouds.items():

            ## return height if cloud above 3 oktas or on last height to check
            if key > 3 / 8 or height[0] == 950:

                heights.append(height[1])

                break

            ## continue to next height/value pair
            else:

                continue

    return heights



if __name__ == "__main__":
    main()

print("finished")
