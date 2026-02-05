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
from datetime import datetime
from pathlib import Path

import iris
import numpy as np
import pandas as pd
import seaborn as sns

## import functions module
import functions as func
from constants import BAS_PATH, SHP_FILE, SITES

warnings.filterwarnings("ignore")

# Seaborn settings - makes it look nice
sns.set_style("darkgrid")

## file path
CLOUD_F_PATH = "/data/scratch/lewis.davies/bas/cloud_m_level"
CLOUD_F_NAME = "*.grib"
CLOUD_FILES = sorted(glob.glob(f"{CLOUD_F_PATH}/{CLOUD_F_NAME}"))


## runs other functions
def main():
    """ 
    Runs all functions within file.
    """
    
    ## process data from era5 files
    process_data()
    
    #for combo in COMBOS:
        
        #func.make_heatmap_plots(combo, 'cloud_base')
        
    ## make bar plots
    func.make_bar_plots('cloud_base', 0)
    func.make_bar_plots('cloud_base', 1)
    

## processes files into csv files
def process_data():
    """
    Processes era5 files into stats csv files.
    """
    
    ## loops through sites
    for index, site in enumerate(SITES):

        ## path for csv files
        csv_file = Path(f"csv_ouputs/cloud_base_{site}_stats_ml.csv")

        ## checks if csv file exists
        if csv_file.is_file():

            print(f"csv_ouputs/cloud_base_{site}_stats_ml.csv", 'exists')
            
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
    
    ## create big dataframe for all stats
    big_df = pd.DataFrame()

    ## loops through cloud files
    for file in CLOUD_FILES:
        
        t1 = datetime.now()
        ## load cube
        cloud_cube = iris.load(file)[0]
        
        ## convert timezone of cube
        func.convert_timezone(cloud_cube)
        
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
        con_123 = iris.Constraint(model_level_number=123)
        con_124 = iris.Constraint(model_level_number=124)
        con_125 = iris.Constraint(model_level_number=125)
        con_126 = iris.Constraint(model_level_number=126)
        con_127 = iris.Constraint(model_level_number=127)
        con_128 = iris.Constraint(model_level_number=128)
        con_129 = iris.Constraint(model_level_number=129)
        con_130 = iris.Constraint(model_level_number=130)
        con_131 = iris.Constraint(model_level_number=131)
        con_132 = iris.Constraint(model_level_number=132)
        con_133 = iris.Constraint(model_level_number=133)
        con_134 = iris.Constraint(model_level_number=134)
        con_135 = iris.Constraint(model_level_number=135)
        con_136 = iris.Constraint(model_level_number=136)
        con_137 = iris.Constraint(model_level_number=137)

        ## Separate cubes by pressure levels
        cube_123 = site_cube.extract(con_123)    
        cube_124 = site_cube.extract(con_124)
        cube_125 = site_cube.extract(con_125)
        cube_126 = site_cube.extract(con_126)
        cube_127 = site_cube.extract(con_127)
        cube_128 = site_cube.extract(con_128)
        cube_129 = site_cube.extract(con_129)
        cube_130 = site_cube.extract(con_130)
        cube_131 = site_cube.extract(con_131)
        cube_132 = site_cube.extract(con_132)
        cube_133 = site_cube.extract(con_133)
        cube_134 = site_cube.extract(con_134)
        cube_135 = site_cube.extract(con_135)
        cube_136 = site_cube.extract(con_136)
        cube_137 = site_cube.extract(con_137)

        # Put data into dataframe
        df = pd.DataFrame(
            {
                "Date and Time": time_grid.flatten(),
                "Latitude": lat_grid.flatten(),
                "Longitude": lon_grid.flatten(),
                "Cloud_123": cube_123.data.flatten(),
                "Cloud_124": cube_124.data.flatten(),
                "Cloud_125": cube_125.data.flatten(),
                "Cloud_126": cube_126.data.flatten(),
                "Cloud_127": cube_127.data.flatten(),
                "Cloud_128": cube_128.data.flatten(),
                "Cloud_129": cube_129.data.flatten(),
                "Cloud_130": cube_130.data.flatten(),
                "Cloud_131": cube_131.data.flatten(),
                "Cloud_132": cube_132.data.flatten(),
                "Cloud_133": cube_133.data.flatten(),
                "Cloud_134": cube_134.data.flatten(),
                "Cloud_135": cube_135.data.flatten(),
                "Cloud_136": cube_136.data.flatten(),
                "Cloud_137": cube_137.data.flatten(),
            }
        )

        # Drop masked values
        df = df.dropna(subset=["Cloud_123"])
        df = df.dropna(subset=["Cloud_124"])
        df = df.dropna(subset=["Cloud_125"])
        df = df.dropna(subset=["Cloud_126"])
        df = df.dropna(subset=["Cloud_127"])
        df = df.dropna(subset=["Cloud_128"])
        df = df.dropna(subset=["Cloud_129"])
        df = df.dropna(subset=["Cloud_130"])
        df = df.dropna(subset=["Cloud_131"])
        df = df.dropna(subset=["Cloud_132"])
        df = df.dropna(subset=["Cloud_133"])
        df = df.dropna(subset=["Cloud_134"])
        df = df.dropna(subset=["Cloud_135"])
        df = df.dropna(subset=["Cloud_136"])
        df = df.dropna(subset=["Cloud_137"])

        # Add in month, day and hour columns
        df["Year"] = [d_time.year for d_time in df["Date and Time"]]
        df["Month"] = [d_time.month for d_time in df["Date and Time"]]
        df["Day"] = [d_time.day for d_time in df["Date and Time"]]
        df["Hour"] = [d_time.hour for d_time in df["Date and Time"]]

        # Concatenate to big dataframe
        big_df = pd.concat([big_df, df], ignore_index=True)
        
        t2 = datetime.now()
        
        print('single loop: ', t2-t1)
        print(big_df)

    heights = get_cloud_base(big_df)

    print(heights)
    print(big_df)

    big_df["cloud_base"] = heights

    ## sort values by date
    big_df = big_df.sort_values(by="Date and Time")

    ## save dataframe as csv
    big_df.to_csv(
        f"{BAS_PATH}/csv_ouputs/cloud_base_{site}_stats_ml.csv"
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

    ## loop through big dataframe
    for index, row in big_df.iterrows():

        ## model levels and cloud heights to check
        clouds = {
            row["Cloud_137"]: [137, 32],
            row["Cloud_136"]: [136, 101],
            row["Cloud_135"]: [135, 176],
            row["Cloud_134"]: [134, 259],
            row["Cloud_133"]: [133, 349],
            row["Cloud_132"]: [132, 448],
            row["Cloud_131"]: [131, 556],
            row["Cloud_130"]: [130, 673],
            row["Cloud_129"]: [129, 800],
            row["Cloud_128"]: [128, 940],
            row["Cloud_127"]: [127, 1095],
            row["Cloud_126"]: [126, 1263],
            row["Cloud_125"]: [125, 1443],
            row["Cloud_124"]: [124, 1640],
            row["Cloud_123"]: [123, 1856],
        }

        ## loop through clouds dict
        for key, height in clouds.items():

            ## return height if cloud above 3 oktas or on last height to check
            if key > 3 / 8 or height[0] == 123:

                heights.append(height[1])

                break

            ## continue to next height/value pair
            else:

                continue

    return heights



if __name__ == "__main__":
    main()

print("finished")
