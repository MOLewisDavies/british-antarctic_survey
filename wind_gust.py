""" 

script to extract wind gust data from era5 files

functions in file:

main(): runs other functions 
get_stats(): gets data from era5 cube
mask_cube(): masks cube based on shapefile

"""

import iris
import glob
from pathlib import Path
import seaborn as sns
import numpy as np
import pandas as pd

## wind gust era5 files
gust_f_path = "/data/scratch/lewis.davies/bas/wind"
gust_f_name = "gust*.grib"
gust_files =  sorted(glob.glob(f"{gust_f_path}/{gust_f_name}"))

## site names
SITES = ["theron_hills", "provider", "argentina_range"]

## shapefile path
filepath = "/home/users/lewis.davies/british_antarctic_survey/East Station Met AOI"
shape_name = '*.shp'
SHAPES = glob.glob(f'{filepath}/{shape_name}')

## runs other functions
def main():
    
    ## get shapefile
    for shape_file in SHAPES:
        
        ## loops through sites
        for index, site in enumerate(SITES):

            ## path for csv files
            csv_file = Path(f"csv_ouputs/wind_gust_{site}_stats.csv")

            ## checks if csv file exists
            if csv_file.is_file():
              
                continue

            # ...if it doesn't - produces data
            else:
            
                ## get stats from cubes
                get_stats(index, site, shape_file)


def get_stats(index, site, shape_file):
    
    ## empty stats
    stats = {"full_date": [], "month": [], "day": [], 
             "hour": [], "grid_square": [], "wind_gust": []}
    
    
    for file in gust_files:
        
        gust_cube = iris.load_cube(file)
        
        site_cube = mask_cube(index, gust_cube, shape_file)
        
        for site_cube in site_cube.slices_over('time'):
            
            ## get time values
            time_coord = site_cube.coord('time')
            full_dt = time_coord.units.num2pydate(time_coord.points)[0]
            month_str = full_dt.strftime("%m")
            day_str = full_dt.strftime("%d")
            hour_str = full_dt.strftime("%H")
            
            ## loop through flattened data array 
            for grid, gust in enumerate(site_cube.data.flatten()):
                
                ## if value is masked, do not add to dict
                if type(gust) == np.ma.core.MaskedConstant:
                    
                    pass
                
                else:
                    
                    ## add values to dictionary
                    print(full_dt)
                    stats["full_date"].append(full_dt)
                    stats["month"].append(month_str)
                    stats["day"].append(day_str)
                    stats["hour"].append(hour_str)
                    stats["grid_square"].append(grid)
                    stats["wind_gust"].append(gust)
    
    
    ## turn dict into dataframe   
    site_df = pd.DataFrame(stats)
    
    ## sort values by date
    site_df = site_df.sort_values(by='full_dt')

    ## save dataframe as csv
    site_df.to_csv(f'/home/users/lewis.davies/british_antarctic_survey/analysis_scripts/csv_ouputs/wind_gust_{site}_stats.csv')

## mask cube using shapefile polygons
def mask_cube(index, cube, shape_file):
    
    ## define shapefile reader
    shape_reader = Reader(shape_file)
    
    ## get list of polygons from file
    shape_con = shape_reader.geometries()
    
    ## list polygons 
    shape_cons = list(shape_con)
    
    ## get polygon that matches site region
    shape_con = shape_cons[index]
    
    ## mask cube based on shapefile polygon
    shape_cube = mask_cube_from_shapefile(cube, shape_con)
    
    ## return masked cube
    return shape_cube


if __name__ == "__main__":
    main()
print('finished')
