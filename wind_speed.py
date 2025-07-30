""" 

script to extract data from era5 u and v wind files and 
calculate wind speed

functions in file:

main(): run other functions
get_stats(): gets data from era5 cube
calculate_wind_speed(): calculates wind speed
mask_cube(): masks cube based on shapefile

"""

from datetime import datetime
import iris
import iris.coord_categorisation as ic
import numpy as np
import pandas as pd
import glob
from iris.util import mask_cube_from_shapefile
from cartopy.io.shapereader import Reader


## file path
DATA_DIR = "/data/scratch/lewis.davies/bas/wind"

## file name structure
U_FILE_NAME = 'uwind_*.grib'
V_FILE_NAME = 'vwind_*.grib'

## files
U_WIND_FILES = sorted(glob.glob(f'{DATA_DIR}/{U_FILE_NAME}'))
V_WIND_FILES = sorted(glob.glob(f'{DATA_DIR}/{V_FILE_NAME}'))

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
            csv_file = Path(f"csv_ouputs/wind_speed_{site}_stats.csv")

            ## checks if csv file exists
            if csv_file.is_file():
                
                
                continue


            # ...if it doesn't - produces data
            else:
            
                ## get stats from cubes
                get_stats(index, site, shape_file)

    
    return

## get stats from era5 data
def get_stats(index, site, shape_file):

    ## empty stats
    stats = {"full_date": [], "month": [], "day": [], 
             "hour": [], "grid_square": [], "wind_speed": []}
    
    ## loop through u and v wind files together
    for u_file, v_file in zip(U_WIND_FILES, V_WIND_FILES):
        
        ## load u and v files
        u_cube = iris.load_cube(u_file)
        v_cube = iris.load_cube(v_file)

        ## calculate new wind speed cube from u and v
        ws_cube = calculate_wind_speed(u_cube, v_cube)
        
        ## mask wind speed cube based on site
        ws_site_cube = mask_cube(index, cube, shape_file)
        
        ## loop through times in each cube
        for ws_site_cube in ws_site_cube.slices_over("time"):

            ## define time strings
            time_coord = ws_cube.coord('time')
            full_dt = time_coord.units.num2pydate(time_coord.points)[0]
            month_str = full_dt.strftime("%m")
            day_str = full_dt.strftime("%d")
            hour_str = full_dt.strftime("%H")

            ## loop through flattened data array
            for grid, wind_speed in enumerate(ws_cube.flatten()):

                ## if value is masked, do not add to dict
                if type(data) == np.ma.core.MaskedConstant:
                    
                    pass
                
                else:
                    
                    ## add values to dictionary
                    stats["full_date"].append(full_dt)
                    stats["month"].append(month_str)
                    stats["day"].append(day_str)
                    stats["hour"].append(hour_str)
                    stats["grid_square"].append(grid)
                    stats["wind_speed"].append(wind_speed)

    ## turn dict into dataframe
    site_df = pd.DataFrame(stats)

    ## sort values by date
    site_df = site_df.sort_values(by='full_dt')
    
    ## save dataframe as csv
    site_df.to_csv(f'/home/users/lewis.davies/british_antarctic_survey/analysis_scripts/csv_ouputs/wind_speed_{site}_stats_test.csv')


## calculate wind speed from u and v cubes
def calculate_wind_speed(u_cube, v_cube): 
    
    # multiply u and v components 
    wind_speed_data = np.sqrt(u_cube.data**2 + v_cube.data**2)

    #copy to u_cube to retain other coordinates
    ws_cube = u_cube.copy(wind_speed_data)

    #change cube name
    ws_cube.standard_name = 'wind_speed'

    # convert wind speed units to knots

    ws_cube.convert_units('knots')

    return ws_cube 


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
