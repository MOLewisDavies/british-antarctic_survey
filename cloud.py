""" 

script to extract era5 cloud cover data

functions in file:

main(): runs other functions 
get_stats(): gets data from era5 cube
get_cloud_base(): calculates lowest cloud base
get_values(): gets non-maksed values
mask_cube(): masks cube based on shapefile

"""

import iris
import glob
from pathlib import Path
import seaborn as sns
import numpy as np
import pandas as pd
from iris.util import mask_cube_from_shapefile
from cartopy.io.shapereader import Reader

## file path
CLOUD_F_PATH = "data/scratch/lewis.davies/bas/cloud"
CLOUD_F_NAME = "cloud_base_*"
CLOUD_FILES = sorted(glob.glob(f"{CLOUD_F_PATH}/{CLOUD_F_NAME}"))

## pressure levels (hpa) and heights (ft)
LEVELS = {1000: 350, 975: 1100, 950: 1850}

## site names
SITES = ["theron_hills", "provider", "argentina_range"]

## site shapefile
filepath = "/home/users/lewis.davies/british_antarctic_survey/East Station Met AOI"
shape_name = '*.shp'
SHAPES = glob.glob(f'{filepath}/{shape_name}')


## runs other functions
def main():

    ## get site shapefile
    for shape_file in SHAPES:
        
        ## loops through sites
        for index, site in enumerate(SITES):

            ## path for csv files
            csv_file = Path(f"csv_ouputs/cloud_base_{site}_stats.csv") 
            
            ## checks if csv file exists
            if csv_file.is_file():
            
                continue
            
            ## ...if it doesn't - produces data
            else:
                
                ## get data from era5
                get_stats(index, site, shape_file)
        
    return


##  get data from era 5
def get_stats(index, site, shape_file):

    ## empty stats
    stats = {"full_date": [], "month": [], "day": [], 
                 "hour": [], "latitude": [], "longitude": [],
                 'cloud_base': []}
    
    ## loops through cloud files     
    for file in CLOUD_FILES:
                
        ## load cube
        cube = iris.load_cube(file)
                
        ## mask cube to polygon based on site
        site_press_cube = mask_cube(index, cube, shape_file)
             
        ## loop through cubes by time   
        for cube in site_press_cube.slices_over('time'):
            
            ## gt time values
            time_coord = site_cube.coord('time')
            full_dt = time_coord.units.num2pydate(time_coord.points)[0]
            month_str = full_dt.strftime("%m")
            day_str = full_dt.strftime("%d")
            hour_str = full_dt.strftime("%H")
            
            for cube in cube.slices_over(['longitude', 'latitude']):
                
                ## spressure level constraints 
                con_950 = iris.Constraint(pressure=950)
                con_975 = iris.Constraint(pressure=975)
                con_1000 = iris.Constraint(pressure=1000)
                    
                ## Separate cubes by pressure levels        
                cube_950 = cube.extract(con_950)
                cube_975 = cube.extract(con_975)
                cube_1000 = cube.extract(con_1000)
                
                ## get cube data        
                cube_950_data = cube_950.data.flatten()
                cube_975_data = cube_975.data.flatten()
                cube_1000_data = cube_1000.data.flatten()
      
                ## get height from cube data 
                cloud_base = get_height(cube_950_data, cube_975_data, cube_1000_data)
                lat = cube.coord('latitude').points[0]
                lon = cube.coord('longitude').points[0]
                
                ## move to next iteration if value is masked
                if cloud_base = 'masked value':
                    
                    continue
                
                ## add values to dictionary
                else:
                                
                    print(full_dt)
                    stats["full_date"].append(full_dt)
                    stats["month"].append(month_str)
                    stats["day"].append(day_str)
                    stats["hour"].append(hour_str)
                    stats["longitude"].append(lon)
                    stats["latitude"].append(lat)
                    stats['cloud_base'].append(cloud_base)
                    
                 
    ## turn dict into dataframe   
    site_df = pd.DataFrame(stats)
    
    ## sort values by date
    site_df = site_df.sort_values(by='full_dt')

    ## save dataframe as csv
    site_df.to_csv(f'/home/users/lewis.davies/british_antarctic_survey/analysis_scripts/csv_ouputs/cloud_base_{site}_stats.csv')


## finds lowest cloud base height
def get_cloud_base(values, heights):
    
    ## loop through height dict
    for height, value in heights.items():
        
        ## return height if cloud above 5 oktas
        if value >= 5/8:
            
            return height 
        
        ## return height if on last height to check
        elif height == 1850:
            
            return height
        
        ## continue to next height/value pair
        else: 
            
            continue
  

## get non masked values from array  
def get_height(cube_950_data, cube_975_data, cube_1000_data):
    
    ## loop through data values
    for val_950, val_975, val_1000 in zip(cube_950_data, cube_975_data, cube_1000_data):
        
        ## heights to check
        heights = {350: val_1000, 1150: val_975, 1850: val_950}
        
        ## if  
        if type(val_950) == np.ma.core.MaskedConstant:
            
            height = 'masked value'
            
            return height
                
        else:
                                
            height = get_cloud_base(values, heights)
                                
            return height  
    

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
