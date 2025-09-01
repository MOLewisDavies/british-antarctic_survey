""" 

script to analyse precipitaion data from era5

threshold for precipiation used = 0.2mm/hr

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
from iris.util import mask_cube_from_shapefile
from cartopy.io.shapereader import Reader

# Seaborn settings - makes it look nice
sns.set_style('darkgrid')

## precipitation files
precip_f_path = "/data/scratch/lewis.davies/bas/precip"
precip_f_name = "*.grib"
precip_files = sorted(glob.glob(f"{precip_f_path}/{precip_f_name}"))

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
            csv_file = Path(f"csv_ouputs/precip_{site}_stats_grid_test.csv") 
            
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
             "hour": [], "grid_square": [], 'latitude': [], 
             'longitude': [], "precip": []}
    
    ## loop through precipitation files
    for file in precip_files:

        ## load file
        precip_cube = iris.load_cube(file)

        ## constrain cube to site location
        site_cube = mask_cube(index, precip_cube, shape_file)

        ## loopp through cube times
        for site_cube in site_cube.slices_over('time'):
            
            ## get time values
            time_coord = site_cube.coord('time')
            full_dt = time_coord.units.num2pydate(time_coord.points)[0]
            month_str = full_dt.strftime("%m")
            day_str = full_dt.strftime("%d")
            hour_str = full_dt.strftime("%H")
            
            values = []
            lat_values = []
            lon_values = []
            
            for site_cube in site_cube.slices_over(['longitude', 'latitude']):
                
                lats = site_cube.coord('latitude').points
                lons = site_cube.coord('longitude').points
                
                ## loop through flattened data array 
                for lat, lon, precip in zip(lats, lons, site_cube.data.flatten()):

                    
                    ## if value is masked, do not add to dict
                    if type(precip) == np.ma.core.MaskedConstant:
                        
                        pass
                    
                    else:
                        
                        ## add values to list
                        values.append(precip)
                        lat_values.append(lat)
                        lon_values.append(lon)
                    

                    
            for lon, lat, precip in zip(lon_values, lat_values, values):
                
                print(full_dt)
                stats["full_date"].append(full_dt)
                stats["month"].append(month_str)
                stats["day"].append(day_str)
                stats["hour"].append(hour_str)
                stats["longitude"].append(lon)
                stats["latitude"].append(lat)
                stats["precip"].append(precip)   
    
    ## turn dict into dataframe   
    site_df = pd.DataFrame(stats)

    ## sort values by date
    site_df = site_df.sort_values(by='full_date')

    ## save dataframe as csv
    site_df.to_csv(f'/home/users/lewis.davies/british_antarctic_survey/analysis_scripts/csv_ouputs/precip_{site}_stats.csv')


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


