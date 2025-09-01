""" 

script to extract era5 temperature and water content data 
and calculate visibility

functions in file:

main(): runs other functions 
get_stats(): gets data from era5 cube
conver_wc_units(): converts water content units
calculate_visibility(): calculate visibility approximation
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
TEMP_F_PATH = "/data/scratch/lewis.davies/bas/temperature"
VIS_F_PATH = "/data/scratch/lewis.davies/bas/visibility"

LWC_F_NAME = "liq_water_*"
IWC_F_NAME = "ice_water_*"
TEMP_F_NAME = "*.grib"

TEMP_FILES = sorted(glob.glob(f"{TEMP_F_PATH}/{TEMP_F_NAME}"))
LWC_FILES = sorted(glob.glob(f"{LWC_F_PATH}/{LWC_F_NAME}"))
IWC_FILES = sorted(glob.glob(f"{IWC_F_PATH}/{IWC_F_NAME}"))

## site shapefile
filepath = "/home/users/lewis.davies/british_antarctic_survey/East Station Met AOI"
shape_name = '*.shp'
SHAPES = glob.glob(f'{filepath}/{shape_name}')


def main():
    
    ## get shapefile
    for shape_file in SHAPES:
        
        ## loops through sites
        for index, site in enumerate(SITES):

            ## path for csv files
            csv_file = Path(f"csv_ouputs/visibility_{site}_stats.csv")

            ## checks if csv file exists
            if csv_file.is_file():
                
                
                continue


            # ...if it doesn't - produces data
            else:
            
                ## get stats from cubes
                get_stats(index, site, shape_file)
                
                
                
## extract data from era5 files   
def get_stats(index, site, shape_file):
    
    ## empty stats
    stats = {"full_date": [], "month": [], "day": [], 
             "hour": [], 'latitude': [], 'longitude': [],
             "vis": []}
    
    ## loop through lwc, iwc and temp files together
    for temp_file, lwc_file, iwc_file in zip(TEMP_FILES, LWC_FILES, IWC_FILES):
        
        ## loadd cubes
        temp_cube = iris.load_cube(temp_file)
        lwc_cube = iris.load_cube(lwc_file)
        iwc_cube = iris.load_cube(iwc_file)
        
        ## calculate visibility using cubes
        vis_cube = calculate_visibility(temp_cube, lwc_cube, iwc_cube)
        
        ## mask cube for correct site using shape_file
        site_cube = mask_cube(index, vis_cube, shape_file)

        for site_cube in site_cube.slices_over('time')
            
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
                for lat, lon, vis in zip(lats, lons, site_cube.data.flatten()):

                    
                    ## if value is masked, do not add to dict
                    if type(vis) == np.ma.core.MaskedConstant:
                        
                        pass
                    
                    else:
                        
                        ## add values to list
                        values.append(precip)
                        lat_values.append(lat)
                        lon_values.append(lon)
            
            
            
            ## loop through flattened data array
            for lon, lat, vis in zip(lon_values, lat_values, values):
                    
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
                    stats["longitude"].append(lon)
                    stats["latitude"].append(lat)
                    stats["vis"].append(vis)
    
    
    ## turn dict into dataframe   
    site_df = pd.DataFrame(stats)
    
    ## sort values by date
    site_df = site_df.sort_values(by='full_dt')

    ## save dataframe as csv
    site_df.to_csv(f'/home/users/lewis.davies/british_antarctic_survey/analysis_scripts/csv_ouputs/visibility_{site}_stats.csv')   
    
    return


def convert_wc_units(temp_cube, wc):
    
    R = 287.05
    
    air_density = 100000/(R * temp_cube.data) ## convert pressure units

    converted_data = (wc.data * air_density)*1000 ## converts kg/kg to gm-3
    
    return converted_data


def calculate_visibility(temp_cube, lwc_cube, iwc_cube):
    
    lwc_data = convert_wc_units(temp_cube, lwc_cube)
    iwc_data = convert_wc_units(temp_cube, iwc_cube)
    
    N = 50 
    e = 0.02

    liq_ext_coeff = -(math.log(e)/1.002e3)*((lwc_data*N)**0.6473)
    ice_ext_coeff = (163.0e-3)*iwc_data
    rain_ext_coeff = (5.0e-3)*(lwc_data**0.75)
    snow_ext_coeff = (4.0e-3)*(iwc_data**0.78)
    
    cld_ext_coeff = liq_ext_coeff + ice_ext_coeff + snow_ext_coeff + rain_ext_coeff
    air_ext_coeff = (math.log(e))/(10**5)
    
    total_ext_coeff = - air_ext_coeff + cld_ext_coeff 
    
    visibility_data = -(math.log(e))/total_ext_coeff
    
    vis_cube = temp_cube.copy(visibility_data)
    
    vis_cube.standard_name = 'visibility'
    
    return vis_cube


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
