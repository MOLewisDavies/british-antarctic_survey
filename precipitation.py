""" 

script to analyse precipitaion data from era5

threshold for precipiation useed = 0.2mm/hr

functions in file:

main(): runs other functions 

 """

import iris
import glob
from pathlib import Path
import seaborn as sns
import numpy as np
import pandas as pd

# Seaborn settings - makes it look nice
sns.set_style('darkgrid')

## thoughts on file outline:

# read in precip files 
# constrain to polygon map areas
# count 2 hour periods where precip <0.2mm/hr for light, 

# plots: 
# no. of days per week with flying window?
# no. of days per month with flying window?
# "heat map" of polygon with precipitation values??? averages???


precip_f_path = "/data/scratch/lewis.davies/bas/precip"
precip_f_name = "*.grib"
precip_files = sorted(glob.glob(f"{precip_f_path}/{precip_f_name}"))

## site names
SITES = ["theron_hills", "provider", "argentina_range"]


## runs other functions
def main():

    ## loops through sites
    for site in SITES:
        
        if site == "theron_hills":
            x = np.arange(-30, -25.5, 0.25)
            y = np.arange(-79.25, -78, 0.25)

        elif site == "provider":
            x = np.arange(-30, -26, 0.25)
            y = np.arange(-80.5, -80.0, 0.25)

        elif site == "argentina_range":
            x = np.arange(-45, -40.0, 0.25)
            y = np.arange(-83.0, -81.0, 0.25)

        ## path for csv files
        csv_file = Path(f"csv_ouputs/precip_{site}_stats.csv") 
        
        ## checks if csv file exists
        if csv_file.is_file():
            
            
            continue


        # ...if it doesn't - produces data
        else:
        
            get_stats(site, x, y)
        
    return


##  get data from era 5
def get_stats(site, x, y):

    stats = {"full_date": [], "month": [], "day": [], 
                 "hour": [], "grid_square": [], "precip": []}
    
    """  ## number of grid points for each site
    
    if site == "theron_hills":
        
        grid_points = 21

    elif site == "provider":
        
        grid_points = 11

    elif site == "argentina_range":
        
        grid_points = 58 """
    
    
    for file in precip_files:

        ## load file
        precip_cube = iris.load_cube(file)

        ## constrain cube to site location
        area_con = iris.Constraint(longitude = x, latitude = y)
        site_cube = precip_cube.extract(area_con)

        ## loopp through cube times
        for site_cube in site_cube.slices_over('time'):
            
            ## get time values
            time_coord = site_cube.coord('time')
            full_dt = time_coord.units.num2pydate(time_coord.points)[0]
            month_str = full_dt.strftime("%m")
            day_str = full_dt.strftime("%d")
            hour_str = full_dt.strftime("%H")
            
            for grid, precip in enumerate(site_cube.data.flatten()):
                
                print(full_dt, grid)
                
                stats["full_date"].append(full_dt)
                stats["month"].append(month_str)
                stats["day"].append(day_str)
                stats["hour"].append(hour_str)
                stats["grid_square"].append(grid)
                stats["precip"].append(precip)
                        
                        

    site_df = pd.DataFrame(stats)
    
    site_df = site_df.sort_values(by='full_dt')

    site_df.to_csv(f'/home/users/lewis.davies/british_antarctic_survey/analysis_scripts/csv_ouputs/precip_{site}_stats_test.csv')


if __name__ == "__main__":
    main()


