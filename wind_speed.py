from datetime import datetime
import iris
import iris.coord_categorisation as ic
import numpy as np
import pandas as pd
import glob

DATA_DIR = "/data/scratch/lewis.davies/bas/wind"

U_FILE_NAME = 'uwind_*.grib'
V_FILE_NAME = 'vwind_*.grib'

U_WIND_FILES = sorted(glob.glob(f'{DATA_DIR}/{U_FILE_NAME}'))
V_WIND_FILES = sorted(glob.glob(f'{DATA_DIR}/{V_FILE_NAME}'))

## site names
SITES = ["theron_hills", "provider", "argentina_range"]


def main():
    
    ## loops through sites
    for site in SITES:
        
        if site == "theron_hills":
            x = np.arange(-30.25, -25.5, 0.25)
            y = np.arange(-79.25, -78.5, 0.25)

        elif site == "provider":
            x = np.arange(-29.75, -26.5, 0.25)
            y = np.arange(-80.5, -80.0, 0.25)

        elif site == "argentina_range":
            x = np.arange(-44.25, -40.0, 0.25)
            y = np.arange(-83.0, -81.5, 0.25)

        ## path for csv files
        csv_file = Path(f"csv_ouputs/wind_speed_{site}_stats.csv")

        ## checks if csv file exists
        if csv_file.is_file():
            
            
            continue


        # ...if it doesn't - produces data
        else:
        
            get_stats(site, x, y)

    
    return


def get_stats(site, x, y):

    stats = {"full_date": [], "month": [], "day": [], 
             "hour": [], "grid_square": [], "wind_speed": []}
    
    for u_file, v_file in zip(U_WIND_FILES, V_WIND_FILES):

        u_cube = iris.load_cube(u_file)
        v_cube = iris.load_cube(v_file)

        site_con = iris.Constraint(longitude = x, latitude = y)

        u_site_cube = u_cube.extract(site_con)
        v_site_cube = v_cube.extract(site_con)

        ws_cube = calculate_wind_speed(u_site_cube, v_site_cube)

        for ws_cube in ws_cube.slices_over("time"):

            time_coord = ws_cube.coord('time')
            full_dt = time_coord.units.num2pydate(time_coord.points)[0]
            month_str = full_dt.strftime("%m")
            day_str = full_dt.strftime("%d")
            hour_str = full_dt.strftime("%H")

            for grid, wind_speed in enumerate(ws_cube.flatten()):

                stats["full_date"].append(full_dt)
                stats["month"].append(month_str)
                stats["day"].append(day_str)
                stats["hour"].append(hour_str)
                stats["grid_square"].append(grid)
                stats["wind_speed"].append(wind_speed)

    site_df = pd.DataFrame(stats)

    site_df = site_df.sort_values(by='full_dt')
    
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


if __name__ == "__main__":
    main()
