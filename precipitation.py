"""

script to analyse precipitaion data from era5

threshold for precipiation used = 0.1mm/hr

functions in file:

main(): runs other functions
process_data(): process era 5 files
get_stats(): gets data from era5 cube
run_all_points_heatmap_func(): creates heatmap for all era5 grid points
make_box_whisker_plots(): makes box plots
make_box_whisker_east_west(): make box plots for east and west of sites

"""

import glob
import pickle
import warnings
from pathlib import Path

import cartopy
import cartopy.crs as ccrs
import iris
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.collections import PatchCollection

## import fuctions and constants
import functions as func
from constants import BAS_PATH, SHP_FILE, SITES, YEARS, MONTHS, MONTH_DAYS

warnings.filterwarnings("ignore")

# Seaborn settings - makes it look nice
sns.set_style("darkgrid")

## precipitation files
PRECIP_F_PATH = "/data/scratch/lewis.davies/bas/precip"
PRECIP_F_NAME = "*.grib"
PRECIP_FILES = sorted(glob.glob(f"{PRECIP_F_PATH}/{PRECIP_F_NAME}"))


## runs other functions
def main():

    """ 
    Runs all functions within file.
    """
    
    #process_data()
    func.make_heatmap_plots('Precip')  
    
    
    #get_stats_all_points()
    #run_all_points_heatmap_func()

    # make_heatmap_plots()
    # make_box_whisker_plots('TOP THREE', 'AVERAGE')
    # make_box_whisker_plots('TOP THREE', 'EVERY')
    # make_box_whisker_plots('ALL', 'AVERAGE')
    # make_box_whisker_plots('ALL', 'EVERY')
    # make_box_whisker_eastload sight specific _west('EAST', 'AVERAGE')
    # make_box_whisker_east_west('EAST', 'EVERY')
    # make_box_whisker_east_west('WEST', 'AVERAGE')
    # make_box_whisker_east_west('WEST', 'EVERY')

    return


## run get stats functions
def process_data():
    """
    Processes era5 files into stats csv files.
    """

    ## loops through sites
    for index, site in enumerate(SITES):

        ## path for csv files
        csv_file = Path(f"csv_ouputs/Precip_{site}_stats.csv")

        ## checks if csv file exists
        if csv_file.is_file():

            print(f"Precip_{site}_stats.csv exists")

        # ...if it doesn't - produces data
        else:

            ## get stats from cubes
            get_stats(index, site, SHP_FILE)


## get data from era 5
def get_stats_all_points():
    """ 
    Creates dataframe of all grid points in ourly era5 data.
    
    """


    big_df = pd.DataFrame()

    for file in PRECIP_FILES:

        ## load file
        precip_cube = iris.load_cube(file)

        # Extract coordinates
        time_points = precip_cube.coord("time").points
        time = precip_cube.coord("time").units.num2date(time_points)
        lat = precip_cube.coord("latitude").points
        lon = precip_cube.coord("longitude").points

        # Create full coordinate meshgrid
        time_grid, lat_grid, lon_grid = np.meshgrid(time, lat, lon, indexing="ij")

        # Put data into dataframe
        df = pd.DataFrame(
            {
                "Date and Time": time_grid.flatten(),
                "Latitude": lat_grid.flatten(),
                "Longitude": lon_grid.flatten(),
                "Precip": precip_cube.data.flatten(),
            }
        )

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
        f"{BAS_PATH}/csv_ouputs/precip_stats_all_points.csv"
    )


## get data from era 5
def get_stats(index, site, shape_file):
    """ 
    Creates dataframe of hourly era5 weather data, separated into grid
    point.
    
    Args:
        
        index (int): index of shape file in list
        site (str): Antarctica site.
        shape_file (path): path to shape file
    """
    
    # Empty dataframe to concatenate to
    big_df = pd.DataFrame()

    ## loop through precipitation files
    for file in PRECIP_FILES:

        ## load file
        precip_cube = iris.load_cube(file)

        ## constrain cube to site location
        site_cube = func.mask_cube(index, precip_cube, shape_file)

        # Extract coordinates
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
                "Precip": site_cube.data.flatten(),
            }
        )

        # Drop masked values
        df = df.dropna(subset=["Precip"])

        # Add in month, day and hour columns
        df["Year"] = [d_time.year for d_time in df["Date and Time"]]
        df["Month"] = [d_time.month for d_time in df["Date and Time"]]
        df["Day"] = [d_time.day for d_time in df["Date and Time"]]
        df["Hour"] = [d_time.hour for d_time in df["Date and Time"]]

        # Concatenate to big dataframe
        big_df = pd.concat([big_df, df], ignore_index=True)

    ## sort values by date
    big_df = big_df.sort_values(by="Date and Time")

    # print(big_df)

    ## save dataframe as csv
    big_df.to_csv(
        f"{BAS_PATH}/csv_ouputs/Precip_{site}_stats.csv"
    )


## make all points heatmap
def run_all_points_heatmap_func():
    """ 
    Runs count_flyable_days_all_points, to make_all_points_heatmap.
    
    """

    ## create area dictioanry
    area_dict = {
        "lat": [],
        "lon": [],
        "jan": [],
        "feb": [],
        "oct": [],
        "nov": [],
        "dec": [],
    }

    ## load dataframe
    data_df = pd.read_csv("/data/scratch/lewis.davies/bas/csv_ouputs/precip_stats_all_points.csv")

    ## create flyable file
    flyable_file = count_flyable_days_all_points(data_df)

    ## open flyable daya dict file
    with open(flyable_file, "rb") as file:

        flyable_days_dict = pickle.load(file)

    ## turn into dataframe
    flyable_df = pd.DataFrame(flyable_days_dict)

    print(flyable_df)

    ## loop through dataframe rows
    for index, row in flyable_df.iterrows():

        ## loop through dicitonay keys
        for var in area_dict.keys():

            ## append value from row to dictionary for area heatmap
            area_dict[f"{var}"].append(row[f"{var}"].item())

    ## create dataframe from area dictionary
    area_df = pd.DataFrame(area_dict)

    ## make heatmap of whole area
    make_all_points_heatmap(area_df)


## make heatmap of whole area
def make_all_points_heatmap(area_df):
    """ 
    Makes heatmap for all grid points in data.
    
    ARGS:
    
        area_df (DataFrame): weather data
    """

    # month list for creating plots
    months = ["total", "jan", "feb", "oct", "nov", "dec"]

    ## create total column
    area_df["total"] = (
            area_df["jan"]
            + area_df["feb"]
            + area_df["oct"]
            + area_df["nov"]
            + area_df["dec"]
        )
    
    
    ## get longitude and latitudes
    lons, lats = area_df["lon"], area_df["lat"]

    ## set plot aspect
    lon_width = max(lons) - min(lons)
    lat_height = max(lats) - min(lats)
    plot_aspect = lat_height / lon_width

    for month in months:
        
        # set up axes
        fig = plt.figure(figsize=(12, 12 * plot_aspect))
        ax = plt.axes(projection=ccrs.PlateCarree())

        # Set extent of plot
        ax.set_extent(
            [min(lons) - 0.5, max(lons) + 0.5, min(lats) - 0.5, max(lats) + 0.5],
            ccrs.PlateCarree(),
        )

        # Draw coastlines and background stuff
        ax.add_feature(cartopy.feature.LAND)
        ax.add_feature(cartopy.feature.OCEAN)
        ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.2)
        ax.add_feature(cartopy.feature.BORDERS, linestyle=":", linewidth=0.3)
        ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
        ax.add_feature(cartopy.feature.RIVERS)

        ## define precip data
        precip_data = (area_df[f"{month}"] / MONTH_DAYS[f'{month}'])*100

        ## create normalised data column
        area_df["norm"] = (area_df[f"{month}"] - area_df[f"{month}"].min()) / (
            area_df[f"{month}"].max() - area_df[f"{month}"].min()
        )
        precip_norm = area_df["norm"]

        ## get colour map
        cmap = cm.get_cmap("viridis")

        ## create patch collection for making heatmap
        patches = [
            plt.Rectangle(
                (lon - 0.125, lat - 0.125), 0.25, 0.25, alpha=1, facecolor=cmap(precip)
            )
            for lon, lat, precip in zip(lons, lats, precip_norm)
        ]
        patch_collection = PatchCollection(patches, match_original=True)

        ## create scalar mappable for colourbar
        norm = plt.Normalize(precip_data.min(), precip_data.max())
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

        ## add heatmap patch collection to axes
        ax.add_collection(patch_collection)

        ## add title
        ax.set_title(f"{month.capitalize()} Flyable Days by Grid Point: All Sites")

        ## add colourbar
        plt.colorbar(
            sm,
            ax=ax,
            label="No. of Flyable Days",
            ticks=[precip_data.min(), precip_data.median(), precip_data.max()],
        )

        ## save figure
        fig.savefig(
            f"/{BAS_PATH}/plots/{month}_all_data_precip_heatmap.png"
        )


## count flyable days per year+month for box plots
def count_flyable_days_yearly(data_df, site):
    """ 
    Creates flyable file to be used for making box plots.
    
    ARGS:
    
        data_df (DataFrame): weather data.
        site (str): Antarctica site.
        
    Returns:
    
        .txt file: pickled dictionary file path.
    """
    
    
    ## define flyable file path
    flyable_file = Path(f"csv_ouputs/precip_{site}_fly_days_years.txt")

    ## check if file exists
    if flyable_file.is_file():

        ## if exists, return file
        return flyable_file

    ## otherwise
    else:

        ## set up flyable dictionary
        flyable_days_dict = {
            "lat": [],
            "lon": [],
            "year": [],
            "month": [],
            "flyable_days": [],
        }

        ## get list of unique lat lon pairs
        lat_lon_list = data_df[["Latitude", "Longitude"]].drop_duplicates()

        ## convert into individual lat lon lists
        lats = lat_lon_list["Latitude"].tolist()
        lons = lat_lon_list["Longitude"].tolist()

        ## loop through years
        for year in YEARS:

            ## subset dataframe for year
            year_df = data_df.loc[data_df["Year"] == year]

            ## loop through lat and lon lists
            for lat, lon in zip(lats, lons):

                ## constrain by latitude and longitude
                lat_con = year_df["Latitude"] == lat
                lon_con = year_df["Longitude"] == lon
                lat_lon_df = year_df.loc[lat_con & lon_con]

                ## loop through months
                for m_name, m_number in MONTHS.items():

                    ## get correct number of days
                    if m_number == 1 or m_number == 10 or m_number == 12:

                        ## define length of specific month
                        days = np.arange(1, 32, 1)

                    ## get correct number of days
                    elif m_number == 2:

                        ## define length of specific month
                        days = np.arange(1, 29, 1)

                    ## else:
                    else:

                        ## define length of specific month
                        days = np.arange(1, 31, 1)

                    ## constrain by month
                    month_df = lat_lon_df.loc[data_df["Month"] == m_number]
                    
                    flyable_days = 0
                    
                    ## loop through days in month
                    for day in days:

                        ## subset by day
                        day_df = month_df.loc[month_df["Day"] == day]

                        ## get True or False for a flyable day
                        flyable = func.is_day_flyable(day_df, 'Precip')

                        ## check if day is flyable after analysis
                        if flyable:

                            ## add 1 to total
                            flyable_days += 1

                    ## add total for month and year to dict
                    flyable_days_dict["month"].append(m_number)
                    flyable_days_dict["flyable_days"].append(flyable_days)
                    flyable_days_dict["year"].append(year)
                    flyable_days_dict["lat"].append(lat)
                    flyable_days_dict["lon"].append(lon)

        ## pickle out dictionary
        with open(f"csv_ouputs/precip_{site}_fly_days_years.txt", "wb") as file:

            pickle.dump(flyable_days_dict, file, protocol=2)

    ## return flyable file
    return flyable_file


## count flyable days for all grid points in era5 files
def count_flyable_days_all_points(data_df):
    """ 
    Creates flyable days file for all points in era5 data files..
    
    ARGS:
    
        data_df (DataFrame): weather data.
        
    Returns:
    
        .txt file: pickled dictionary file path.
    """
    
    
    ## define flyable file path
    flyable_file = Path("csv_ouputs/precip_all_points_fly_days.txt")

    ## check if file exists
    if flyable_file.is_file():

        ## if exists, return file
        return flyable_file

    ## otherwise
    else:

        ## set up flyable dictionary
        flyable_days_dict = {
            "lat": [],
            "lon": [],
            "jan": [],
            "feb": [],
            "oct": [],
            "nov": [],
            "dec": [],
        }

        ## get list of unique lat lon pairs
        lat_lon_list = data_df[["Latitude", "Longitude"]].drop_duplicates()

        ## convert into individual lat lon lists
        lats = lat_lon_list["Latitude"].tolist()
        lons = lat_lon_list["Longitude"].tolist()

        ## loop through lists
        for lat, lon in zip(lats, lons):

            ## constrain data df by lat and lon values
            lat_con = data_df["Latitude"] == lat
            lon_con = data_df["Longitude"] == lon
            lat_lon_df = data_df.loc[lat_con & lon_con]

            ## loop through months
            for m_name, m_number in MONTHS.items():

                ## get correct number of days
                if m_number == 1 or m_number == 10 or m_number == 12:

                    ## define length of specific month
                    days = np.arange(1, 32, 1)

                ## get correct number of days
                elif m_number == 2:

                    ## define length of specific month
                    days = np.arange(1, 29, 1)

                ## else:
                else:

                    ## define length of specific month
                    days = np.arange(1, 31, 1)

                ## constrain by month
                month_df = lat_lon_df.loc[data_df["Month"] == m_number]

                ## define value for flyable list
                flyable_days = 0

                ## loop through years
                for year in YEARS:

                    ## subset dataframe for year
                    year_month_df = month_df.loc[month_df["Year"] == year]

                    ## loop through days in month
                    for day in days:

                        ## subset by day
                        day_df = year_month_df.loc[year_month_df["Day"] == day]

                        ## get True or False for a flyable day
                        flyable = func.is_day_flyable(day_df, 'Precip')

                        ## check if day is flyable after analysis
                        if flyable:

                            ## add 1 to total
                            flyable_days += 1

                print(m_name, flyable_days)
                ## add total for month and year to dict
                flyable_days_dict[m_name].append(flyable_days)

            ## append to dict
            flyable_days_dict["lat"].append(lat)
            flyable_days_dict["lon"].append(lon)

        ## pickle out dictionary
        with open("csv_ouputs/precip_all_points_fly_days.txt", "wb") as file:

            pickle.dump(flyable_days_dict, file, protocol=2)

    ## return flyable file
    return flyable_file


## makes box plot for subsets of grid points
def make_box_whisker_plots(grid_points, d_type):
    """ 
    Makes box and whisker plots for  different sets of grid points.
    
    ARGS:
    
        grid_points (str): grid_point subset (ALL or TOP THREE)
        d_type (str): analysis method (EVERY or AVERAGE)
    """
    
    
    ## set up figure
    fig, ax = plt.subplots(figsize=(10, 6))

    ## set up empty df list
    df_list = []

    ## loop through sites
    for site in SITES:

        ## load data
        data_df = pd.read_csv(f"csv_ouputs/precip_{site}_stats.csv")

        ## count flyable days and pickle file
        flyable_file = count_flyable_days_yearly(data_df, site)

        ## open flyable days file
        with open(flyable_file, "rb") as file:

            ## load flyable dictionary
            flyable_days_dict = pickle.load(file)

        ## load flyable
        flyable_df = pd.DataFrame(flyable_days_dict)
        flyable_df["site"] = site

        ## loopp through years
        for year in YEARS:

            ## constrain by year
            year_df = flyable_df.loc[flyable_df["year"] == year]

            ## loop through month name and number
            for m_name, m_number in MONTHS.items():

                ## constrain by month number
                year_month_df = year_df.loc[year_df["month"] == m_number]

                ## check grid_points input
                if grid_points == "ALL":

                    ## run if input is EVERY
                    if d_type == "EVERY":

                        ## add df to list
                        df_list.append(year_month_df)

                    ## run if input is AVERAGE
                    elif d_type == "AVERAGE":

                        ## set up average dictionary
                        average_dict = {
                            "year": [],
                            "month": [],
                            "flyable_days": [],
                            "site": [],
                        }

                        ## create list of flyable day values
                        flyable_days = year_month_df["flyable_days"]

                        ## append data to average dictionary
                        average_dict["year"].append(year)
                        average_dict["month"].append(m_number)
                        average_dict["flyable_days"].append(np.mean(flyable_days))
                        average_dict["site"].append(site)

                        ## create dataframe from dict
                        average_df = pd.DataFrame(average_dict)

                        ## append dictionary to list
                        df_list.append(average_df)

                ## check grid points input
                elif grid_points == "TOP THREE":

                    ## run if EVERY
                    if d_type == "EVERY":

                        ## get largest 3 flyable values from dataframe
                        largest_df = year_month_df.nlargest(3, "flyable_days")

                        ## append df to list
                        df_list.append(largest_df)

                    ## run if AVERAGE
                    if d_type == "AVERAGE":

                        ## set up average dictionary
                        average_dict = {
                            "year": [],
                            "month": [],
                            "flyable_days": [],
                            "site": [],
                        }

                        ## get largest 3 values from dict
                        largest_df = year_month_df.nlargest(3, "flyable_days")
                        largest_values = largest_df["flyable_days"]

                        ## append values to dictionary
                        average_dict["year"].append(year)
                        average_dict["month"].append(m_number)
                        average_dict["flyable_days"].append(np.mean(largest_values))
                        average_dict["site"].append(site)

                        ## create dictionary
                        average_df = pd.DataFrame(average_dict)

                        ## append dataframe to list
                        df_list.append(average_df)

    ## concat dataframe list
    area_df = pd.concat(df_list)

    ## create box plot
    sns.boxplot(data=area_df, x="month", y="flyable_days", hue="site", ax=ax)

    ## set title
    ax.set_title(
        f"Flyable Days Per Month - {grid_points} grid point(s) - {d_type}",
        fontsize=18,
        fontweight="bold",
    )

    ## save figure
    fig.savefig(
        f"{BAS_PATH}/plots/monhtly_precip_whisker_plot_{grid_points}_{d_type}.png"
    )


## makes box plot for east/west area split
def make_box_whisker_east_west(AREA, d_type):
    """ 
    Makes box and whisker plots for East and West grid points.
    
    ARGS:
    
        AREA (str): analysis area (EAST or WEST)
        d_type (str): analysis method (EVERY or AVERAGE)
    """


    ## set up plot
    fig, ax = plt.subplots(figsize=(10, 6))

    ## empty datframe list for concatenation
    df_list = []

    ## loop through sites
    for site in SITES:

        ## load file
        data_df = pd.read_csv(f"csv_ouputs/precip_{site}_stats.csv")

        ## count flyable days per year+month
        flyable_file = count_flyable_days_yearly(data_df, site)

        ## open flyable days dict file
        with open(flyable_file, "rb") as file:

            ## opens created/saved flyable file
            flyable_days_dict = pickle.load(file)

        ## turn dict into dataframe
        flyable_df = pd.DataFrame(flyable_days_dict)

        ## create site column
        flyable_df["site"] = site

        ## loop through years
        for year in YEARS:

            ## constrain to year
            year_df = flyable_df.loc[flyable_df["year"] == year]

            ## loop through months
            for m_name, m_number in MONTHS.items():

                ## constrain by month number
                year_month_df = year_df.loc[year_df["month"] == m_number]

                ## get latitude list
                latitudes = year_month_df["lat"]

                ## check area input variable
                if AREA == "WEST":

                    ## constrain to western side of site
                    ## smallest 1/2 of lats = western half of shape
                    west_df = year_month_df.nsmallest(round(len(latitudes) / 2), "lat")

                    ## create list of flyable days
                    flyable_days = west_df["flyable_days"]

                    ## check d_type input variable
                    if d_type == "EVERY":

                        ## add df to dataframe list
                        df_list.append(west_df)

                    ## check d_type input variable
                    elif d_type == "AVERAGE":

                        ## set up dict for average values
                        average_dict = {
                            "year": [],
                            "month": [],
                            "flyable_days": [],
                            "site": [],
                        }

                        ## append values to average dict
                        average_dict["year"].append(year)
                        average_dict["month"].append(m_number)
                        average_dict["flyable_days"].append(np.mean(flyable_days))
                        average_dict["site"].append(site)

                        ## turn dict into dataframe
                        average_df = pd.DataFrame(average_dict)

                        ## append dataframe to list
                        df_list.append(average_df)

                ## check inpput variable
                elif AREA == "EAST":

                    ## constrain to eastern side of site
                    ## largest 1/2 of lats = western half of shape
                    east_df = year_month_df.nlargest(round(len(latitudes) / 2), "lat")

                    ## create list of flyable days
                    flyable_days = east_df["flyable_days"]

                    ## check area input variable
                    if d_type == "EVERY":

                        ## add dataframe to list
                        df_list.append(east_df)

                    ## check area input variable
                    elif d_type == "AVERAGE":

                        ## set up dict for average values
                        average_dict = {
                            "year": [],
                            "month": [],
                            "flyable_days": [],
                            "site": [],
                        }

                        ## append values to average dict
                        average_dict["year"].append(year)
                        average_dict["month"].append(m_number)
                        average_dict["flyable_days"].append(np.mean(flyable_days))
                        average_dict["site"].append(site)

                        ## turn dict into dataframe
                        average_df = pd.DataFrame(average_dict)

                        ## append dataframe to list
                        df_list.append(average_df)

    ## concat list into one dataframe
    area_df = pd.concat(df_list)

    ## create box plot
    sns.boxplot(data=area_df, x="month", y="flyable_days", hue="site", ax=ax)

    ## set title
    ax.set_title(
        f"Flyable Days Per Month - {AREA}ERN grid point(s) - {d_type}",
        fontsize=18,
        fontweight="bold",
    )

    ## save figure
    fig.savefig(
        f"{BAS_PATH}/plots/monhtly_precip_whisker_plot_{AREA}_{d_type}.png"
    )


if __name__ == "__main__":
    main()
print("finished")
