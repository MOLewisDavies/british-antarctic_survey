"""

script to extract data from era5 u and v wind files and
calculate wind speed

functions in file:

main(): runs other functions
process_data(): runs get stats loop
make heatmap plots(): runs functions for producing plots
make_area_heatmap(): makes heatmap of entire region()
make_site_heatmap(): makes site specific heatmaps
get_stats(): gets data from era5 cube
calculate_wind_speed(): calculates wind speed
mask_cube(): masks cube based on shapefile
count_flyable_days(): creates table of flyable  days per month
is_day_flyable(): checks data against thresholds

"""

import glob
import heapq
import pickle
import warnings
from datetime import datetime, timedelta
from pathlib import Path

import cartopy
import cartopy.crs as ccrs
import iris
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns
from cartopy.io.shapereader import Reader
from cartopy.mpl.gridliner import LongitudeFormatter
from iris.util import mask_cube_from_shapefile
from matplotlib.axis import Axis
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

warnings.filterwarnings("ignore")

# Seaborn settings - makes it look nice
sns.set_style("darkgrid")

MONTHS = {"jan": 1, "feb": 2, "oct": 10, "nov": 11, "dec": 12}

YEARS = np.arange(1993, 2024, 1)

## file path
DATA_DIR = "/data/scratch/lewis.davies/bas/wind"

## file name structure
U_FILE_NAME = "uwind_*.grib"
V_FILE_NAME = "vwind_*.grib"

## files
U_WIND_FILES = sorted(glob.glob(f"{DATA_DIR}/{U_FILE_NAME}"))
V_WIND_FILES = sorted(glob.glob(f"{DATA_DIR}/{V_FILE_NAME}"))

## site names
SITES = ["theron_hills", "provider", "argentina_range"]

## shapefile path
filepath = "/home/users/lewis.davies/british_antarctic_survey/East Station Met AOI"
shape_name = "*.shp"
SHAPES = glob.glob(f"{filepath}/{shape_name}")


## runs other functions
def main():

    process_data()
    make_heatmap_plots()

    return


## make plots
def make_heatmap_plots():

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

    ## loop through SITES
    for site in SITES:

        ## load dataframe
        data_df = pd.read_csv(f"csv_ouputs/wind_speed_{site}_stats.csv")

        ## create flyable file
        flyable_file = count_flyable_days(data_df, site)

        ## open flyable daya dict file
        with open(flyable_file, "rb") as file:

            flyable_days_dict = pickle.load(file)

        ## turn into dataframe
        flyable_df = pd.DataFrame(flyable_days_dict)

        ## loop through dataframe rows
        for index, row in flyable_df.iterrows():

            ## loop through dicitonay keys
            for var in area_dict.keys():

                ## append value from row to dictionary for area heatmap
                area_dict[f"{var}"].append(row[f"{var}"].item())

        ## make heatmap for sites based on months and totals
        make_site_heatmap(flyable_df, "month", site)
        make_site_heatmap(flyable_df, "total", site)

    ## create dataframe from area dictionary
    area_df = pd.DataFrame(area_dict)

    ## make heatmap of whole area
    make_area_heatmap(area_df)


## make heatmap of whole area
def make_area_heatmap(area_df):

    ## calculate total days for percentage calc
    days = (
        max(area_df["jan"])
        + max(area_df["feb"])
        + max(area_df["oct"])
        + max(area_df["nov"])
        + max(area_df["dec"])
    )

    ## create total column
    area_df["total"] = (
        (
            area_df["jan"]
            + area_df["feb"]
            + area_df["oct"]
            + area_df["nov"]
            + area_df["dec"]
        )
        / days
    ) * 100

    ## get longitude and latitudes
    lons, lats = area_df["lon"], area_df["lat"]

    ## set plot aspect
    lon_width = max(lons) - min(lons)
    lat_height = max(lats) - min(lats)
    plot_aspect = lat_height / lon_width

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
    wind_speed_data = area_df["total"]

    ## create normalised data column
    area_df["norm"] = (area_df["total"] - area_df["total"].min()) / (
        area_df["total"].max() - area_df["total"].min()
    )
    wind_speed_norm = area_df["norm"]

    ## get colour map
    cmap = cm.get_cmap("viridis")

    ## add gridpoint crosses to plot
    # sns.scatterplot(data = flyable_df, x=lons, y=lats, color='red', size=200, marker='+', ax=ax)

    ## create patch collection for making heatmap
    patches = [
        plt.Rectangle(
            (lon - 0.125, lat - 0.125), 0.25, 0.25, alpha=1, facecolor=cmap(wind)
        )
        for lon, lat, wind in zip(lons, lats, wind_speed_norm)
    ]

    patch_collection = PatchCollection(patches, match_original=True)

    ## create scalar mappable for colourbar
    norm = plt.Normalize(wind_speed_data.min(), wind_speed_data.max())
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

    ## add heatmap patch collection to axes
    ax.add_collection(patch_collection)

    ## add title
    ax.set_title(f"Total Flyable Days by Grid Point: All Sites")

    ## add colourbar
    plt.colorbar(
        sm,
        ax=ax,
        label="No. of Flyable Days",
        ticks=[
            wind_speed_data.min(),
            wind_speed_data.median(),
            wind_speed_data.max(),
        ],
    )

    ## save figure
    fig.savefig(
        f"/home/users/lewis.davies/british_antarctic_survey/analysis_scripts/plots/total_all_sites_wind_speed_heatmap.png"
    )


## make heatmap for each site individually
def make_site_heatmap(flyable_df, d_type, site):

    ## check if plot type is by month or all days
    if d_type == "month":

        # month list for creating plots
        months = ["jan", "feb", "oct", "nov", "dec"]

    else:

        months = ["total"]

    ## create total column
    flyable_df["total"] = (
        flyable_df["jan"]
        + flyable_df["feb"]
        + flyable_df["oct"]
        + flyable_df["nov"]
        + flyable_df["dec"]
    )

    ## get longitude and latitudes
    lons, lats = flyable_df["lon"], flyable_df["lat"]

    ## set plot aspect
    lon_width = max(lons) - min(lons)
    lat_height = max(lats) - min(lats)
    plot_aspect = lat_height / lon_width

    ## loop through months
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
        wind_speed_data = flyable_df[f"{month}"]

        ## create normalised data column
        flyable_df["norm"] = (flyable_df[f"{month}"] - flyable_df[f"{month}"].min()) / (
            flyable_df[f"{month}"].max() - flyable_df[f"{month}"].min()
        )
        wind_speed_norm = flyable_df["norm"]

        ## get colour map
        cmap = cm.get_cmap("viridis")

        ## add gridpoint crosses to plot
        # sns.scatterplot(data = flyable_df, x=lons, y=lats, color='red', size=200, marker='+', ax=ax)

        ## create patch collection for making heatmap
        patches = [
            plt.Rectangle(
                (lon - 0.125, lat - 0.125), 0.25, 0.25, alpha=1, facecolor=cmap(wind)
            )
            for lon, lat, wind in zip(lons, lats, wind_speed_norm)
        ]
        patch_collection = PatchCollection(patches, match_original=True)

        ## create scalar mappable for colourbar
        norm = plt.Normalize(wind_speed_data.min(), wind_speed_data.max())
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

        ## add heatmap patch collection to axes
        ax.add_collection(patch_collection)

        ## format site name for plot title
        if site == "theron_hills":

            form_site = "Theron Hills Site"

        elif "provider" in site:
            form_site = "Provider Site"

        else:
            form_site = "Argentina Range Site"

        ## add title
        ax.set_title(f"{month.capitalize()} Flyable Days by Grid Point: {form_site}")

        ## add colourbar
        plt.colorbar(
            sm,
            ax=ax,
            label="No. of Flyable Days",
            ticks=[
                wind_speed_data.min(),
                wind_speed_data.median(),
                wind_speed_data.max(),
            ],
        )

        ## save figure
        fig.savefig(
            f"/home/users/lewis.davies/british_antarctic_survey/analysis_scripts/plots/{month}_{site}_wind_speed_heatmap.png"
        )


## run get stats functions
def process_data():

    ## get shapefile
    for shape_file in SHAPES:

        ## loops through sites
        for index, site in enumerate(SITES):

            ## path for csv files
            csv_file = Path(f"csv_ouputs/wind_speed_{site}_stats.csv")

            ## checks if csv file exists
            if csv_file.is_file():

                print(f"wind_speed_{site}_stats.csv exists")

            # ...if it doesn't - produces data
            else:

                ## get stats from cubes
                get_stats(index, site, shape_file)


## get stats from era5 data
def get_stats(index, site, shape_file):

    big_df = pd.DataFrame()

    ## loop through u and v wind files together
    for u_file, v_file in zip(U_WIND_FILES, V_WIND_FILES):

        ## load u and v files
        u_cube = iris.load_cube(u_file)
        v_cube = iris.load_cube(v_file)

        ## calculate new wind speed cube from u and v
        ws_cube = calculate_wind_speed(u_cube, v_cube)

        ## mask wind speed cube based on site
        ws_site_cube = mask_cube(index, ws_cube, shape_file)

        # Extract coordinates
        time_points = ws_site_cube.coord("time").points
        time = ws_site_cube.coord("time").units.num2date(time_points)
        lat = ws_site_cube.coord("latitude").points
        lon = ws_site_cube.coord("longitude").points

        # Create full coordinate meshgrid
        time_grid, lat_grid, lon_grid = np.meshgrid(time, lat, lon, indexing="ij")

        # Put data into dataframe
        df = pd.DataFrame(
            {
                "Date and Time": time_grid.flatten(),
                "Latitude": lat_grid.flatten(),
                "Longitude": lon_grid.flatten(),
                "Wind Speed": ws_site_cube.data.flatten(),
            }
        )

        # Drop masked values
        df = df.dropna(subset=["Wind Speed"])

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
        f"/home/users/lewis.davies/british_antarctic_survey/analysis_scripts/csv_ouputs/wind_speed_{site}_stats.csv"
    )


## calculate wind speed from u and v cubes
def calculate_wind_speed(u_cube, v_cube):

    # multiply u and v components
    wind_speed_data = np.sqrt(u_cube.data**2 + v_cube.data**2)

    # copy to u_cube to retain other coordinates
    ws_cube = u_cube.copy(wind_speed_data)

    # change cube name
    ws_cube.standard_name = "wind_speed"

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


## count flyable days per year+month for box plots
def count_flyable_days_yearly(data_df, site):

    ## define flyable file path
    flyable_file = Path(f"csv_ouputs/wind_speed_{site}_fly_days_years.txt")

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

                    ## loop through days in month
                    for day in days:

                        ## subset by day
                        day_df = month_df.loc[month_df["Day"] == day]

                        ## get True or False for a flyable day
                        flyable = is_day_flyable(day_df, days)

                        ## check if day is flyable after analysis
                        if flyable == True:

                            ## add 1 to total
                            flyable_days += 1

                    ## add total for month and year to dict
                    flyable_days_dict["month"].append(m_number)
                    flyable_days_dict["flyable_days"].append(flyable_days)
                    flyable_days_dict["year"].append(year)
                    flyable_days_dict["lat"].append(lat)
                    flyable_days_dict["lon"].append(lon)

        ## pickle out dictionary
        with open(f"csv_ouputs/wind_speed_{site}_fly_days_years.txt", "wb") as file:

            pickle.dump(flyable_days_dict, file, protocol=2)

    ## return flyable file
    return flyable_file


## count flyable days for heatmap plots
def count_flyable_days(data_df, site):

    ## define flyable file path
    flyable_file = Path(f"csv_ouputs/wind_speed_{site}_fly_days.txt")

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
                        flyable = is_day_flyable(day_df, days)

                        ## check if day is flyable after analysis
                        if flyable == True:

                            ## add 1 to total
                            flyable_days += 1

                ## add total for month and year to dict
                flyable_days_dict[m_name].append(flyable_days)

            ## append to dict
            flyable_days_dict["lat"].append(lat)
            flyable_days_dict["lon"].append(lon)

        ## pickle out dictionary
        with open(f"csv_ouputs/wind_speed_{site}_fly_days.txt", "wb") as file:

            pickle.dump(flyable_days_dict, file, protocol=2)

    ## return flyable file
    return flyable_file


## checks if day is flyable against threshold
def is_day_flyable(day_df, days):

    ## define is day flyable
    is_day_flyable = False

    ## loopp through dataframe rows
    for index, row in day_df.iterrows():

        ## check if row is hour 23
        if row["Hour"] == 23:

            ## skip row as there's no pair --- change later to include following morning's hour
            continue

        ## otherwise
        else:

            ## get date string
            date_str = row["Date and Time"]

            ## convert to datetime
            date = datetime.strptime(date_str, "%Y-%m-%d %H:%M:%S")

            ## add 1 hour to date string
            consec_date = date + timedelta(hours=1)

            ## convert to date string
            consec_date_str = str(consec_date)

            ## get precip for each hour of pair
            wind_speed_1 = row["Wind Speed"]
            wind_speed_2 = day_df.loc[
                day_df["Date and Time"] == consec_date_str, "Wind Speed"
            ]

            ## check if both precip values are below threshold
            if (wind_speed_1 <= 10) and (wind_speed_2.item() <= 10):

                ## return true
                return True

            ## otherwise
            else:

                ## return false
                return False


if __name__ == "__main__":
    main()
print("finished")
