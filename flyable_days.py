"""

script to calculate number of flyable days from all BAS criteria.

functions in file():

main(): runs functions
concat_data(): combines individual dataframes
count_flyable_days(): counts flyable days and saves as txt file
is_day_flyable(): determines if day has flyable window
make_heatmap_plots(): runs functions for making heatmaps
make_area_heatmap(): makes heatmap of all 3 sites
make_site_heatmap(): makes heatmaps of each site individually


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

## set seaborn style
sns.set_style("darkgrid")

## csv file path
CSV_PATH = (
    "/home/users/lewis.davies/british_antarctic_survey/analysis_scripts/csv_ouputs"
)

OUTPUT_DIR = "/home/users/lewis.davies/british_antarctic_survey/analysis_scripts/plots"

MONTHS = {"jan": 1, "feb": 2, "oct": 10, "nov": 11, "dec": 12}

YEARS = np.arange(1993, 2024, 1)

VARIABLES = ["Wind Speed", "Gust", "cloud_base", "Visibility", "Precip"]

## site names
SITES = ["theron_hills", "provider", "argentina_range"]

LIMITS = ["twin_otter"]


## runs other functions
def main():

    ## makes heatmap plots
    make_heatmap_plots()

    return


## creates one dataframe from variable frames
def concat_data(CLOUD, PRECIP, WIND_SPD, WIND_GUST, VIS):

    dataframes = [CLOUD, PRECIP, WIND_SPD, WIND_GUST, VIS]

    total_df = pd.concat(dataframes, axis=1)

    total_df = total_df.loc[:, ~total_df.columns.duplicated()]

    print(total_df.columns)

    total_df = total_df.dropna()

    return total_df


## counts 2 hour periods
def count_flyable_days(site, limit, data_df):

    ## define flyable file path
    flyable_file = Path(f"csv_ouputs/all_vars_{limit}_{site}_fly_days.txt")

    ## check if file exists
    if flyable_file.is_file():

        ## if exists, return file
        return flyable_file

    else:

        if limit == "perfect":

            ## define upper thresholds
            wind_gust_limit = 36
            wind_speed_limit = 50
            vis_limit = 5000
            cld_base_limit = 1800
            precip_limit = 0.0001

        else:

            ## define lower thresholds
            wind_gust_limit = 10
            wind_speed_limit = 10
            vis_limit = 5000
            cld_base_limit = 1800
            precip_limit = 0.0005

        ## list thresholds
        limits = [
            wind_speed_limit,
            wind_gust_limit,
            cld_base_limit,
            vis_limit,
            precip_limit,
        ]

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
                        flyable = is_day_flyable(data_df, day_df, days)

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
        with open(f"csv_ouputs/all_vars_{limit}_{site}_fly_days.txt", "wb") as file:

            pickle.dump(flyable_days_dict, file, protocol=2)

    ## return flyable file
    return flyable_file


## check if day is flyable
def is_day_flyable(total_df, day_df, limits):

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
            hour_1 = row
            hour_2 = day_df.loc[day_df["Date and Time"] == consec_date_str]

            ## define all limits as being unmet
            all_limits_met = []

            ## loop through variables and corresponding limit
            for var, limit in zip(VARIABLES, limits):

                ## if speed, gust or precip limit is exceeded if greater than
                if var in ["Wind Speed", "Gust", "Precip"]:

                    ## checks if either hour breaks limit
                    if (hour_1[var] <= limit) and (hour_2[var].item() <= limit):

                        ## checks if either hour breaks limit
                        all_limits_met.append(1)

                    else:

                        all_limits_met.append(0)

                else:

                    ## checks if either hour breaks limit
                    if (hour_1[var] > limit) and (hour_2[var].item() > limit):

                        ## checks if either hour breaks limit
                        ## checks if either hour breaks limit
                        all_limits_met.append(1)

                    else:

                        all_limits_met.append(0)

            # print(all_limits_met)

            if all(x == 1 for x in all_limits_met) == True:

                return True

            else:

                continue

    ## day is not flyable, return false
    return False


## runs functions needed for making heatmaps
def make_heatmap_plots():

    for limit in LIMITS:

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

        for site in SITES:

            ## csv file paths
            CLOUD_FILE = f"{CSV_PATH}/cloud_base_{site}_stats.csv"
            PRECIP_FILE = f"{CSV_PATH}/precip_{site}_stats.csv"
            WIND_SPD_FILE = f"{CSV_PATH}/wind_speed_{site}_stats.csv"
            WIND_GUST_FILE = f"{CSV_PATH}/wind_gust_{site}_stats.csv"
            VIS_FILE = f"{CSV_PATH}/visibility_{site}_stats.csv"

            ## load csv files
            cloud_df = pd.read_csv(CLOUD_FILE, index_col=0)
            precip_df = pd.read_csv(PRECIP_FILE, index_col=0)
            wind_spd_df = pd.read_csv(WIND_SPD_FILE, index_col=0)
            wind_gust_df = pd.read_csv(WIND_GUST_FILE, index_col=0)
            vis_df = pd.read_csv(VIS_FILE, index_col=0)

            ## concat to single dataframe
            total_df = concat_data(
                cloud_df, precip_df, wind_spd_df, wind_gust_df, vis_df
            )

            ## count flyable days and return file
            flyable_file = count_flyable_days(site, limit, total_df)

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

            print(flyable_df)

            ## make heatmap for sites based on months and totals
            make_site_heatmap(flyable_df, "month", site, limit)
            make_site_heatmap(flyable_df, "total", site, limit)

        ## create dataframe from area dictionary
        area_df = pd.DataFrame(area_dict)

        ## make heatmap of whole area
        make_area_heatmap(area_df, limit)


## make heatmap of all 3 sites together
def make_area_heatmap(area_df, limit):

    months = ["toal"]

    days = (
        max(area_df["jan"])
        + max(area_df["feb"])
        + max(area_df["oct"])
        + max(area_df["nov"])
        + max(area_df["dec"])
    )
    print(days)

    ## create total column
    area_df["total"] = (
        (
            area_df["jan"]
            + area_df["feb"]
            + area_df["oct"]
            + area_df["nov"]
            + area_df["dec"]
        )
        / 4681
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

    flyable_data = area_df["total"]

    ## create normalised data column
    area_df["norm"] = (area_df["total"] - area_df["total"].min()) / (
        area_df["total"].max() - area_df["total"].min()
    )
    flyable_norm = area_df["norm"]

    cmap = cm.get_cmap("viridis")

    ## create patch collection for making heatmap
    patches = [
        plt.Rectangle(
            (lon - 0.125, lat - 0.125), 0.25, 0.25, alpha=1, facecolor=cmap(days)
        )
        for lon, lat, days in zip(lons, lats, flyable_norm)
    ]

    patch_collection = PatchCollection(patches, match_original=True)

    ## create scalar mappable for colourbar
    norm = plt.Normalize(flyable_data.min(), flyable_data.max())
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
            flyable_data.min(),
            flyable_data.median(),
            flyable_data.max(),
        ],
    )

    ## save figure
    fig.savefig(
        f"/home/users/lewis.davies/british_antarctic_survey/analysis_scripts/plots/total_all_sites_all_vars_heatmap.png"
    )

    return


## make heatmap for each site individually
def make_site_heatmap(flyable_df, d_type, site, limit):

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
        flyable_data = flyable_df[f"{month}"]

        ## create normalised data column
        flyable_df["norm"] = (flyable_df[f"{month}"] - flyable_df[f"{month}"].min()) / (
            flyable_df[f"{month}"].max() - flyable_df[f"{month}"].min()
        )
        flyable_norm = flyable_df["norm"]

        ## get colour map
        cmap = cm.get_cmap("viridis")

        ## add gridpoint crosses to plot
        # sns.scatterplot(data = flyable_df, x=lons, y=lats, color='red', size=200, marker='+', ax=ax)

        ## create patch collection for making heatmap
        patches = [
            plt.Rectangle(
                (lon - 0.125, lat - 0.125), 0.25, 0.25, alpha=1, facecolor=cmap(days)
            )
            for lon, lat, days in zip(lons, lats, flyable_norm)
        ]
        patch_collection = PatchCollection(patches, match_original=True)

        ## create scalar mappable for colourbar
        norm = plt.Normalize(flyable_data.min(), flyable_data.max())
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
                flyable_data.min(),
                flyable_data.median(),
                flyable_data.max(),
            ],
        )

        ## save figure
        fig.savefig(
            f"/home/users/lewis.davies/british_antarctic_survey/analysis_scripts/plots/{month}_{site}_all_vars_heatmap.png"
        )


if __name__ == "__main__":
    main()
