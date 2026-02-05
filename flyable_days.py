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

import pickle
import warnings
from datetime import datetime, timedelta
from pathlib import Path
import numpy as np
import functions as func

import cartopy
import cartopy.crs as ccrs
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.collections import PatchCollection

from constants import (BAS_PATH, COMBOS, MONTH_DAYS, MONTHS, SITES, VARIABLES,
                       YEARS, SCR_PATH)

warnings.filterwarnings("ignore")

## csv file path
CSV_PATH = f"{BAS_PATH}/csv_ouputs"

## output directory for plots
OUTPUT_DIR = f"{BAS_PATH}/plots"


## runs other functions
def main():
    """ 
    Runs all functions within file.
    """
    
    ## makes heatmap plots
    make_heatmap_plots()

    return


## creates one dataframe from variable frames
def concat_data(cloud, precip, wind_spd, wind_gust, vis):
    """ 
    Concats weather variable dataframes into one dataframe.
    
    ARGS:
    
        cloud (DataFrame): cloud base data.
        precip (DataFrame): precipitation data.
        wind_spd (DataFrame): wind speed data.
        wind_gust (DataFrame): wind gust data.
        vis (DataFrame): visibility data.
        
    RETURNS:
    
        DataFrame: total concatenated dataframe
    """
    
    ## list of dataframes from iris files
    dataframes = [cloud, precip, wind_spd, wind_gust, vis]

    ## concat dataframes together
    total_df = pd.concat(dataframes, axis=1)

    ## remove duplicate columns
    total_df = total_df.loc[:, ~total_df.columns.duplicated()]

    ## drop na data
    total_df = total_df.dropna()

    ## return total dataframe
    return total_df


## counts 2 hour periods
def count_flyable_days(site, combo, data_df):
    """ 
    Counts flyable days per grid point for each Antarctica site.
    
    Args:
    
        site (str): Antarctica site.
        combo (str): specific combination of weather criteria.
        data_df (DataFrame): concatenated weather data
        
    Returns:

        .txt file:  pickled dictionary of flyable days
    """
    
    ## define flyable file path
    flyable_file = Path(f"csv_ouputs/all_vars_{combo}_{site}_fly_days.txt")

    ## check if file exists
    if flyable_file.is_file():

        print(flyable_file, 'exists!')
        ## if exists, return file
        return flyable_file

    else:

        if combo == "perfect":

            ## define upper thresholds
            wind_gust_limit = 30
            wind_speed_limit = 20
            vis_limit = 10000
            cld_base_limit = 1800
            precip_limit = 0.5

        else:

            ## define Twin Otter thresholds
            wind_gust_limit = 30
            wind_speed_limit = 20
            vis_limit = 5000
            cld_base_limit = 1000
            precip_limit = 0.5

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

                days = m_number[1]

                ## constrain by month
                month_df = lat_lon_df.loc[data_df["Month"] == m_number[0]]

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
                        flyable = is_day_flyable(day_df, limits)

                        ## check if day is flyable after analysis
                        if flyable:

                            ## add 1 to total
                            flyable_days += 1

                ## add total for month and year to dict
                flyable_days_dict[m_name].append(flyable_days)

            ## append to dict
            flyable_days_dict["lat"].append(lat)
            flyable_days_dict["lon"].append(lon)

        ## pickle out dictionary
        with open(f"csv_ouputs/all_vars_{combo}_{site}_fly_days.txt", "wb") as file:

            pickle.dump(flyable_days_dict, file, protocol=2)

    ## return flyable file
    return flyable_file


## check if day is flyable
def is_day_flyable(day_df, limits):
    """ 
    Checks all weather variables against criteria for 2 hour windows 
    where all are met.
    
    Args:
    
        day_df (DataFrame): weather data for 1 day
        limits (dict): weather parameter criteria.
        
    Returns:
    
        True or False
    """
    
    day_df.replace('', np.nan, inplace=True)
    day_df.dropna(inplace=True)
    
    
    ## loopp through dataframe rows
    for index, row in day_df.iterrows():

        ## check if row is hour 23
        if (row["Hour"] == 23) | ((row["Hour"] == 20) & (row['Month'] == 2)):

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
                if var in ["Wind_Speed", "Gust", "Precip"]:
                    
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

            if all(x == 1 for x in all_limits_met):

                return True


    ## day is not flyable, return false
    return False


## runs functions needed for making heatmaps
def make_heatmap_plots():
    """ 
    Runs concat_data, count_flyable_days, is_day_flyable to make_site_heatmap & 
    make_area_heatmap.
    
    """
    
    big_dict = {
            "site": [],
            "combo": [],
            "lat": [],
            "lon": [],
            "jan": [],
            "feb": [],
            "oct": [],
            "nov": [],
            "dec": [],
        }
    
    ## loop through combinations
    for combo in COMBOS:

        ## loop through sites
        for site in SITES:

            ## lines 339 - 361: comment out after *fly_days.txt files made

            ## csv file paths
            cloud_file = f"{SCR_PATH}/cloud_base_{site}_stats_ml.csv"
            precip_file = f"{SCR_PATH}/Precip_{site}_stats.csv"
            wind_spd_file = f"{SCR_PATH}/Wind_Speed_{site}_stats.csv"
            wind_gust_file = f"{SCR_PATH}/Gust_{site}_stats.csv"
            vis_file = f"{SCR_PATH}/Visibility_{site}_stats.csv"

            ## load csv files
            cloud_df = pd.read_csv(cloud_file, index_col=0)
            precip_df = pd.read_csv(precip_file, index_col=0)
            wind_spd_df = pd.read_csv(wind_spd_file, index_col=0)
            wind_gust_df = pd.read_csv(wind_gust_file, index_col=0)
            vis_df = pd.read_csv(vis_file, index_col=0)

            ## concat to single dataframe
            total_df = concat_data(
                cloud_df, precip_df, wind_spd_df, wind_gust_df, vis_df
            )

            total_df = func.create_flying_season(total_df)
            
            ## count flyable days and return file
            flyable_file = count_flyable_days(site, combo, total_df)
            
            flyable_file = Path(f"{SCR_PATH}/csv_ouputs/csv_ouputs/all_vars_{combo}_{site}_fly_days.txt")

            ## open flyable daya dict file
            with open(flyable_file, "rb") as file:

                flyable_days_dict = pickle.load(file)

            ## turn into dataframe
            flyable_df = pd.DataFrame(flyable_days_dict)

            ## loop through dataframe rows
            for index, row in flyable_df.iterrows():

                ## loop through dicitonay keys
                for var in list(big_dict.keys())[2:]:

                    ## append value from row to dictionary for area heatmap
                    big_dict[f"{var}"].append(row[f"{var}"].item())

                big_dict['combo'].append(combo)
                big_dict['site'].append(site)
            

            ## make heatmap for sites based on months and totals
            ## make_site_heatmap(flyable_df, site, combo)
            ## make_site_heatmap(flyable_df, site, combo)

    ## create dataframe from area dictionary
    area_df = pd.DataFrame(big_dict)
        
    ## make heatmap of whole area
    make_area_heatmap(area_df)


## make heatmap of all 3 sites together
def make_area_heatmap(area_df):
    """ 
    Makes heatmap of all 3 sites on one plot for each month and all
    months together.
    
    Args:
    
        area_df (DataFrame): dataframe of flyable days by grid point.
        combo (str):  specific combination of weather criteria.
    """

    ## create total column
    area_df["total"] = (
            area_df["jan"]
            + area_df["feb"]
            + area_df["oct"]
            + area_df["nov"]
            + area_df["dec"])


    area_df["lon"] = area_df["lon"].apply(lambda x: x-360)

    ## get longitude and latitudes
    lons, lats = area_df["lon"], area_df["lat"]

    ## set plot aspect
    lon_width = max(lons) - min(lons)
    lat_height = max(lats) - min(lats)
    plot_aspect = lat_height / lon_width
    
    ## create normalised data column
    area_df["norm"] = (area_df["total"] - area_df["total"].min()) / (
        area_df["total"].max() - area_df["total"].min()
    )
    
    cmap = cm.get_cmap("viridis")
    
    flyable_data = (area_df["total"] / MONTH_DAYS["total"])*100
    
    ## create scalar mappable for colourbar
    norm = plt.Normalize(flyable_data.min(), flyable_data.max())
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

    for combo in COMBOS:
        
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
        
        flyable_df = area_df[area_df['combo'] == combo]
        
        flyable_norm = flyable_df['norm']
        
        ## create patch collection for making heatmap
        patches = [
            plt.Rectangle(
                (lon - 0.125, lat - 0.125), 0.25, 0.25, 
                alpha=1, facecolor=cmap(days)
            )
            for lon, lat, days in zip(lons, lats, flyable_norm)
        ]
        
        patch_collection = PatchCollection(patches, match_original=True)
        
        ## add heatmap patch collection to axes
        ax.add_collection(patch_collection)
        
        ## add title
        ax.set_title("Flyable Days by Grid Point for Flyable Season",
                     fontsize=18,
                     fontweight="bold")

        gls = ax.gridlines(draw_labels = True, alpha = 0.2)
        
        gls.top_labels=False   # suppress top labels
        gls.right_labels=False # suppress right labels
        
        ## add colourbar
        plt.colorbar(
            sm,
            ax=ax,
            label="Flyable Days (%)",
            ticks=[
                80, 85, 86, 87, 88, 89, 90, 91, 92, 93, 
                94, 95, 96, 97, 98, 99, 100,
            ],
        )

        ## save figure
        fig.savefig(
            f"{BAS_PATH}/plots/{combo}_all_sites_all_vars_heatmap.png",
            bbox_inches = 'tight'
        )
    

## make heatmap for each site individually
def make_site_heatmap(flyable_df, site, combo):
    """ 
    Makes heatmap of all individual sites for each month and all months.
    
    Args:
    
        flyable_df (DataFrame): flyable days by grid point and site.
        site (str): Antarctica site.
        combo (str): specific combination of weather criteria.
    """
    
    # month list for creating plots
    months = ["total", "jan", "feb", "oct", "nov", "dec"]

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
        flyable_data = (flyable_df[f"{month}"] / MONTH_DAYS[f'{month}'])*100

        ## create normalised data column
        flyable_df["norm"] = (flyable_df[f"{month}"] - flyable_df[f"{month}"].min()) / (
            flyable_df[f"{month}"].max() - flyable_df[f"{month}"].min()
        )
        
        flyable_norm = flyable_df["norm"]

        ## get colour map
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
            f"{BAS_PATH}/plots/{month}_{site}_{combo}__all_vars_heatmap.png"
        )


if __name__ == "__main__":
    main()
