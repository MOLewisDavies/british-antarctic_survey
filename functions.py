""" 

functions used by programmes in BAS project
import into specific files

IMPORTANT: Flyable Days file uses slightly different make_heatmap_plots
           so these are for separate parameters only.

"""

import pickle
from datetime import datetime, timedelta
from pathlib import Path
import numpy as np

import cartopy
import cartopy.crs as ccrs
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import pandas as pd
from cartopy.io.shapereader import Reader
from iris.util import mask_cube_from_shapefile
from matplotlib.collections import PatchCollection
import seaborn as sns

from constants import (MONTH_DAYS, MONTHS, SITES, YEARS, BAS_PATH, 
                       LIMITS)


## make plots
def make_heatmap_plots(combo, var):

    """ 
    Runs count_flyable_days, to make_site_heatmap & make_area_heatmap.
    
    Args:

        combo (str): specific combination of weather parameters.
        var (str): weather parameter.
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

    ## loop through SITES
    for site in SITES:

        ## load dataframe
        data_df = pd.read_csv(f"csv_ouputs/{var}_{site}_stats.csv")

        data_df = create_flying_season(data_df)
        
        ## create flyable file
        flyable_file = count_flyable_days(combo, data_df, site, var)

        ## open flyable daya dict file
        with open(flyable_file, "rb") as file:

            flyable_days_dict = pickle.load(file)

        ## turn into dataframe
        flyable_df = pd.DataFrame(flyable_days_dict)

        ## loop through dataframe rows
        for index, row in flyable_df.iterrows():

            ## loop through dicitonay keys
            for key in area_dict.keys():

                ## append value from row to dictionary for area heatmap
                area_dict[f"{key}"].append(row[f"{key}"].item())

        ## make heatmap for sites based on months and totals
        make_site_heatmap(flyable_df, site, combo, var)

    ## create dataframe from area dictionary
    area_df = pd.DataFrame(area_dict)

    ## make heatmap of whole area
    make_area_heatmap(area_df, combo, var)
    

## make heatmap of all sites
def make_area_heatmap(area_df, combo, var):

    """ 
    Makes heatmap of all 3 sites on one plot for each month and all
    months.
    
    Args:
    
        area_df (DataFrame): dataframe of flyable days by grid point.
        combo (str): specific combination of weather parameters.
        var (str): weather parameter.
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
            [min(lons)-0.5, max(lons) + 0.5, min(lats) - 0.5, max(lats) + 0.5],
            ccrs.PlateCarree(),
        )

        # Draw coastlines and background stuff
        ax.add_feature(cartopy.feature.LAND)
        ax.add_feature(cartopy.feature.OCEAN)
        ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.2)
        ax.add_feature(cartopy.feature.BORDERS, linestyle=":", linewidth=0.3)
        ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
        ax.add_feature(cartopy.feature.RIVERS)

        ## define data + convert to percentage
        flyable_data = (area_df[f"{month}"] / MONTH_DAYS[f'{month}'])*100

        ## create normalised data column
        area_df["norm"] = (area_df[f"{month}"] - area_df[f"{month}"].min()) / (
            area_df[f"{month}"].max() - area_df[f"{month}"].min()
        )
        flyable_norm = area_df["norm"]

        ## get colour map
        cmap = cm.get_cmap("viridis")


        ## create patch collection for making heatmap
        patches = [
            plt.Rectangle(
                (lon - 0.125, lat - 0.125), 
                0.25, 0.25, alpha=1, facecolor=cmap(val)
            )
            for lon, lat, val in zip(lons, lats, flyable_norm)
        ]
        patch_collection = PatchCollection(patches, match_original=True)

        ## create scalar mappable for colourbar
        norm = plt.Normalize(flyable_data.min(), flyable_data.max())
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

        ## add heatmap patch collection to axes
        ax.add_collection(patch_collection)

        ## add title
        ax.set_title(f"{month.capitalize()} Flyable Days by Grid: All Sites")

        ## add colourbar
        plt.colorbar(
            sm,
            ax=ax,
            label="Flyable Days (%)",
            ticks=[50, 
                   flyable_data.median(), 
                   100]
        )
        
        ## save figure
        fig.savefig(
            f"{BAS_PATH}/plots/{var}_{combo}_{month}_all_sites_heatmap.png",
            bbox_inches = 'tight'
        )
        
        
## make heatmap for each site individually
def make_site_heatmap(flyable_df, site, combo, var):

    """ 
    Makes heatmap of all individual sites for each month and all months.
    
    Args:
    
        flyable_df (DataFrame): flyable days by grid point and site.
        site (str): Antarctica site.
        combo (str): specific combination of weather parameters
        var (str): weather parameter.
    """

    # month list for creating plots
    months = ["total", "jan", "feb", "oct", "nov", "dec"]

    ## create total column
    flyable_df["total"] = (
        flyable_df["jan"]
            + flyable_df["feb"]
            + flyable_df["oct"]
            + flyable_df["nov"]
            + flyable_df["dec"])

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

        ## define data + convert to percentage
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
                (lon - 0.125, lat - 0.125), 0.25, 0.25, 
                alpha=1, 
                facecolor=cmap(val)
            )
            for lon, lat, val in zip(lons, lats, flyable_norm)
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
            ticks=[flyable_data.min(), 
                   flyable_data.median(), 
                   flyable_data.max()],
        )

        ## save figure
        fig.savefig(
            f"{BAS_PATH}/plots/{var}_{combo}_{month}_{site}_heatmap.png"
        )
        
        
## count flyable days for heatmap plots
def count_flyable_days(combo, data_df, site, var):
    
    """ 
    Counts flyable days per grid point for each Antarctica site.
    
    Args:

        combo (str): specific combination of weather parameters
        data_df (DataFrame): weather parameter stats.
        site (str): Antarctica site.
        var (str): weather parameter.
        
    Returns:

        .txt file:  pickled dictionary of flyable days
    """
    
    ## define flyable file path
    flyable_file = Path(f"csv_ouputs/{var}_{combo}_{site}_fly_days.txt")

    ## get list of unique lat lon pairs
    lat_lon_list = data_df[["Latitude", "Longitude"]].drop_duplicates()
    
    ## convert into individual lat lon lists
    lats = lat_lon_list["Latitude"].tolist()
    lons = lat_lon_list["Longitude"].tolist()
    
    ## check if file exists
    if flyable_file.is_file():

        print(flyable_file, 'exists')
        
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
                        flyable = is_day_flyable(day_df, combo, var)

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
        with open(f"csv_ouputs/{var}_{combo}_{site}_fly_days.txt", "wb") as file:

            pickle.dump(flyable_days_dict, file, protocol=2)

    ## return flyable file
    return flyable_file


## checks if day is flyable against threshold
def is_day_flyable(day_df, combo, var):

    """ 
    Checks weather data against criteria for 2 hour windows.
    
    Args:
    
        day_df (DataFrame): weather data for 1 day
        var (str): weather parameter
        combo (str): specific combination of weather parameters
        
    Returns:
    
        True or False
    """
    
    if combo == 'perfect':
    
        limit = LIMITS[f'{var}'][1]
        
    else:
        
        limit = LIMITS[f'{var}'][0]
    
    
    ## loopp through dataframe rows
    for index, row in day_df.iterrows():

        ## check if row is hour 20 - utc-3 time
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
            consec_str = str(consec_date)

            ## get data for each hour of pair
            hour_1 = row[f"{var}"]
            hour_2 = day_df.loc[day_df["Date and Time"] == consec_str, f"{var}"]

            if var in ["Wind_Speed", "Gust", "Precip"]:
                
                ## check if both values are below threshold
                if (hour_1 <= limit) and (hour_2.item() <= limit):

                    ## return true
                    return True

            else:

                ## checks if either hour breaks limit
                if (hour_1 > limit) and (hour_2.item() > limit):
                    
                    ## return true
                    return True
               
                
## mask cube using shapefile polygons
def mask_cube(index, cube, shape_file):

    """ 
    Masks iris cube based on shape file.
    
    Args:
    
        index (int): index of shape file in list
        cube (iris cube): weather cube to be masked
        shape_file (path): path to shape file
        
    Returns:
    
        iris cube: masked weather data cube
    """

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


## creates thirty whole flying seasons from 31 years of data
def create_flying_season(data_df):
    """ 
    Drops Dec from 92, Jan, Feb, Sept from 93 and Oct, Nov, Dec from 23 to 
    create 30 full flying seasons.
    
    Args:
    
        data_df (DataFrame): weather variable data.
        
    Returns:
    
        DataFrame: dataframe with rows dropped.
    """
    
    
    ## conditions for dropping each month
    
    dec92_con = (data_df['Year'] == 1992) & (data_df['Month'] == 12)
    sep93_con = (data_df['Year'] == 1993) & (data_df['Month'] == 9)
    jan_con = (data_df['Year'] == 1993) & (data_df['Month'] == 1)
    feb_con = (data_df['Year'] == 1993) & (data_df['Month'] == 2)
    oct_con = (data_df['Year'] == 2023) & (data_df['Month'] == 10)
    nov_con = (data_df['Year'] == 2023) & (data_df['Month'] == 11)
    dec_con = (data_df['Year'] == 2023) & (data_df['Month'] == 12)
    
    ## list as a set of "or" conditions
    cons = (jan_con|feb_con|oct_con|nov_con|dec_con|dec92_con|sep93_con)
    
    ## drop rows from dataframe
    data_df = data_df.drop(data_df[cons].index)
    
    ## return dataframe
    return data_df


## make box plots
def make_box_plot(var, d_type, grid_points):
    
    """ 
    Creates box plot using weather data basedon grid point and analysis 
    method.
    
    ARGS:
    
        var (str): Weather variable.
        d_type (str): Analysis method (takes EVERY or AVERAGE)
        grid_points (str): No. of gtrid points (takes ALL or TOP THREE)
    
    """
    
    ## set up axes
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ## empty dataframe list
    df_list = []
    
    ## loop through sites
    for site in SITES:
        
        ## load dataframe
        data_df = pd.read_csv(f"csv_ouputs/{var}_{site}_stats.csv")
        
        ## add site column
        data_df['site'] = site
        
        ## group dataframe by year and month to subset dataframe
        for data, grouped_df in data_df.groupby(['Year', 'Month']):
            
            ## get month name from dict using month number in dataframe
            m_name = [key for key, val in MONTHS.items() if val[0] == data[1]]
            
            ## change month numbers to names in dataframe 
            grouped_df.loc[grouped_df['Month'] == data[1], 'Month'] = m_name[0]
            
            ## if grid points is all
            if grid_points == 'ALL':
                
                ## checks analysis method
                if d_type == "EVERY":
                    
                    ## append dataframe to list
                    df_list.append(grouped_df)
                    
                ## checks analysis method    
                elif d_type == 'AVERAGE':
                    
                    ## set up average dictionary
                    average_dict = {
                        "Year": [],
                        "Month": [],
                        f"{var}": [],
                        "site": [],
                    }
                    
                    ## get weather values
                    values = grouped_df[f'{var}']
                    
                    ## append data to average dictionary
                    average_dict["Year"].append(data[0])
                    average_dict["Month"].append(m_name[0])
                    average_dict[f"{var}"].append(np.mean(values))
                    average_dict["site"].append(site)
                    
                    ## create dataframe from dict
                    average_df = pd.DataFrame(average_dict)

                    ## append dictionary to list
                    df_list.append(average_df)
            
            ## if grid points is top three
            elif grid_points == "TOP THREE":
                
                ## gets dataframe of 3 largest values
                largest_df = grouped_df.nlargest(3, f'{var}')
                
                ## checks analysis method
                if d_type == "EVERY":
                    
                    ## append to list
                    df_list.append(largest_df)
                    
                ## checks analysis method
                elif d_type == "AVERAGE":
                    
                    ## get weather values
                    values = largest_df[f'{var}']
                    
                    ## set up average dictionary
                    average_dict = {
                        "Year": [],
                        "Month": [],
                        f"{var}": [],
                        "site": [],
                    }
                    
                     ## append data to average dictionary
                    average_dict["Year"].append(data[0])
                    average_dict["Month"].append(m_name[0])
                    average_dict[f"{var}"].append(np.mean(values))
                    average_dict["site"].append(site)
                    
                    ## create dataframe from dict
                    average_df = pd.DataFrame(average_dict)

                    ## append dictionary to list
                    df_list.append(average_df)
                          
    ## concat dataframe list
    area_df = pd.concat(df_list)

    
    ## create box plot
    sns.boxplot(data=area_df, x="Month", y=var, 
                order = ['oct', 'nov', 'dec', 'jan', 'feb'], hue="site",
                showfliers = True, ax=ax)

    ##  add limits to box plots as lines
    ax.axhline(y=LIMITS[f'{var}'][0], ls='dashed', color='r', alpha=0.5,
               label = 'Twin Otter')
    
    ax.axhline(y=LIMITS[f'{var}'][1], ls='dashed', color='g', alpha=0.5,
               label = 'Perfect')
    
    ## set title
    ax.set_title(
        f"{var} values - {grid_points} grid point(s) - {d_type}",
        fontsize=18,
        fontweight="bold",
    )

    ## save figure
    fig.savefig(
        f"{BAS_PATH}/plots/monhtly_{var}_whisker_plot_{grid_points}_{d_type}.png",
        bbox_inches = 'tight'
    )
                    
                    
                
                
            

        
    
