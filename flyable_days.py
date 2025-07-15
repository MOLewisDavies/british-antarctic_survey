""" 

script to calculate number of flyable days from all BAS criteria.


"""

## package imports
import matplotlib.pyplot as plt 
import pandas as pd
import numpy as np
import seaborn as sns

## set seaborn style
sns.set_style('darkgrid')

## csv file path
CSV_PATH = '/home/users/lewis.davies/british_antarctic_survey/analysis_scripts/csv_ouputs'

OUTPUT_DIR = "/home/users/lewis.davies/british_antarctic_survey/analysis_scripts/plots"

MONTHS = ['jan', 'feb', 'oct', 'nov', 'dec']

YEARS = np.arange(1993, 2024, 1)

VARIABLES = ['wind_speed', 'wind_gust', 'cld_base', 
             'vis', 'total_precip']

## site names
SITES = ["theron_hills", "provider", "argentina_range"]

LIMITS = ['lower', 'upper']


## runs other functions
def main():
    
    for site in SITES:

        ## csv files
        CLOUD_FILE = f'{CSV_PATH}/cloud_base_{site}_stats.csv'
        PRECIP_FILE = f'{CSV_PATH}/precip_{site}_stats.csv'
        WIND_SPD_FILE = f'{CSV_PATH}/wind_speed_{site}_stats.csv'
        WIND_GUST_FILE = f'{CSV_PATH}/wind_gust_{site}_stats.csv'
        VIS_FILE = f'{CSV_PATH}/vis_{site}_stats.csv'

        total_df = concat_data(CLOUD_FILE, PRECIP_FILE, 
                               WIND_SPD_FILE, WIND_GUST_FILE, 
                               VIS_FILE)
        
        for limit in LIMITS:
            
            flyable_days_dict = count_flyable_windows(limit, total_df, site)

            make_line_plot(flyable_days_dict, site)

    return

## matches each data variable to datetime for comparison
def concat_data(CLOUD_FILE, PRECIP_FILE, WIND_SPD_FILE, WIND_GUST_FILE, VIS_FILE):
    
    cloud_df = pd.DataFrame(CLOUD_FILE)
    w_spd_df = pd.DataFrame(WIND_SPD_FILE)
    w_gust_df = pd.DataFrame(WIND_GUST_FILE)
    precip_df = pd.DataFrame(PRECIP_FILE)
    vis_df = pd.DataFrame(VIS_FILE)

    total_df = pd.merge(cloud_df, w_gust_df, 
                        w_spd_df, precip_df, vis_df, 
                        on=['full_date', 'month', 'day', 'hour'])
    
    total_df = total_df.dropna()

    return total_df

## counts 2 hour periods
def count_flyable_windows(limit, total_df):

    if limit = 'upper':
        
        ## define upper thresholds 
        wind_gust_limit = 36
        wind_speed_limit = 50
        vis_limit = 5000
        cld_base_limit = 1800
        precip_limit = 0.0001

    else: 

        ## define lower thresholds 
        wind_gust_limit = 27
        wind_speed_limit = 36
        vis_limit = 1000
        cld_base_limit = 1800
        precip_limit = 0.0001


    ## list thresholds
    limits = [wind_speed_limit, wind_gust_limit, 
              cld_base_limit, vis_limit, precip_limit]

    ## empty dict for collecting month: flyable day counts
    flyable_days_dict = {'jan': [], 
                         'feb': [],
                         'oct': [],
                         'nov': [],
                         'dec': []}

    ## loop through years
    for year in YEARS:

        ## loop through months
        for month in MONTHS:
            
            ## flyable day total
            flyable_days = 0 

            ## get correct number of days
            if month in ('jan', 'oct', 'dec'):

                days = np.arange(1, 32, 1)

            elif month == 'feb':

                days = np.arange(1, 29, 1)

            else:

                days = np.arange(1, 31, 1)
            
            ## define dataframe constraints
            year_con = (total_df[year] == year)
            month_con = (total_df[year] == month)

            ## subset dataframe for year and month
            year_month_df = total_df.loc[year_con & month_con]

            ## loop through days in month
            for day in days:
                
                ## subset by day
                day_df = year_month_df.loc[year_month_df[day] == day]

                ## get True or False for a flyable day
                flyable = is_day_flyable(total_df, day_df, limits)
                
                ## check if day is flyable after analysis
                if flyable == True:
                    
                    ## add 1 to total
                    flyable_days =+ 1

                ## if not, move to next day
                else:

                    continue
            
            ## add total for month and year to dict
            flyable_days_dict[month].append(flyable_days)
    
    return flyable_days_dict

## check if day is flyable
def is_day_flyable(total_df, day_df, limits):

    ## define day as not flyable
    is_day_flyable = False

    ## loop through df rows
    for row in day_df.iterrows():
        

        date = row['full_date']
        consec_date = (date + timedelta(hours = 1))
        grid = row['grid_square']

        ## get every unique consecutive 2 hour window for each grid square
        hour_1 = row
        hour_2 = total_df.iloc[[total_df['full_date'] == consec_date] and [total_df['grid_square'] == grid]]

        ## define all limits as being unmet
        all_limits_met = False

        ## loop through variables and corresponding limit
        for var, limit in zip(VARIABLES, limits):
                        
            ## if speed, gust or precip limit is exceeded if greater than
            if var in ('wind_speed', 'wind_gust', 'precip'):
                            
                ## checks if either hour breaks limit
                if (hour_1[var] > limit) or (hour_2[var] > limit):
                                
                    ## checks if either hour breaks limit
                    all_limits_met = True

                else:
                            
                    ## checks if either hour breaks limit
                    if (hour_1[var] < limit) or (hour_2[var] < limit):
                                
                        ## checks if either hour breaks limit
                        all_limits_met = True

                        ## check if limits met
                if all_limits_met == True:
                            
                    ## add 1 counter
                    is_day_flyable = True

                    ## day is flyable, return True
                    return is_day_flyable

                else:
                            
                    ## continue if not flyable yet
                    continue
    
    ## day is not flyable, return false
    return is_day_flyable


def make_box_plot(data):

    data_df = pd.DataFrame(data)
    
    fig, ax = plt.subplots(figsize=(10, 6))

    box_plot = sns.boxplot(data=data_df)

    ##### function needs finishing #####


    

    return


def make_line_plot(data, site):

    data_df = pd.DataFrame.from_dict(data.items(), columns = ['months', 'flyable_days'])
    data_df = data_df.reset_index()

    fig, ax = plt.subplots()
    
    sns.lineplot(data=data, x = data_df['months'], y = data_df['flyable_days'])

    ax.set_xlabel('Months')
    ax.set_ylabel('No. of flyable days')
    ax.set_title('No. of Flyable Days Per Month (any grid point)')
    
    plt.tight_layout()
    fig.savefig(f'{OUTPUT_DIR}/flyable_days_all_grid_points_{site}.png')
    plt.close()
    
    
def make_heatmap():
    return