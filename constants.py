import numpy as np

## month dictionary
MONTHS = {"jan": [1, np.arange(1, 32, 1)], "feb": [2, np.arange(1, 29, 1)],
          "oct": [10, np.arange(1, 32, 1)], "nov": [11, np.arange(1, 31, 1)], 
          "dec": [12, np.arange(1, 32, 1)]}

M_FULL_NAMES = {"total": "Flyable Season", "jan": "January", "feb": "February", "oct": "October",
                "nov": "November", "dec": "December"}

## years list
YEARS = np.arange(1993, 2024, 1)

## site names
SITES = ["theron_hills", "provider", "argentina_range"]

## path to bas directory
BAS_PATH = "/home/users/lewis.davies/british_antarctic_survey/analysis_scripts"

SCR_PATH = "/data/scratch/lewis.davies/bas"

## site shapefile
SHP_FILE = f"{BAS_PATH}/bas_sites/bas_sites.shp"

## max possible days for each month type
MONTH_DAYS = {'total': 4681, 'jan': 961, 'feb': 868, 
              'oct': 961, 'nov': 930, 'dec': 958}

## weather variables
VARIABLES = ["Wind_Speed", "Gust", "cloud_base", "Visibility", "Precip"]

## limit combinations   
COMBOS = ["perfect", "twin_otter"]

SEASONS = ['summer', 'winter', 'flying']

## limits for individual files [twin_otter, perfect]
LIMITS = {'Precip': [0.1, 0.5], 'Wind_Speed': [20, 30, 50], 
               'cloud_base': [1000, 1800], 'Visibility': [1000, 5000, 10000], 
               'Gust': [30, 50]}
