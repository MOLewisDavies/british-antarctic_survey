import numpy as np

## month dictionary
MONTHS = {"jan": 1, "feb": 2, "oct": 10, "nov": 11, "dec": 12}

## years list
YEARS = np.arange(1993, 2024, 1)

## site names
SITES = ["theron_hills", "provider", "argentina_range"]

## path to bas directory
BAS_PATH = "/home/users/lewis.davies/british_antarctic_survey/analysis_scripts"

## site shapefile
SHP_FILE = f"{BAS_PATH}/bas_sites/bas_sites.shp"

## max possible days for each month type
MONTH_DAYS = {'total': 4681, 'jan': 961, 'feb': 868, 
              'oct': 961, 'nov': 930, 'dec': 961}

## weather variables
VARIABLES = ["Wind_Speed", "Gust", "cloud_base", "Visibility", "Precip"]

## limit combonations
COMBOS = ["perfect", "twin_otter"]
