import pathlib
import pandas as pd
import geopandas
import rioxarray
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from calendar import monthrange,month_name
from datetime import datetime
import argparse
from itertools import product
from tqdm import tqdm


#--------------------------------------------------------------------
def day_range(year=2024,month=1,day=None):
   """
   Determine the yearly day number for the first and last days
   of the selected month. Take into account leap years.
   """
   if(day is None):
      last_day  = monthrange(year, month)[1]

      first_day = datetime(year, month, 1).timetuple().tm_yday
      last_day  = datetime(year, month, last_day).timetuple().tm_yday

      return first_day,last_day
   else:
      return (datetime(year, month, day).timetuple().tm_yday,)

#--------------------------------------------------------------------
def read_data(quantity="mean temperature",year=2024):
   # Load the raster and reproject it to EPSG:3067
   url  = "https://www.nic.funet.fi/index/geodata/ilmatiede/"

   match quantity:
      case "mean temperature":
         url += f"10km_daily_mean_temperature/geotiff/tday_{year}.tif"
      case "minimum temperature":
         url += f"10km_daily_minimum_temperature/geotiff/tmin_{year}.tif"
      case "maximum temperature":
         url += f"10km_daily_maximum_temperature/geotiff/tmax_{year}.tif"

   temp = rioxarray.open_rasterio(url).rio.reproject("EPSG:3067")

   return temp

#--------------------------------------------------------------------
def extract_month(temp_data,crs,year=2024,month=1,day=None):
   days = day_range(year,month,day)

   if(len(days)==2):
      temp_daily = temp_data[days[0]:days[1]]
      temp_daily = xr.where(temp_daily >= -1e30, temp_daily, np.nan)
      temp_daily = temp_daily.where(temp_daily == temp_daily)
      temp_daily = temp_daily.rio.write_crs(temp_data.rio.crs)

      temp_avg = temp_daily.mean(dim='band')
      temp_avg = temp_avg.rio.write_crs(temp_data.rio.crs)

      return temp_daily,temp_avg
   else:
      temp_daily = temp_data[days[0]]
      temp_daily = xr.where(temp_daily >= -1e30, temp_daily, np.nan)
      temp_daily = temp_daily.where(temp_daily == temp_daily)
      temp_daily = temp_daily.rio.write_crs(temp_data.rio.crs)

      return temp_daily,temp_daily

#--------------------------------------------------------------------
def temperature_by_municipality(municipalities,temp_raster,
         quantity="mean temperature"):

   muni_temp = np.array([])
   for gml_id in municipalities['GML_ID']:
      municipality = municipalities[municipalities['GML_ID'] == gml_id]
      try:
         clipped_raster = temp_raster.rio.clip(municipality.geometry,
                                 municipality.crs)
         match quantity:
            case "mean temperature":
               muni_temp = np.append(muni_temp,clipped_raster.mean().values)
            case "minimum temperature":
               muni_temp = np.append(muni_temp,clipped_raster.min().values)
            case "maximum temperature":
               muni_temp = np.append(muni_temp,clipped_raster.max().values)
      except:
         muni_temp = np.append(muni_temp,np.nan)

   return muni_temp

#--------------------------------------------------------------------
def extract_temperatures(municipalities,quantity="mean temperature",
         year=2024,month=1,day=None):
   temp_data = read_data(quantity,year)
   assert temp_data.rio.crs == municipalities.crs, \
               "CRS mismatch between the raster and the GeoDataFrame."

   temp_daily,temp_avg = extract_month(temp_data,temp_data.rio.crs,
            year,month,day)
   assert temp_data.rio.crs == temp_daily.rio.crs == temp_avg.rio.crs, \
               "CRS mismatch between the rasters."

   muni_temp = temperature_by_municipality(municipalities,temp_daily,
            quantity)

   return muni_temp

#--------------------------------------------------------------------
def plot_temperature_evolution(temperatures,city,
         year_start=1961,year_end=2024,month=1):

   # Read the temperature for the selected city
   year      = np.array(temperatures["year"]).astype(int)
   mean_temp = np.array(temperatures[f"mean_temp_{city}"]).astype(float)
   min_temp  = np.array(temperatures[f"min_temp_{city}"]).astype(float)
   max_temp  = np.array(temperatures[f"max_temp_{city}"]).astype(float)

   colnames = temperatures.columns
   y_min    = temperatures.loc[:,colnames.str.contains("min_")].min(axis=None)
   y_max    = temperatures.loc[:,colnames.str.contains("max_")].max(axis=None)

   ###########################
   ### Add exponential fit ###
   ###########################

   fig = plt.figure(figsize=(10,7))

   plt.scatter(year,mean_temp, marker='o', s=8,color='black')
   plt.scatter(year,min_temp,  marker='o', s=8,color='blue')
   plt.scatter(year,max_temp,  marker='o', s=8,color='red')
#   plt.fill_between(year,mean_temp,max_temp,
#         color='xkcd:deep sky blue',alpha=0.4,linewidth=0)
#   plt.fill_between(year,min_temp,mean_temp,
#         color='xkcd:bright red',alpha=0.4,linewidth=0)

   plt.title(f"Temperature in {month_name[month]} in {city}")
   plt.xlabel("Year")
   plt.ylabel("Temperature (Celsius degrees)")
   plt.xlim([year_start-1,year_end+1])
   plt.ylim([y_min-1,y_max+1])

   return fig

#--------------------------------------------------------------------



parser = argparse.ArgumentParser()
parser.add_argument("--month", type=int, default=1,
                    help="Month (number) of interest.")

args  = parser.parse_args()
month = args.month


# Paths definition
NOTEBOOK_PATH  = pathlib.Path().resolve()
DATA_DIRECTORY = NOTEBOOK_PATH / "data"
FIG_DIRECTORY  = NOTEBOOK_PATH / "figures"


# Read the Finland municipalities (as defined in 2021)
municipalities = geopandas.read_file(DATA_DIRECTORY /
         "finland_municipalities_2021.gpkg").to_crs("EPSG:3067")
         #.to_crs("EPSG:4326")


# Select the municipalities of interest
cities = ["Helsinki", "Jyväskylä", "Oulu", "Rovaniemi", "Inari"]

indexes = [i for i in municipalities.index 
               if municipalities["NAMEFIN"][i] in cities]
municipalities = municipalities.iloc[indexes]


# Initialisation temperature DataFrame for the selected month
temperatures = pd.DataFrame(columns=["year"]+[ f"{quantity}_temp_{city}"
            for city,quantity in product(cities,["mean","min","max"])])


for year in tqdm(range(1961,2025,1)):
   mean_temp = extract_temperatures(municipalities,"mean temperature",
                     year=year,month=month)
   min_temp  = extract_temperatures(municipalities,"minimum temperature",
                     year=year,month=month)
   max_temp  = extract_temperatures(municipalities,"maximum temperature",
                     year=year,month=month)


   d = {"year": year}
   for city in cities:
      idx = np.where(municipalities["NAMEFIN"]==city)
      d[f"mean_temp_{city}"] = [float(mean_temp[idx][0])]
      d[f"min_temp_{city}"]  = [float(min_temp[idx][0])]
      d[f"max_temp_{city}"]  = [float(max_temp[idx][0])]

   new_row      = pd.DataFrame(d)
   temperatures = pd.concat([temperatures, new_row], ignore_index=True)



for city in cities:
   fig = plot_temperature_evolution(temperatures,city,
                  year_start=1961,year_end=2024,month=month)
   fig.savefig(FIG_DIRECTORY / "Temperatures_evolution" /
               f"Temperature_{month_name[month]}_{city}.png")


