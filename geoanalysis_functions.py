import xarray as xr
import numpy as np
import rioxarray
from calendar import monthrange
from datetime import datetime


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
def temperature_by_municipality(municipalities,temp_raster):
   muni_temp = np.array([])
   for gml_id in municipalities['GML_ID']:
      municipality = municipalities[municipalities['GML_ID'] == gml_id]
      try:
         clipped_raster = temp_raster.rio.clip(municipality.geometry,
                                 municipality.crs)
         muni_temp = np.append(muni_temp,clipped_raster.mean().values)
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

   muni_temp = temperature_by_municipality(
            municipalities,temp_daily)

   return muni_temp,temp_avg

#--------------------------------------------------------------------


