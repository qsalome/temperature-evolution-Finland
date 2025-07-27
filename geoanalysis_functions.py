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
      case "precipitations":
         url += f"10km_daily_precipitation/geotiff/rrday_{year}.tif"

   data = rioxarray.open_rasterio(url).rio.reproject("EPSG:3067")

   return data

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
def extract_month(data,crs,quantity="mean temperature",
         year=2024,month=1,day=None):
   days = day_range(year,month,day)

   if(len(days)==2):
      data_daily = data[days[0]:days[1]+1]
      data_daily = xr.where(data_daily >= -1e30, data_daily, np.nan)
      data_daily = data_daily.where(data_daily == data_daily)
      data_daily = data_daily.rio.write_crs(data.rio.crs)

      if(quantity=="precipitations"):
         monthly = data_daily.sum(dim='band')
         monthly = monthly.rio.write_crs(data.rio.crs)
      else:
         monthly = data_daily.mean(dim='band')
         monthly = monthly.rio.write_crs(data.rio.crs)

      return data_daily,monthly
   else:
      data_daily = data[days[0]]
      data_daily = xr.where(data_daily >= -1e30, data_daily, np.nan)
      data_daily = data_daily.where(data_daily == data_daily)
      data_daily = data_daily.rio.write_crs(data.rio.crs)

      return data_daily,data_daily

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
            quantity,year,month,day)
   assert temp_data.rio.crs == temp_daily.rio.crs == temp_avg.rio.crs, \
               "CRS mismatch between the rasters."

   muni_temp = temperature_by_municipality(municipalities,temp_daily)

   return muni_temp,temp_avg

#--------------------------------------------------------------------
def precipitations_by_municipality(municipalities,rain_raster):
   muni_rain = np.array([])
   for gml_id in municipalities['GML_ID']:
      municipality = municipalities[municipalities['GML_ID'] == gml_id]
      try:
         clipped_raster = rain_raster.rio.clip(municipality.geometry,
                                 municipality.crs)
         rain_avg  = clipped_raster.mean().values#/100
         muni_rain = np.append(muni_rain,rain_avg)
      except:
         muni_rain = np.append(muni_rain,np.nan)

   return muni_rain

#--------------------------------------------------------------------
def extract_precipitations(municipalities,quantity="precipitations",
         year=2024,month=1,day=None):
   rain_data = read_data(quantity,year)
   assert rain_data.rio.crs == municipalities.crs, \
               "CRS mismatch between the raster and the GeoDataFrame."

   rain_daily,rain_monthly = extract_month(rain_data,rain_data.rio.crs,
            quantity,year,month,day)
   assert rain_data.rio.crs == rain_daily.rio.crs == rain_monthly.rio.crs, \
               "CRS mismatch between the rasters."

   muni_rain = precipitations_by_municipality(municipalities,rain_monthly)

   return muni_rain,rain_monthly

#--------------------------------------------------------------------


