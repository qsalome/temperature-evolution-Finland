import xarray as xr
import numpy as np
import rioxarray
from calendar import monthrange
from datetime import datetime


#--------------------------------------------------------------------
def read_data(quantity="mean temperature",year=2024):
   """
   Read, reproject and return the requested raster data.
   
   Parameters
   ----------
   quantity: str
         data set to be read
   year: int
         year of interest for the data set

   Returns
   -------
   xarray.core.dataarray.DataArray
        3d raster data reprojected to EPSG:3067
   """
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
   
   Parameters
   ----------
   year: int
         year of interest
   month: int
         month of interest
   day: int
         if specify, returns only one value

   Returns
   -------
   tuple
         Yearly day numbers for the first and last days of the month,
         or for the specified day.
   """
   if(day is None):
      last_day  = monthrange(year, month)[1]

      first_day = datetime(year, month, 1).timetuple().tm_yday
      last_day  = datetime(year, month, last_day).timetuple().tm_yday

      return first_day,last_day
   else:
      return (datetime(year, month, day).timetuple().tm_yday,)

#--------------------------------------------------------------------
def extract_month(data,quantity="mean temperature",
         year=2024,month=1,day=None):
   """
   Extract a subset of the data corresponding to the period of interest.
   
   Parameters
   ----------
   data: xarray.core.dataarray.DataArray
         Input raster data for the full year
   quantity: str
         data set to consider
   year: int
         year of interest
   month: int
         month of interest
   day: int
         day of interest (if specified)

   Returns
   -------
   data_daily: xarray.core.dataarray.DataArray
         3d daily raster data for the period of interest (month or day).
   monthly: xarray.core.dataarray.DataArray
         2d raster data of the average or sum over the period of interest.
   """
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
   """
   Extract the area of the raster data corresponding to each municipality
   and calculate the mean value.
   
   Parameters
   ----------
   municipalities: geopandas.geodataframe.GeoDataFrame
         Municipalities of interest
   temp_raster: xarray.core.dataarray.DataArray
         2d or 3d raster data to analyse

   Returns
   -------
   numpy.ndarray
         Average value within each municipalities.
   """
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
   """
   Read the raster data, extract the period of interest and calculate
   the mean value within each municipality.
   
   Parameters
   ----------
   municipalities: geopandas.geodataframe.GeoDataFrame
         Municipalities of interest
   quantity: str
         data set to consider
   year: int
         year of interest
   month: int
         month of interest
   day: int
         day of interest (if specified)

   Returns
   -------
   muni_temp: numpy.ndarray
         Average monthly temperature within each municipalities.
   temp_avg: xarray.core.dataarray.DataArray
         2d raster data over the period of interest
   """
   temp_data = read_data(quantity,year)
   assert temp_data.rio.crs == municipalities.crs, \
               "CRS mismatch between the raster and the GeoDataFrame."

   temp_daily,temp_avg = extract_month(temp_data,quantity,year,month,day)
   assert temp_data.rio.crs == temp_daily.rio.crs == temp_avg.rio.crs, \
               "CRS mismatch between the rasters."

   muni_temp = temperature_by_municipality(municipalities,temp_daily)

   return muni_temp,temp_avg

#--------------------------------------------------------------------
def precipitations_by_municipality(municipalities,rain_raster):
   """
   Extract the area of the raster data corresponding to each municipality
   and calculate the mean value.
   
   Parameters
   ----------
   municipalities: geopandas.geodataframe.GeoDataFrame
         Municipalities of interest
   rain_raster: xarray.core.dataarray.DataArray
         2d or 3d raster data to analyse

   Returns
   -------
   numpy.ndarray
         Average value within each municipalities.
   """
   muni_rain = np.array([])
   for gml_id in municipalities['GML_ID']:
      municipality = municipalities[municipalities['GML_ID'] == gml_id]
      try:
         clipped_raster = rain_raster.rio.clip(municipality.geometry,
                                 municipality.crs)
         rain_avg  = clipped_raster.mean().values
         muni_rain = np.append(muni_rain,rain_avg)
      except:
         muni_rain = np.append(muni_rain,np.nan)

   return muni_rain

#--------------------------------------------------------------------
def extract_precipitations(municipalities,quantity="precipitations",
         year=2024,month=1,day=None):
   """
   Read the raster data, extract the period of interest and calculate
   the mean value within each municipality.
   
   Parameters
   ----------
   municipalities: geopandas.geodataframe.GeoDataFrame
         Municipalities of interest
   quantity: str
         data set to consider
   year: int
         year of interest
   month: int
         month of interest
   day: int
         day of interest (if specified)

   Returns
   -------
   muni_rain: numpy.ndarray
         Average monthly precipitation within each municipalities.
   temp_avg: xarray.core.dataarray.DataArray
         2d raster data over the period of interest
   """
   rain_data = read_data(quantity,year)
   assert rain_data.rio.crs == municipalities.crs, \
               "CRS mismatch between the raster and the GeoDataFrame."

   rain_daily,rain_monthly = extract_month(rain_data,quantity,year,month,day)
   assert rain_data.rio.crs == rain_daily.rio.crs == rain_monthly.rio.crs, \
               "CRS mismatch between the rasters."

   muni_rain = precipitations_by_municipality(municipalities,rain_monthly)

   return muni_rain,rain_monthly

#--------------------------------------------------------------------


