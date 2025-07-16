import pathlib
import geopandas
import rioxarray
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from calendar import monthrange,month_name
from datetime import datetime
import argparse


#--------------------------------------------------------------------
def day_range(year=2024,month=1):
   """
   Determine the yearly day number for the first and last days
   of the selected month. Take into account leap years.
   """
   last_day  = monthrange(year, month)[1]

   first_day = datetime(year, month, 1).timetuple().tm_yday
   last_day  = datetime(year, month, last_day).timetuple().tm_yday

   return first_day,last_day

#--------------------------------------------------------------------
def read_data(path,year):
   # Read the Finland municipalities (as defined in 2021)
   municipalities = geopandas.read_file(path /
            "finland_municipalities_2021.gpkg").to_crs("EPSG:3067")
            #.to_crs("EPSG:4326")

   # Load the raster and reproject it to EPSG:3067
   url = (
      "https://www.nic.funet.fi/index/geodata/ilmatiede/"
      "10km_daily_mean_temperature/geotiff/"
      f"tday_{year}.tif"
   )

   temp = rioxarray.open_rasterio(url).rio.reproject("EPSG:3067")

#   temp = rioxarray.open_rasterio(path /
#            "10km_daily_mean_temperature/geotiff" /
#            f"tday_{year}.tif").rio.reproject("EPSG:3067")

   return municipalities,temp

#--------------------------------------------------------------------
def extract_month(temp_data,year,month,crs):
   days = day_range(year,month)

   temp_daily = temp_data[days[0]:days[1]]
   temp_daily = xr.where(temp_daily >= -1e30, temp_daily, np.nan)
   temp_daily = temp_daily.where(temp_daily == temp_daily)
   temp_daily = temp_daily.rio.write_crs(temp_data.rio.crs)

   temp_avg = temp_daily.mean(dim='band')
   temp_avg = temp_avg.rio.write_crs(temp_data.rio.crs)

   return temp_daily,temp_avg

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
def plot_temperature_maps(municipalities,temp_avg,year,month):
   # Reproject to EPSG:4326
   municipalities = municipalities.to_crs("EPSG:4326")
   temp_avg = temp_avg.rio.reproject("EPSG:4326")
   assert temp_avg.rio.crs == municipalities.crs, \
            "CRS mismatch which can cause problems when plotting."
   minx, miny, maxx, maxy = municipalities.total_bounds

   fig, axs = plt.subplots(1, 2, figsize=(10,7))

   im = temp_avg.plot(ax=axs[0], cmap='terrain', add_colorbar=True)
   municipalities.boundary.plot(ax=axs[0], color='k', alpha=0.25)

   colorbar = im.colorbar
   colorbar.set_label('Mean temperature (Celsius degrees)')
   axs[0].set_title('10x10 km grid')
   axs[0].set_xlabel("Eastern longitude (degrees)")
   axs[0].set_ylabel("Latitude (degrees)")
   axs[0].set_xlim(minx, maxx)
   axs[0].set_ylim(miny, maxy)



   municipalities.plot(ax=axs[1], column="temperature", cmap="terrain",
          linewidth=0,  legend=True,
          legend_kwds={"label": "Mean temperature (Celsius degrees)"}
          )
   axs[1].set_title('By municipality')
   axs[1].set_xlabel("Eastern longitude (degrees)")
#   axs[1].set_ylabel("Latitude (degrees)")
   axs[1].set_xlim(minx, maxx)
   axs[1].set_ylim(miny, maxy)

   fig.suptitle(f"{month_name[month]} {year}")

   return fig

#--------------------------------------------------------------------



parser = argparse.ArgumentParser()
parser.add_argument("--year",  type=int,
                    help="Year of interest.")
parser.add_argument("--month", type=int,
                    help="Month (number) of interest.")

args  = parser.parse_args()
year  = args.year
month = args.month


# Paths definition
NOTEBOOK_PATH  = pathlib.Path().resolve()
DATA_DIRECTORY = NOTEBOOK_PATH / "data"
FIG_DIRECTORY  = NOTEBOOK_PATH / "figures"


municipalities,temp_data = read_data(DATA_DIRECTORY,year)
assert temp_data.rio.crs == municipalities.crs, \
            "CRS mismatch between the raster and the GeoDataFrame."


temp_daily,temp_avg = extract_month(temp_data,year,month,temp_data.rio.crs)
assert temp_data.rio.crs == temp_daily.rio.crs == temp_avg.rio.crs, \
            "CRS mismatch between the rasters."


muni_temp = temperature_by_municipality(municipalities,
         temp_daily)
municipalities["temperature"] = muni_temp


fig = plot_temperature_maps(municipalities,temp_avg,year,month)
fig.savefig(FIG_DIRECTORY / 
            f"Mean_temperature_{month_name[month]}_{year}.png")


