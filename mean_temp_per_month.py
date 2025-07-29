import pathlib
import argparse
import geopandas
import rioxarray
import matplotlib.pyplot as plt
from calendar import month_name
from datetime import datetime
from geoanalysis_functions import extract_month,temperature_by_municipality


#--------------------------------------------------------------------
def read_data(path,year=2024):
   """
   Read the data for the year of interest.
   
   Parameters
   ----------
   path: str
         Path of the data folder
   year: int
         year of interest

   Returns
   -------
   municipalities: geopandas.geodataframe.GeoDataFrame
         Municipalities with geospatial information reprojected to EPSG:3067
   temp: xarray.core.dataarray.DataArray
         3d raster data reprojected to EPSG:3067
   """
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
def plot_temperature_maps(municipalities,temp_avg,vmin=None,vmax=None,
         year=2024,month=1,day=None):
   """
   Read the data for the year of interest.
   
   Parameters
   ----------
   municipalities: geopandas.geodataframe.GeoDataFrame
         Municipalities of interest
   temp_avg: xarray.core.dataarray.DataArray
         2d raster data over the period of interest
   year: int
         year of interest
   month: int
         month of interest
   day: int
         day of interest (if specified)

   Returns
   -------
   matplotlib.figure.Figure
         Map of the 2d raster data with the contours of the municipalities,
         and map of the mean temperature within each municipality.
   """
   # Reproject to EPSG:4326
   municipalities = municipalities.to_crs("EPSG:4326")
   temp_avg = temp_avg.rio.reproject("EPSG:4326")
   assert temp_avg.rio.crs == municipalities.crs, \
            "CRS mismatch which can cause problems when plotting."
   minx, miny, maxx, maxy = municipalities.total_bounds

   fig, axs = plt.subplots(1, 2, figsize=(10,7))

   im = temp_avg.plot(ax=axs[0], cmap='terrain', vmin=vmin,vmax=vmax,
            add_colorbar=True)
   municipalities.boundary.plot(ax=axs[0], color='k', alpha=0.25)

   colorbar = im.colorbar
   colorbar.set_label('Mean temperature (Celsius degrees)')
   axs[0].set_title('10x10 km grid')
   axs[0].set_xlabel("Eastern longitude (degrees)")
   axs[0].set_ylabel("Latitude (degrees)")
   axs[0].set_xlim(minx, maxx)
   axs[0].set_ylim(miny, maxy)



   municipalities.plot(ax=axs[1], column="temperature", cmap="terrain",
          vmin=vmin,vmax=vmax, linewidth=0,  legend=True,
          legend_kwds={"label": "Mean temperature (Celsius degrees)"}
          )
   axs[1].set_title('By municipality')
   axs[1].set_xlabel("Eastern longitude (degrees)")
#   axs[1].set_ylabel("Latitude (degrees)")
   axs[1].set_xlim(minx, maxx)
   axs[1].set_ylim(miny, maxy)


   if(day is None):
      fig.suptitle(f"{month_name[month]} {year}")
   else:
      fig.suptitle(f"{day:02d} {month_name[month]} {year}")

   return fig

#--------------------------------------------------------------------



parser = argparse.ArgumentParser()
parser.add_argument("--year",  type=int, default=2024,
                    help="Year of interest.")
parser.add_argument("--month", type=int, default=1,
                    help="Month (number) of interest.")
parser.add_argument("--day", type=int,
                    help="Day of interest.")
parser.add_argument("--vmin", type=float,
                    help="Minimum of colorscale.")
parser.add_argument("--vmax", type=float,
                    help="Maximum of colorscale.")

args  = parser.parse_args()
year  = args.year
month = args.month
day   = args.day
vmin  = args.vmin
vmax  = args.vmax


# Paths definition
NOTEBOOK_PATH  = pathlib.Path().resolve()
DATA_DIRECTORY = NOTEBOOK_PATH / "data"
FIG_DIRECTORY  = NOTEBOOK_PATH / "figures"


municipalities,temp_data = read_data(DATA_DIRECTORY,year)
assert temp_data.rio.crs == municipalities.crs, \
            "CRS mismatch between the raster and the GeoDataFrame."


temp_daily,temp_avg = extract_month(temp_data,year=year,month=month,day=day)
assert temp_data.rio.crs == temp_daily.rio.crs == temp_avg.rio.crs, \
            "CRS mismatch between the rasters."


muni_temp = temperature_by_municipality(municipalities,temp_daily)
municipalities["temperature"] = muni_temp


fig = plot_temperature_maps(municipalities,temp_avg,vmin=vmin,vmax=vmax,
         year=year,month=month,day=day)
if(day is None):
   fig.savefig(FIG_DIRECTORY / 
               f"Mean_temperature_{month_name[month]}_{year}.png")
else:
   fig.savefig(FIG_DIRECTORY / 
               f"Mean_temperature_{day:02d}_{month_name[month]}_{year}.png")


