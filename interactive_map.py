import pathlib
import folium
from folium.plugins import GroupedLayerControl
import base64

from geocube.vector import vectorize

import branca
from branca.element import MacroElement
from jinja2 import Template

import argparse
import pandas as pd
import geopandas
import rioxarray
import xarray as xr
import numpy as np
from calendar import monthrange,month_name
from datetime import datetime
from tqdm import tqdm
from itertools import product
import lmfit
import fit_functions
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.cm as cm
from folium import IFrame


#--------------------------------------------------------------------
class BindColormap(MacroElement):
    """Binds a colormap to a given layer.

    Parameters
    ----------
    colormap : branca.colormap.ColorMap
        The colormap to bind.
    """
    def __init__(self, layer, colormap):
        super(BindColormap, self).__init__()
        self.layer = layer
        self.colormap = colormap
        self._template = Template(u"""
        {% macro script(this, kwargs) %}
            {{this.colormap.get_name()}}.svg[0][0].style.display = 'block';
            {{this._parent.get_name()}}.on('overlayadd', function (eventLayer) {
                if (eventLayer.layer == {{this.layer.get_name()}}) {
                    {{this.colormap.get_name()}}.svg[0][0].style.display = 'block';
                }});
            {{this._parent.get_name()}}.on('overlayremove', function (eventLayer) {
                if (eventLayer.layer == {{this.layer.get_name()}}) {
                    {{this.colormap.get_name()}}.svg[0][0].style.display = 'none';
                }});
        {% endmacro %}
        """)

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

   muni_temp = temperature_by_municipality(municipalities,temp_daily)

   return muni_temp,temp_avg

#--------------------------------------------------------------------
def new_image_layer(raster,interactive_map,name="",visible=False,
         year=2024,month=1,day=None):
   raster = raster.rio.reproject("EPSG:4326")
   gdf    = vectorize(raster.astype("float32"))
   gdf["id"] = gdf.index.astype(str)
   boundList = [x for x in raster.rio.bounds()]
   bounds=[[boundList[1],boundList[0]], [boundList[3],boundList[2]]]

   minmax = float(raster.min().values),float(raster.max().values)
   print(minmax)
   raster.attrs['_FillValue']=-9999
   raster = raster.fillna(value=-9999)

#   layer = folium.raster_layers.ImageOverlay(
#         image=raster.values,
#         bounds=[[boundList[1],boundList[0]], [boundList[3],boundList[2]]],
#         opacity=0.6,
#         name=name,
#         colormap=cm.terrain,
#         vmin=minmax[0], vmax=minmax[1],#lambda x: get_color(x)
#         show=visible
#   )

   layer = folium.Choropleth(
         geo_data=gdf,
         data=gdf,
         columns=("id","_data"),
         key_on="feature.id",

         bins=9,
         fill_color="YlOrRd",
         line_weight=0,
         opacity=0.6,
         legend_name="Mean temperature",
         name=name,
         show=visible,

#         highlight=True
   )

   return layer
#--------------------------------------------------------------------
def new_polygon_layer(gdf,name="",image=None,
         year=2024,month=1,day=None):

   if image is not None:
      popup_html = np.array([])
      for idx in range(len(gdf)):
         city = gdf["NAMEFIN"][idx]
         imname  = image.name.format(f"{month_name[month]}",city)
         encoded = base64.b64encode(
                     open(image.with_name(imname), 'rb').read()).decode()
         string = f'<img src="data:image/png;base64,{encoded}" width="300">'
         popup_html = np.append(popup_html,string)
      gdf['popup_html'] = popup_html

      popup = folium.GeoJsonPopup(
         fields=["popup_html"],
         aliases=[""],
         labels=True,
         sticky=False,
         localize=True,
      )
   else:
      popup=None

   # Define custom tooltip with HTML
   tooltip = folium.features.GeoJsonTooltip(
      fields=("NAMEFIN",f"TEMP_EVO_{month_name[month]}",),
      aliases=("Municipality:","1961-2024 (°C):",),
      labels=True,
      sticky=False,
      localize=True,
      )

   layer = folium.GeoJson(
      gdf,
      style_function=lambda feature: {
            "fillColor": "transparent",
            "color": "black",
            "weight": 2
            },
      opacity=0.8,
      name=name,
      show=False,
      tooltip=tooltip,
      popup=popup
   )

   return layer

#--------------------------------------------------------------------



# Paths definition
NOTEBOOK_PATH  = pathlib.Path().resolve()
DATA_DIRECTORY = NOTEBOOK_PATH / "data"
FIG_DIRECTORY  = NOTEBOOK_PATH / "figures"
HTML_DIRECTORY = NOTEBOOK_PATH / "html"


# Initial map
interactive_map = folium.Map(
    max_bounds=True,
    location=[64.5, 27.1],
    zoom_start=5,
    min_lat=59.7,
    max_lat=70.1,
    min_lon=21.2,
    max_lon=31.6,
    interactive=True
)


YlGn_cmap = branca.colormap.LinearColormap(
    colors=branca.colormap.linear.YlGn_09.colors,
    vmin=0,
    vmax=state_data.Unemployment.max(),  # setting max value for legend
    caption='Unemployment Rate (%)'
)


YlGn_cmap = branca.colormap.LinearColormap(
    colors=branca.colormap.linear.YlGn_09.colors,
    vmin=0,
    vmax=state_data.Unemployment.max(),  # setting max value for legend
    caption='Unemployment Rate (%)'
)


# Read the Finland municipalities (as defined in 2021)
municipalities = geopandas.read_file(DATA_DIRECTORY /
         "finland_municipalities_temperature_evolution.gpkg"
         ).to_crs("EPSG:3067")


temp_layers = {}
muni_layers = {}
k1,k2   = 0,0
visible = True

for month in [1,7]:
   for year in tqdm([1961,1981,2001,2021]):
      mean_temp,map = extract_temperatures(municipalities,"mean temperature",
                        year=year,month=month)


      folium_layer = new_image_layer(map,interactive_map.crs,
               f"{month_name[month]} {year}",
               visible,year,month)
      temp_layers.update({f"{k1}": folium_layer})
      temp_layers[f"{k1}"].add_to(interactive_map)

      visible = False
      k1 = k1+1


   imname = FIG_DIRECTORY / "Temperatures_evolution"
   imname = imname / "Temperature_{}_{}.png"
   folium_layer = new_polygon_layer(municipalities,
               f"{month_name[month]}",image=imname,
               year=year,month=month)
   muni_layers.update({f"{k2}": folium_layer})
   muni_layers[f"{k2}"].add_to(interactive_map)

   k2 = k2+1



folium.LayerControl(collapsed=False).add_to(interactive_map)

list_layers = [temp_layers[key] for key in temp_layers.keys()]
GroupedLayerControl(
   groups={'Mean temperature': list_layers},
   collapsed=False,
   exclusive_groups=True
).add_to(interactive_map)
list_layers = [muni_layers[key] for key in muni_layers.keys()]
GroupedLayerControl(
   groups={'Municipalities': list_layers},
   collapsed=False,
   exclusive_groups=True
).add_to(interactive_map)

interactive_map.save(HTML_DIRECTORY /
            "map_temperature_municipalities.html")


