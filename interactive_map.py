import pathlib
import folium
from folium.plugins import GroupedLayerControl
import base64

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


def get_color(x):
    decimals = 2
    x = np.around(x, decimals=decimals)
    ls = np.linspace(0,1,10**decimals+1)
    if x==-9999:
        return (0, 0, 0, 0)
    elif x != -9999:
        return cm.get_cmap('terrain')(ls)[np.argwhere(ls==x)]
    else:
        raise ValueError()
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
    interactive=False
)


# Read the Finland municipalities (as defined in 2021)
municipalities = geopandas.read_file(DATA_DIRECTORY /
         "finland_municipalities_2021.gpkg").to_crs("EPSG:3067")
         #.to_crs("EPSG:4326")


# Select the municipalities of interest
cities = ["Helsinki", "Jyväskylä", "Oulu", "Rovaniemi", "Inari"]

indexes = [i for i in municipalities.index 
               if municipalities["NAMEFIN"][i] in cities]
municipalities = municipalities.iloc[indexes]

images = np.array([])
for city in municipalities['NAMEFIN']:
   images = np.append(images,city)

temp_layers = {}
muni_layers = {}
k=0


svg = """
<object data="data:image/jpg;base64,{}" width="{}" height="{} type="image/svg+xml">
</object>""".format

width, height, fat_wh = 200, 140, 1.3


for year in tqdm([1961,1981,2001,2021]):
   for month in [1,7]:
      mean_temp,map = extract_temperatures(municipalities,"mean temperature",
                        year=year,month=month)

      map = map.rio.reproject("EPSG:4326")
      minmax = float(map.min().values),float(map.max().values)
      map.attrs['_FillValue']=-9999
      map = map.fillna(value=-9999)

      if((year == 1961) & (month == 1)):
         visible = True
      else:
         visible = False

      boundList = [x for x in map.rio.bounds()]
      temp_layers.update({f"{k}": 
          folium.raster_layers.ImageOverlay(
             image=map.values,
             bounds=[[boundList[1],boundList[0]], [boundList[3],boundList[2]]],
             opacity=0.6,
             name=f"{month_name[month]} {year}",
             colormap=cm.terrain,
             vmin=minmax[0], vmax=minmax[1],#lambda x: get_color(x)
             show=visible
          )
      })
      temp_layers[f"{k}"].add_to(interactive_map)


      municipalities["temperature"] = (100*mean_temp).astype('int')/100
#      municipalities["image_url"] = "figures/test.jpg"



#      popup  = folium.Popup(iframe, parse_html = True, max_width=1500)

      popup_html = np.array([])
      for idx in range(len(municipalities)):
         encoded = base64.b64encode(open("figures/test.png", 'rb').read()).decode()
#         iframe = IFrame(svg(encoded.decode('UTF-8'), width, height) , width=width*fat_wh, height=height*fat_wh)
         string = f'<img src="data:image/png;base64,{encoded}" width="300">'

#         municipality = municipalities.iloc[idx]
#         string = f"""<img src="{municipality['image_url']}"
#                  width="200" height="142">
#               """
         popup_html = np.append(popup_html,string)
      municipalities['popup_html'] = popup_html

      popup = folium.GeoJsonPopup(
         fields=["popup_html"],
         aliases=[""],
         labels=True,
         sticky=False,
         localize=True,
      )

      # Define custom tooltip with HTML
      tooltip = folium.features.GeoJsonTooltip(
         fields=("NAMEFIN","temperature",),
         aliases=("Municipality:","Mean temperature:",),
#         labels=True,
#         sticky=True,
         localize=True,
         )
#      tooltip = folium.features.GeoJsonTooltip(
#         fields=("tooltip_html",),
##         aliases=(""),
##         labels=True,
##         sticky=True,
#         localize=True,
#         )

      muni_layers.update({f"{k}":
         folium.GeoJson(
            municipalities,
            style_function=lambda feature: {
                  "fillColor": "transparent",
                  "color": "black",
                  "weight": 2
                  },
            opacity=0.8,
            name=f"{month_name[month]} {year}",
            show=False,
            tooltip=tooltip,
            popup=popup
         )
      })

      
      muni_layers[f"{k}"].add_to(interactive_map)
      k = k+1

#interactive_map

## Read the Finland municipalities (as defined in 2021)
#municipalities = geopandas.read_file(DATA_DIRECTORY /
#         "finland_municipalities_2021.gpkg").to_crs("EPSG:3067")
#         #.to_crs("EPSG:4326")


## Select the municipalities of interest
#cities = ["Helsinki", "Jyväskylä", "Oulu", "Rovaniemi", "Inari"]

#indexes = [i for i in municipalities.index 
#               if municipalities["NAMEFIN"][i] in cities]
#municipalities = municipalities.iloc[indexes]


## Initialisation temperature DataFrame for the selected month
#temperatures = pd.DataFrame(columns=["year"]+[ f"{quantity}_temp_{city}"
#            for city,quantity in product(cities,["mean","min","max"])])


#for year in tqdm(range(1961,2025,1)):
#   mean_temp = extract_temperatures(municipalities,"mean temperature",
#                     year=year,month=month)
#   min_temp  = extract_temperatures(municipalities,"minimum temperature",
#                     year=year,month=month)
#   max_temp  = extract_temperatures(municipalities,"maximum temperature",
#                     year=year,month=month)


#   d = {"year": year}
#   for city in cities:
#      idx = np.where(municipalities["NAMEFIN"]==city)
#      d[f"mean_temp_{city}"] = [float(mean_temp[idx][0])]
#      d[f"min_temp_{city}"]  = [float(min_temp[idx][0])]
#      d[f"max_temp_{city}"]  = [float(max_temp[idx][0])]

#   new_row      = pd.DataFrame(d)
#   temperatures = pd.concat([temperatures, new_row], ignore_index=True)



#for city in cities:
#   fig = plot_temperature_evolution(temperatures,city,
#                  year_start=1961,year_end=2024,month=month)
#   fig.savefig(FIG_DIRECTORY / "Temperatures_evolution" /
#               f"Temperature_{month_name[month]}_{city}.png")

folium.LayerControl(collapsed=False).add_to(interactive_map)

list_layers = [temp_layers[key] for key in temp_layers.keys()]
GroupedLayerControl(
   groups={'Mean temperature': list_layers},
   collapsed=False,
   exclusive_groups=True
).add_to(interactive_map)
list_layers = [muni_layers[key] for key in temp_layers.keys()]
GroupedLayerControl(
   groups={'Municipalities': list_layers},
   collapsed=False,
   exclusive_groups=False
).add_to(interactive_map)

interactive_map.save(HTML_DIRECTORY /
            "map_temperature_municipalities.html")


