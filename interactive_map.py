import pathlib
import geopandas
import numpy as np

import folium
import base64
import branca.colormap
from geocube.vector import vectorize
from folium.plugins import GroupedLayerControl

from calendar import month_name
from geoanalysis_functions import extract_temperatures


from branca.element import MacroElement
from jinja2 import Template


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
def new_image_layer(raster,name="",visible=False,
         year=2024,month=1,day=None):
   """
   Produce a new layer based on raster data
   
   Parameters
   ----------
   raster: xarray.core.dataarray.DataArray
         2d raster data to include
   name: str
         Name of the layer
   visible: boolean
         Define if the layer will be displayed when opening the map
   year: int
         year to be added in name
   month: int
         month to be added in name
   day: int
         day to be added in name

   Returns
   -------
   layer: folium.features.GeoJson
         layer to be added to the interactive map
   cmap: branca.colormap.LinearColormap
         corresponding color map
   """

   raster = raster.rio.reproject("EPSG:4326")
   raster = raster.rename("raster")
   gdf    = vectorize(raster.astype("float32"))
   gdf["id"] = gdf.index.astype(str)


   if(month==1):
      colors=branca.colormap.linear.Blues_05.colors
      invert=[colors[-i-1] for i in range(len(colors))]
      cmap = branca.colormap.LinearColormap(
            colors=invert,
            vmin=-15,vmax=0,
            caption="Mean temperature")
      fill_color="Blues_r"
   else:
      cmap = branca.colormap.LinearColormap(
            colors=branca.colormap.linear.YlOrRd_08.colors,
            vmin=10,vmax=20,
            caption="Mean temperature")
      fill_color="YlOrRd"

#   layer = folium.Choropleth(
#         geo_data=gdf,
#         data=gdf,
#         columns=("id","_data"),
#         key_on="feature.id",

#         bins=10,
#         fill_color=fill_color,
##         colormap=cmap,
#         line_weight=0,
#         opacity=0.6,
#         legend_name="Mean temperature",
#         name=name,
#         show=visible,

##         highlight=True
#   )

   temp_dict = gdf.set_index('id')['raster']

   layer = folium.GeoJson(
      gdf,
      style_function=lambda feature: {
            "fillColor": cmap(temp_dict[feature['id']]),
            "color": "black",
            "weight": 0,
            "fillOpacity": 0.6,
            },
      name=name,
      show=visible
   )

   for child in layer._children:
        if child.startswith("color_map"):
            del layer._children[child]

   return layer,cmap

#--------------------------------------------------------------------
def new_polygon_layer(gdf,name="",image=None,
         year=2024,month=1,day=None):
   """
   Produce a polygon layer based on a GeoDataFrame
   
   Parameters
   ----------
   gdf: geopandas.geodataframe.GeoDataFrame
         Municipalities of interest with temperatures information
   name: str
         Name of the layer
   image: str
         Path of the image to be included in a popup
   year: int
         year to be added in name
   month: int
         month to be added in name
   day: int
         day to be added in name

   Returns
   -------
   folium.features.GeoJson
         layer to be added to the interactive map
   """

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


copyright  = 'Temperature data (c) <a href="https://en.ilmatieteenlaitos.fi/">'
copyright += 'Finnish Meteorological Institute</a> & '
copyright += '<a href="https://paituli.csc.fi/download.html">Paituli</a>, '
copyright += 'Map data (c) <a href="http://www.openstreetmap.org/copyright">'
copyright += 'OpenStreetMap</a> contributors.'

# Initial map
interactive_map = folium.Map(
    max_bounds=True,
    location=[64.5, 27.1],
    zoom_start=5,
    min_lat=59.7,
    max_lat=70.1,
    min_lon=21.2,
    max_lon=31.6,
    interactive=True,
    attr = copyright
)



# Read the Finland municipalities (as defined in 2021)
municipalities = geopandas.read_file(DATA_DIRECTORY /
         "finland_municipalities_weather_evolution.gpkg"
         ).to_crs("EPSG:3067")


temp_layers = {}
cmap_layers = {}
muni_layers = {}
k1,k2   = 0,0
visible = True

for month in [1,7]:
   for year in [1961,1981,2001,2021]:
      mean_temp,map = extract_temperatures(municipalities,"mean temperature",
                        year=year,month=month)


      folium_layer,cmap = new_image_layer(map,f"{month_name[month]} {year}",
               visible,year,month)
      temp_layers.update({f"{k1}": folium_layer})
      cmap_layers.update({f"{k1}": cmap})

#      interactive_map.add_child(cmap_layers[f"{k1}"])
      interactive_map.add_child(temp_layers[f"{k1}"])
#      temp_layers[f"{k1}"].add_to(interactive_map)
#      if(k1>=4):
#         bc = BindColormap(temp_layers[f"{k1}"], cmap_layers["4"])
#      else:
#         bc = BindColormap(temp_layers[f"{k1}"], cmap_layers["0"])
#      interactive_map.add_child(bc)

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

interactive_map.add_child(cmap_layers["0"])
interactive_map.add_child(cmap_layers["4"])


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


