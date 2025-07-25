import pathlib
import argparse
import pandas as pd
import geopandas
import numpy as np
from calendar import month_name
from tqdm import tqdm
from itertools import product
import lmfit
import fit_functions
from geoanalysis_functions import extract_temperatures
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


#--------------------------------------------------------------------
def fit_linear(y,x=None,err=None,method='leastsq'):
   """!
   Minimize the data points with the selected function
   
   @param y: 1D array : data points to fit
   @param params: lmfit method : list of Gaussian
   @param err: Float : rms of spectrum to fit
   @param method: String : minimization method
   
   @return lmfit obj
   """
   if x is None:
      x = np.arange(len(y))

   fit_params = lmfit.create_params(
            coeff={'value': 2, 'min': -10, 'max': 10},
            const={'value': np.mean(y), 'min': -60, 'max': 40})

   fitter = lmfit.Minimizer(fit_functions.linear,
               fit_params,fcn_args=(x,y,err))
   try:
      fit = fitter.minimize(method=method)
   except Exception as mes:
      print("Something wrong with fit: ", mes)
      raise SystemExit

   model = fit_functions.linear(fit.params,x)

   return model

#--------------------------------------------------------------------
def fit_exponential(y,x=None,err=None,method='leastsq'):
   """!
   Minimize the data points with the selected function
   
   @param y: 1D array : data points to fit
   @param params: lmfit method : list of Gaussian
   @param err: Float : rms of spectrum to fit
   @param method: String : minimization method
   
   @return lmfit obj
   """
   if x is None:
      x = np.arange(len(y))

   fit_params = lmfit.create_params(
            amp={'value': 1, 'min': -20, 'max': 20},
            coeff={'value': 1e-2, 'min': -1, 'max': 1},
            const={'value': np.mean(y), 'min': -60, 'max': 40})

   fitter = lmfit.Minimizer(fit_functions.exponential,
               fit_params,fcn_args=(x,y,err))
   try:
      fit = fitter.minimize(method=method)
   except Exception as mes:
      print("Something wrong with fit: ", mes)
      raise SystemExit

   model = fit_functions.exponential(fit.params,x)

   return model

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

   fig,ax = plt.subplots(figsize=(10,7))

   ax.scatter(year,mean_temp, marker='o', s=8,color='black')
   model = fit_linear(mean_temp,x=year)
   ax.plot(year,model,'k--',
               label='Mean temperature: %+.2f'%(model[-1]-model[0]))

   ax.fill_between(year,min_temp,mean_temp,
         color='xkcd:deep sky blue',alpha=0.4,linewidth=0)
#   plt.scatter(year,min_temp,  marker='o', s=8,color='blue')
   model = fit_linear(min_temp,x=year)
   ax.plot(year,model,'b--',
               label='Minimum temperature: %+.2f'%(model[-1]-model[0]))

   ax.fill_between(year,mean_temp,max_temp,
         color='xkcd:bright red',alpha=0.4,linewidth=0)
#   plt.scatter(year,max_temp,  marker='o', s=8,color='red')
   model = fit_linear(max_temp,x=year)
   ax.plot(year,model,'r--',
               label='Maximum temperature: %+.2f'%(model[-1]-model[0]))

   ax.legend(loc='lower right')

   plt.title(f"Temperature in {month_name[month]} in {city}")
   plt.xlabel("Year")
   plt.ylabel("Temperature (Celsius degrees)")
   plt.xlim([year_start-1,year_end+1])
   plt.ylim([y_min-1,y_max+1])

   ax.tick_params(labelright=True,right=True,which='both')
   ax.xaxis.set_major_locator(MultipleLocator(10))
   ax.xaxis.set_minor_locator(MultipleLocator(5))
   ax.yaxis.set_major_locator(MultipleLocator(5))
   ax.yaxis.set_minor_locator(MultipleLocator(1))
   
   fig.tight_layout()

   temp_evo = model[-1]-model[0]
#   temp_evo = (100*(model[-1]-model[0])).astype('int')/100

   return fig,temp_evo

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
   mean_temp,map = extract_temperatures(municipalities,"mean temperature",
                     year=year,month=month)
   min_temp,map  = extract_temperatures(municipalities,"minimum temperature",
                     year=year,month=month)
   max_temp,map  = extract_temperatures(municipalities,"maximum temperature",
                     year=year,month=month)


   d = {"year": year}
   for city in cities:
      idx = np.where(municipalities["NAMEFIN"]==city)
      d[f"mean_temp_{city}"] = [float(mean_temp[idx][0])]
      d[f"min_temp_{city}"]  = [float(min_temp[idx][0])]
      d[f"max_temp_{city}"]  = [float(max_temp[idx][0])]

   new_row      = pd.DataFrame(d)
   temperatures = pd.concat([temperatures, new_row], ignore_index=True)


temp_evo = np.array([])

for city in cities:
   fig,temp = plot_temperature_evolution(temperatures,city,
                        year_start=1961,year_end=2024,month=month)
   temp_evo = np.append(temp_evo,f"{temp:+.2}")
   fig.savefig(FIG_DIRECTORY / "Temperatures_evolution" /
               f"Temperature_{month_name[month]}_{city}.png")

try:
   df = geopandas.read_file(DATA_DIRECTORY /
         "finland_municipalities_temperature_evolution.gpkg"
         ).to_crs("EPSG:3067")
except:
   df = municipalities[["GML_ID","NAMEFIN","NAMESWE","geometry"]]
df[f"TEMP_EVO_{month_name[month]}"] = temp_evo
df.to_file(DATA_DIRECTORY /
      "finland_municipalities_temperature_evolution.gpkg")


