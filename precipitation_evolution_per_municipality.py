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
from geoanalysis_functions import extract_precipitations
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
            coeff={'value': 2, 'min': -100, 'max': 100},
            const={'value': np.mean(y), 'min': 0, 'max': 200})

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
def plot_precipitations_evolution(precipitations,city,
         year_start=1961,year_end=2024,month=1):

   # Read the temperature for the selected city
   year = np.array(precipitations["year"]).astype(int)
   rain = np.array(precipitations[f"rain_mm_{city}"]).astype(float)

   colnames = precipitations.columns
   y_max = precipitations.loc[:,colnames.str.contains("rain_")].max(axis=None)

   ###########################
   ### Add exponential fit ###
   ###########################

   fig,ax = plt.subplots(figsize=(10,7))

   ax.scatter(year,rain, marker='o', s=8,color='black')
   model = fit_linear(rain,x=year)
   ax.plot(year,model,'k--',
               label='Precipitations: %+.2f'%(model[-1]-model[0]))

   ax.legend(loc='lower right')

   plt.title(f"Precipitations in {month_name[month]} in {city}")
   plt.xlabel("Year")
   plt.ylabel("Precipitations (mm)")
   plt.xlim([year_start-1,year_end+1])
   plt.ylim([0,y_max+1])

   ax.tick_params(labelright=True,right=True,which='both')
   ax.xaxis.set_major_locator(MultipleLocator(10))
   ax.xaxis.set_minor_locator(MultipleLocator(5))
   ax.yaxis.set_major_locator(MultipleLocator(10))
   ax.yaxis.set_minor_locator(MultipleLocator(5))

   fig.tight_layout()

   rain_evo = model[-1]-model[0]
#   temp_evo = (100*(model[-1]-model[0])).astype('int')/100

   return fig,rain_evo

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
         "finland_municipalities_weather_evolution.gpkg"
         ).to_crs("EPSG:3067")
         #.to_crs("EPSG:4326")


# Initialisation temperature DataFrame for the selected month
precipitations = pd.DataFrame(columns=["year"]+[ f"rain_mm_{city}"
            for city in municipalities["NAMEFIN"] ])


for year in tqdm(range(1961,2025,1)):
   rain_per_km2,map = extract_precipitations(municipalities,
                     year=year,month=month)

   d = {"year": year}
   for idx in range(len(municipalities)):
      city = municipalities["NAMEFIN"][idx]
      d[f"rain_mm_{city}"] = [float(rain_per_km2[idx])]

   new_row      = pd.DataFrame(d)
   precipitations = pd.concat([precipitations, new_row], ignore_index=True)


rain_evo = np.array([])

for idx in range(len(municipalities)):
   city = municipalities["NAMEFIN"][idx]
   fig,rain = plot_precipitations_evolution(precipitations,city,
                        year_start=1961,year_end=2024,month=month)
   rain_evo = np.append(rain_evo,f"{rain:+.2}")
   fig.savefig(FIG_DIRECTORY / "Precipitations_evolution" /
               f"Precipitations_{month_name[month]}_{city}.png")

#municipalities[f"PRECIPITATION_EVO_{month_name[month]}"] = rain_evo
#municipalities.to_file(DATA_DIRECTORY /
#      "finland_municipalities_weather_evolution.gpkg")


