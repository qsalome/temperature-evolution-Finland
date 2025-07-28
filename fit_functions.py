import numpy as np
import lmfit


#--------------------------------------------------------------------
def linear(pars,x,data=None,eps=None):
    """
    Compute model, residual with and without errors for a linear function
   
   Parameters
   ----------
   pars: lmfit method
         parameters of the function
   x: 1D array
         data points of the independant variable
   data: 1D array
         data points on the dependant variable
   eps: 1D array
         errors on the dependant variable

   Returns
   -------
   1D array
         data points of the model or the residual with or without errors
    """
    parvals = pars.valuesdict()
    coeff = parvals['coeff']
    const = parvals['const']

    model = coeff * x + const

    if data is None:
        return model
    if eps is None:
        return (model-data)
    return (model-data)/eps

#--------------------------------------------------------------------
def exponential(pars,x,data=None,eps=None):
    """
    Compute model, residual with and without errors for an exponential function
   
   Parameters
   ----------
   pars: lmfit method
         parameters of the function
   x: 1D array
         data points of the independant variable
   data: 1D array
         data points on the dependant variable
   eps: 1D array
         errors on the dependant variable

   Returns
   -------
   1D array
         data points of the model or the residual with or without errors
    """
    parvals = pars.valuesdict()
    amp   = parvals['amp']
    coeff = parvals['coeff']
    const = parvals['const']

    model = amp * np.exp(coeff * (x-x[0])) + const

    if data is None:
        return model
    if eps is None:
        return (model-data)
    return (model-data)/eps

#--------------------------------------------------------------------


