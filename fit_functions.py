import numpy as np
import lmfit


#--------------------------------------------------------------------
def linear(pars,x,data=None,eps=None):
    """!
    Compute model, residual with and without errors
    
    @param pars: lmfit method : parameters of the function
    @param x: 1D array : velocity in channel unit
    @param data: 1D array : spectrum / Brightness Temperature
    @param eps: 1D array : errors
    
    @return model if data is None \n
    model - data if eps is None \n
    (model - data) / eps else
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
    """!
    Compute model, residual with and without errors
    
    @param pars: lmfit method : parameters of the function
    @param x: 1D array : velocity in channel unit
    @param data: 1D array : spectrum / Brightness Temperature
    @param eps: 1D array : errors
    
    @return model if data is None \n
    model - data if eps is None \n
    (model - data) / eps else
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


