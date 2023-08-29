#%%
import calc.calc as calc
import importlib
importlib.reload(calc)
import numpy as np

#Set the concentration of the TE ions in mM in the array concTE:
concTE = np.array([1, 5])

#Dictionary with settings for the calculation:

settings = {
    'with_BME': False, #True if BME is present
    'c_BME_percent': 0, #BME concentration in percent
    'save_to_file': True, #True if the results should be saved to a file
    'showAllSol': False, #True if all solutions should be shown
    'showEq': False, #True if the equilibrium equations should be shown
    'showConc': False, #True if the concentrations should be shown
    'showCoeffs': False, #True if the coefficients should be shown
    'basepath' : '', #Path to the folder where the results should be saved
}

calc.calc_ionic_strength_TE(concTE, settings)
#%%