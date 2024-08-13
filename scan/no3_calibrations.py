import pandas as pd
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt

def second_derivative(df, output_calibrated = True, **kwargs):
    window_size = kwargs.get('window_size', 5)
    poly_order = kwargs.get('poly_order',2)
    wavelength = kwargs.get('wavelength', "227.50")
    slope = kwargs.get('slope', 1)
    intercept = kwargs.get('intercept', 0)
    prefix = kwargs.get('prefix', '')
    
    df_2d = df.drop('Date/Time', axis=1).transform(lambda x: savgol_filter(x, window_size, poly_order, deriv=2),axis=1)
    df[prefix+"second_derivative"]  = df_2d[wavelength]
    
    if(output_calibrated):
        df[prefix+'second_derivative_no3'] = df[prefix+"second_derivative"]*slope+intercept
    
    return df

def two_wavelength(df, output_calibrated=True, **kwargs):
    wavelength = kwargs.get('wavelength', '220.00')
    cdom_wavelength = kwargs.get('cdom_wavelength', '275.00')
    alpha = kwargs.get('alpha', 2)
    slope = kwargs.get('slope', 1)
    intercept = kwargs.get('intercept', 0)
    prefix = kwargs.get('prefix', '')
    
    df[prefix+'two_wavelength'] = df[wavelength] - alpha*df[cdom_wavelength]
    
    if(output_calibrated):
        df[prefix+'two_wavelength_no3'] = df[prefix+'two_wavelength']*slope+intercept
        
    return df

def one_wavelength(df, output_calibrated=True, **kwargs):
    wavelength = kwargs.get('wavelength', '220.00')
    slope = kwargs.get('slope', 1)
    intercept = kwargs.get('intercept', 0)
    prefix = kwargs.get('prefix', '')

    df[prefix+'one_wavelength'] = df[wavelength]
    
    if(output_calibrated):
        df[prefix+'one_wavelength_no3'] = df[prefix+'one_wavelength']*slope+intercept
        
    return df
