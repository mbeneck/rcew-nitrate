import pandas as pd
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt

def plot_scan_fp(df, measurement, **kwargs):
    df.T.drop('Date/Time').plot(y=measurement, xlabel = 'Wavelength (nm)', ylabel= 'Absorbance', **kwargs)

def import_scan_fp(path):
    df = pd.read_csv(path, sep='\t', header=1).drop('Status_0', axis=1)
    df['Date/Time'] = pd.to_datetime(df['Date/Time'])
    #df.set_index('Date/Time', inplace=True)
    return df

def apply_calibrations(df, calibrations, calibrations_kwargs):
    for calibration in calibrations:
        df = calibration(df, **calibrations_kwargs[calibration])
    return df