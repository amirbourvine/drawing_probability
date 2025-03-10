import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

"""
This script is fitting various functions to the slope for ratio curve.
first run slope_for_ration.py to create slope_to_ratio.npz, that this script uses.
"""

def sqrt_func(x, a, b):
    return a * np.sqrt(x) + b

def log_func(x, a, b):
    return a * np.log(x) + b

def power_func(x, a, b, c):
    return a * x**b + c

def fit_best_function(x, y):
    """
    Fits different functions to (x, y) data and returns the best fit parameters.
    
    :param x: Array-like, independent variable data.
    :param y: Array-like, dependent variable data.
    :return: Dictionary with best-fit parameters for each function.
    """
    fits = {}
    
    try:
        popt_sqrt, _ = curve_fit(sqrt_func, x, y, maxfev=5000)
        fits['sqrt'] = popt_sqrt
    except:
        fits['sqrt'] = None
    
    try:
        popt_log, _ = curve_fit(log_func, x, y, maxfev=5000)
        fits['log'] = popt_log
    except:
        fits['log'] = None
    
    try:
        popt_power, _ = curve_fit(power_func, x, y, maxfev=5000)
        fits['power'] = popt_power
    except:
        fits['power'] = None
    
    return fits

def plot_fits(x, y, fits):
    """
    Plots the original data along with the fitted functions.
    
    :param x: Array-like, independent variable data.
    :param y: Array-like, dependent variable data.
    :param fits: Dictionary containing the fitted parameters for each function.
    """
    plt.figure(figsize=(8, 6))
    plt.scatter(x, y, label='Original Data', color='black', alpha=0.5)
    
    x_fit = np.linspace(min(x), max(x), 1000)
    
    if fits['sqrt'] is not None:
        plt.plot(x_fit, sqrt_func(x_fit, *fits['sqrt']), label='Square Root Fit', color='blue')
    if fits['log'] is not None:
        plt.plot(x_fit, log_func(x_fit, *fits['log']), label='Log Fit', color='red')
    if fits['power'] is not None:
        plt.plot(x_fit, power_func(x_fit, *fits['power']), label='Power Fit', color='green')
    
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Curve Fitting')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    data = np.load("slope_to_ratio.npz")
    ratios, slopes = data['x'], data['y']

    fits = fit_best_function(ratios, slopes)

    print(fits)

    plot_fits(ratios, slopes, fits)





    
