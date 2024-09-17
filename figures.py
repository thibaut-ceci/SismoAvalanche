import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import computations as cp

def save(path, tight_layout=True):
    if tight_layout:
        plt.tight_layout()
    plt.savefig(path)
    plt.show()


def linear_model(x, a, b):
    return a + b * x



def plot_fitted_curve(catalog, column1, column2, ax, number):
    # Logarithmic transformation of data
    constant = 0.01 # To avoid log(0) and negative values
    log_X = np.log(catalog[column1] + constant)
    log_Y = np.log(catalog[column2] + constant)

    # Fit the linear function to the logarithmically transformed data
    popt, pcov = curve_fit(linear_model, log_X, log_Y)

    # Convert the fitted parameters back to the scale of the original model
    a = np.exp(popt[0])
    b = popt[1]

    # Prepare data for plotting the fitted curve
    X_fit = np.linspace(min(catalog[column1])/300, max(catalog[column1])*500, 400)
    Y_fit = a * X_fit ** b

    # Plotting
    ax[number].scatter(catalog[column1], catalog[column2], label='Original Data')
    ax[number].plot(X_fit, Y_fit, 'r-', label=f'Fitted Curve: {column2} = {a:.2f} * {column1}^{b:.2f}')

    #plot confidence interval
    ax[number].set_xscale('log')
    ax[number].set_yscale('log')
    ax[number].set_xlabel(f"{column1}" + r" (m$^3$)", fontsize=14)
    ax[number].set_ylabel(f"{column2}" + r" (m$^3$)", fontsize=14)
    ax[number].legend()

    return log_X, log_Y, X_fit, Y_fit, popt, pcov, a, b


def plot_fitted_curve_without_incertainties(ESEC_original, catalog_without_incertainties, catalog_filtered, X_fit, Y_fit, Y_fit_lower, Y_fit_upper, Y, a, b, column1, column2, column3, column4, ax, number):
    # Plot the fitted curve
    ax[number].plot(X_fit, Y_fit, 'r-', label=f'Ajustement de courbe : DV = {a:.2f} * volume^{b:.2f}')

    # Plot the original data with incertainties and without incertainties
    ax[number].scatter(catalog_without_incertainties[column1], Y, color='blue', label='Données avec des incertitudes estimées')
    ax[number].errorbar(catalog_filtered[column1], catalog_filtered[column2], xerr=[catalog_filtered[column3], catalog_filtered[column4]], fmt='o', ecolor='gray', color='C0', capsize=2, capthick=0.75, elinewidth=.75)

    # Plot confidence bounds
    ax[number].plot(X_fit, Y_fit_lower, 'r--', alpha=0.5, label='Interval de confianceREFAIRELABEL à 95%')
    ax[number].plot(X_fit, Y_fit_upper, 'r--', alpha=0.5)

    # Compute new incertainties
    catalog_filtered = cp.compute_h_lower_upper(ESEC_original, column1, column3, column4, catalog_without_incertainties, Y, Y_fit_lower, Y_fit_upper, X_fit, catalog_filtered, ax, number)

    ax[number].set_xscale('log')
    ax[number].set_yscale('log')
    ax[number].set_xlabel(f"{column1}" + r" (m$^3$)", fontsize=14)
    ax[number].set_ylabel(f"{column2}" + r" (m$^3$)", fontsize=14)
    ax[number].legend()

    return catalog_filtered#, line_values


def scatter_plot(columnx, columny, labelx, labely):
    plt.figure(figsize=(10, 6))
    plt.scatter(columnx, columny, edgecolors='black')
    plt.xlabel(labelx)
    plt.xscale('log')
    plt.ylabel(labely)
    plt.show()








def plot__all_traces(stream_filtered, scale):
    _, ax = plt.subplots(figsize=(10, 15), dpi=150)
    for trace in stream_filtered:

        time = trace.times()
        time_max = trace.times().max()
        waveform = trace.data
        waveform = waveform - np.mean(waveform)
        waveform = waveform / np.max(np.abs(waveform))
        waveform *= scale
        distance = trace.stats.distance 
        ax.plot(time, waveform + distance, color="k", lw=0.2, rasterized=True)
        ax.set_xlabel("Temps (s)")
        ax.set_ylabel("Distance (km)")
        ax.set_ylim(0, 600)
        ax.set_xmargin(0)
        plt.xticks(np.arange(0, time_max, 10))
        plt.xticks(rotation=90, fontsize=6)

    return ax
