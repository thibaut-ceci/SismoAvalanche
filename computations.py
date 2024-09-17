import datetime as dt
from scipy.signal import spectrogram
import numpy as np
import obspy
import pandas as pd
from scipy.stats import t
import datetime
from datetime import timedelta
from matplotlib import lines as mlines
import matplotlib.pyplot as plt
import glob
import os

import analysis

from tqdm.notebook import tqdm
tqdm.pandas()



def compute_waveform_with_distance(trace, scale):
    waveform = trace - np.mean(trace)
    waveform /= np.max(np.abs(waveform))
    waveform *= scale

    return waveform


def plot_waveform_with_distance2(ESEC_avalanches, trace, add_trace_color, trace_data, event_index, time, distance, ax, color = "k", lw = 0.2, rasterized=True):
    start_time_real = datetime.datetime.strptime(ESEC_avalanches["starttime"][event_index], "%Y_%m_%d %H%M%S")
    time_real = [start_time_real + timedelta(seconds=t) for t in time]
    waveform = compute_waveform_with_distance(trace_data, 10)
    ax.plot(time_real, waveform + distance, color=color, lw=lw, rasterized=rasterized)

    # if add_trace_color == True:
    #     if trace == 14:
    #         ax.plot(time_real, waveform + distance, color="orange", lw=lw, rasterized=rasterized)
    #     if trace == 76:
    #         ax.plot(time_real, waveform + distance, color="orange", lw=lw, rasterized=rasterized)

    return time_real, waveform

def plot_waveform_with_distance(ESEC_avalanches, trace, add_trace_color, trace_data, event_index, time, distance, label1, label2, ax, color = "k", lw = 0.2, rasterized=True):
    start_time_real = datetime.datetime.strptime(ESEC_avalanches["starttime"][event_index], "%Y_%m_%d %H%M%S")
    time_real = [start_time_real + timedelta(seconds=t) for t in time]
    waveform = compute_waveform_with_distance(trace_data, 10)

    if label1 == True:
        ax.plot(time_real, waveform + distance, color=color, lw=lw, rasterized=rasterized, label = "Signal détecté")
        label1 = False
    else:
        ax.plot(time_real, waveform + distance, color=color, lw=lw, rasterized=rasterized)

    if label2 == True:
        if trace == 14:
            ax.plot(time_real, waveform + distance, color="orange", lw=lw, rasterized=rasterized, label = "Fausse détection")
        if trace == 76:
            ax.plot(time_real, waveform + distance, color="orange", lw=lw, rasterized=rasterized)
            label2 = False

    return time_real, waveform, label1, label2


def compute_number_of_stations(ESEC_avalanches, print = False):
    for numero_event in tqdm(ESEC_avalanches["numero"], total=len(ESEC_avalanches)):
        inventory = ESEC_avalanches["inventory"][numero_event]

        num_stations = []

        if len(inventory) > 0: #len(inventory) give network len
            for i in range(len(inventory.networks)):
                num_stations.append(len(inventory.networks[i].stations))

            if print == True:
                print("The number of stations is", sum(num_stations), "for event", str(numero_event))

            ESEC_avalanches.at[numero_event, 'Number of stations'] = sum(num_stations)

        else:
            print("No stations for event", str(numero_event))
            ESEC_avalanches.at[numero_event, 'Number of stations'] = 0

    return ESEC_avalanches


def remove_outliers_in_catalog(catalog, catalog_column, subset):
    # Outlier detection via interquartile range (IQR)
    Q1 = catalog_column.quantile(0.25)
    Q3 = catalog_column.quantile(0.75)
    IQR = Q3 - Q1
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR

    print(f"Lower bound for {subset[1]}: {lower_bound}, Upper bound for {subset[1]}: {upper_bound}")

    outliers = catalog[(catalog_column < lower_bound) | (catalog_column > upper_bound)]

    print(f"Number of outliers in {subset[1]}:", outliers.shape[0])

    ESEC_avalanches = catalog.dropna(subset=subset)
    ESEC_avalanches_filtered = ESEC_avalanches[(catalog_column > lower_bound) & (catalog_column < upper_bound)]

    return ESEC_avalanches_filtered


def confidence_bounds(X_fit, log_X, popt, pcov):
    # Calculate the standard error (square root of the diagonals of the covariance matrix)
    perr = np.sqrt(np.diag(pcov))

    # Calculate the t statistic for a 90% confidence interval assuming a large sample size
    t_val = t.ppf(0.95, len(log_X) - len(popt))  # 0.975 for 95% CI, adjusting df for number of parameters
    # t_val = quantile 0.90

    # Calculate confidence interval bounds for the parameters
    a_conf_interval = np.exp(popt[0] + t_val * perr[0] * np.array([-1, 1])) # Pas compris :( ###############################################
    b_conf_interval = popt[1] + t_val * perr[1] * np.array([-1, 1]) # Pas compris :( ###############################################

    # Upper and lower confidence bounds of the fit
    Y_fit_lower = a_conf_interval[0] * X_fit ** b_conf_interval[0]
    Y_fit_upper = a_conf_interval[1] * X_fit ** b_conf_interval[1]

    return Y_fit_lower, Y_fit_upper


def compute_h_lower_upper(ESEC_original, X, low, high, catalog_without_incertainties, points_y, DX_fit_lower, DX_fit_upper, volumes_fit, catalog_filtered, ax, number):
    # List to store the first and last values of each line
    #line_values = []

    # Create a new column in dfi for the determined uncertainty
    catalog_filtered['determined_uncertainty'] = np.nan

    # Plot horizontal lines connecting each point with DV_fit_lower and DV_fit_upper
    for index, (h, y) in enumerate(zip(catalog_without_incertainties[X], points_y)):
        h_lower = np.interp(y, DX_fit_lower, volumes_fit)
        h_upper = np.interp(y, DX_fit_upper, volumes_fit)
        h_upper = h_upper
        h_lower = h_lower
        
        if index == 0:
            ax[number].plot([h_lower, h_upper], [y, y], 'b', alpha=0.2, label="Incertitudes estimées")
        else:
            ax[number].plot([h_lower, h_upper], [y, y], 'b', alpha=0.2)
        
        # Store the first and last values of the line
        #line_values.append(((h_lower, y), (h_upper, y)))
        
        # Get the original index in the "avalanche" DataFrame
        original_index = catalog_without_incertainties.index[index]
        
        # Assign the uncertainty values to the corresponding values in the "avalanche" DataFrame
        ESEC_original.at[original_index, low] = h_lower
        ESEC_original.at[original_index, high] = h_upper

        # Calculate the determined uncertainty and assign it to the new column
        ESEC_original.at[original_index, 'determined_uncertainty'] = "True"

    return ESEC_original#, line_values


def conversion_du_temps_du_catalogue(trace, start_time_string, add_time):
    start_time = dt.datetime.strptime(start_time_string, '%Y_%m_%d %H%M%S')
    start_time_seconds = (start_time - trace.stats.starttime.datetime).total_seconds()
    start_time_event_shifted = start_time_seconds + add_time

    return start_time_event_shifted


def create_detection_dataframe(ESEC_avalanches, event_index, trimmed_time_starttime, trimmed_time_endtime, distance_all_trace):
    df = pd.DataFrame({'start_time': trimmed_time_starttime, 'end_time': trimmed_time_endtime, 'distance': distance_all_trace})
    df = df.sort_values('distance') #sort by distance
    df = df.reset_index(drop=True) #reset index
    df['detection'] = df['start_time'].apply(lambda x: False if np.isnan(x) else True) #add true if detection is possible. add falsi is detection is not possible
    df['duration'] = df['end_time'] - df['start_time']

    volume = [ESEC_avalanches["volume"][event_index]]

    plt.figure(figsize=(10, 10), dpi=300)

    count_detection = 0
    count_non_detection = 0

    for i in range(len(df["detection"])):
        if df["detection"][i] == True:
            plt.scatter(volume, df["distance"][i], c='blue', marker='o', s=30)
            count_detection += 1
        else:
            plt.scatter(volume, df["distance"][i], c='red', marker='x', s=40)
            count_non_detection += 1

    detection_handle = mlines.Line2D([], [], color='blue', marker='o', linestyle='None', markersize=3, label='détection (' + str(count_detection) + ')')
    no_detection_handle = mlines.Line2D([], [], color='red', marker='x', linestyle='None', markersize=3, label='Mauvaise/absence de détection (' + str(count_non_detection) + ')')

    plt.legend(handles=[detection_handle, no_detection_handle])

    plt.xlabel(r"Volume ($m^3$)")
    plt.ylabel("Distance de la station à l'avalanche (km)")
    plt.ylim(0, 600)


def moyenne_glissante(f, shanon_index, ax, window_size = 40):
    moving_average = np.convolve(shanon_index, np.ones(window_size)/window_size, mode='valid')
    diff = len(f[:len(shanon_index)]) - len(moving_average)
    moving_average_padded = np.append(moving_average, [np.nan]*diff)
    ax[2].plot(f[:len(shanon_index)], moving_average_padded, c="red", label="Moyenne glissante")


def merge_dataframes(dossier = "features/1_fitting/data", name = "curve_parameters*.csv", area_to_save = 'features/1_fitting/data/fitting_df.csv'):
    fichiers_csv = glob.glob(os.path.join(dossier, name))

    dataframes = [pd.read_csv(fichier) for fichier in fichiers_csv]
    energy_2 = pd.concat(dataframes, ignore_index=True)

    energy_2.to_csv(area_to_save, index=False)

    return energy_2





##########################


# def find_the_closest_trace_of_the_event(stream_filtered):
#     distance = 0
#     closest_trace = None
#     for i, trace in enumerate(stream_filtered):
#         if trace.stats.distance < distance:
#             distance = trace.stats.distance
#             closest_trace = trace
#     return distance, closest_trace


# def find_the_closest_trace_of_X_km(stream_filtered, distance):
#     target_distance = distance
#     min_difference = float('inf')
#     closest_trace = None
#     for i, trace in enumerate(stream_filtered):
#         difference = abs(trace.stats.distance - target_distance)
#         if difference < min_difference:
#             min_difference = difference
#             closest_trace = trace
#     return closest_trace.stats.distance, closest_trace

# def find_the_closest_trace_and_index_with_distance(stream_filtered, distance):
#     min_difference = float('inf')
#     closest_trace_index = None

#     for i, trace in enumerate(stream_filtered):
#         current_difference = abs(trace.stats.distance - distance)
#         if current_difference < min_difference:
#             min_difference = current_difference
#             closest_trace_index = i
#             distance = trace.stats.distance

#     return closest_trace_index, stream_filtered[closest_trace_index], distance