"""ESEC catalog management.

This module contains functions to read and filter the catalog obtained by the
ESEC project. A REFAIIIIIIIIIIIIIIIIIIIIIIIRE
"""

import pickle

import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from datetime import timedelta
import matplotlib.dates as mdates
import numpy as np
from obspy.clients.fdsn import Client
import obspy
import pandas as pd
from scipy.signal import welch, spectrogram
from scipy.optimize import curve_fit
import computations as cp
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats
from scipy.signal import spectrogram
import computations as cp
import analysis
import figures


def spectres(trace):
    #Paramètres des traces
    distance = trace.stats.distance
    fs = trace.stats.sampling_rate
    data = trace.data
    
    #Calcul de l'énergie des traces
    window_len = 20
    nperseg = int(window_len * fs)

    frequencies, times, Sxx = spectrogram(data, fs=fs, nperseg=nperseg, noverlap=nperseg * 0.9, nfft=analysis.next2pow(nperseg), scaling="density", mode="psd")

    frequency_band = (0.1, 2) ########## freq high and low to modify
    frequency_band_indices = np.where((frequencies > frequency_band[0]) & (frequencies < frequency_band[1]))[0]
    energy = Sxx[frequency_band_indices, :].sum(axis=0)
    energy_max_index = np.argmax(energy) #index
    time_at_max_energy = times[energy_max_index] #index

    return distance, times, energy, time_at_max_energy, energy_max_index


def seuil(times, energy):
    #Seuil énergétique - détection
    seuil_energy = energy.mean() / 2 #seuil pour sélectionner que la parabole de l'énergie
    first_index = np.argmax(energy > seuil_energy) #premier index de la ligne de détection
    energy_above_seuil_energy = energy > seuil_energy
    last_index = np.where(energy_above_seuil_energy)[0][-1]
    sub_times = times[first_index:last_index+1]
    sub_energy = energy[first_index:last_index+1]
    energy_max_index = np.argmax(sub_energy)

    return sub_times, sub_energy, energy_max_index, first_index, last_index, seuil_energy


def plot(times, energy, distance, seuil_energy, sub_times, sub_energy, i):
    # plt.fill_between(times, energy, zorder=i+1)
    plt.plot(times, energy, label="distance " + str(np.round(distance)) + " km", zorder=i, alpha=0.7) #plot de l'énergie spectrale entière
    plt.plot([times[0], times[-1]], [seuil_energy, seuil_energy], 'g--', label="Seuil de détection")
    plt.plot(sub_times, sub_energy, color="red") #plot de l'énergie spectrale entière (JUSTE LA PARABOLE)
    
    #plot parameters
    plt.xlabel("Temps (s)")
    plt.ylabel("Energy (dB)")
    plt.legend(fontsize=9, bbox_to_anchor=(1.01, 1), loc='upper left')


def compute(ESEC_avalanches, stream, event_index):
    plt.figure()

    distance_all = []
    pente_all = []
    skewness_all = []
    aire_all = []
    starttime_all = []
    endtime_all = []
    energie_trace = []

    for i, trace in enumerate(stream):
        distance, times, energy, time_at_max_energy, energy_max_index = spectres(trace)

        sub_times, sub_energy, energy_max_index, first_index, last_index, seuil_energy = seuil(times, energy)

        plot(times, energy, distance, seuil_energy, sub_times, sub_energy, i)

        starttime_all.append(times[first_index])
        endtime_all.append(times[last_index])
        distance_all.append(distance)
        energie_trace.append(np.max(sub_energy))
        aire_all.append(np.trapz(sub_energy, sub_times))  #calculer l'aire sous la courbe
        skewness_all.append(scipy.stats.skew(sub_energy)) #Calcul de la skewness
        pente_all.append((energy[energy_max_index] - energy[first_index]) / (time_at_max_energy - times[first_index])) #Calcul de la pente

    plt.savefig(f'features/2_energie/pictures/energy_{event_index}.png', bbox_inches='tight')
    plt.show()

    linespace_stream = np.linspace(0, len(stream)-1, len(stream))
    event_index_all = np.linspace(event_index, event_index, len(linespace_stream))

    endtime_all = np.array(endtime_all)
    starttime_all = np.array(starttime_all)

    differences = endtime_all - starttime_all

    #'trace': linespace_stream,
    dataframe_event = pd.DataFrame({'event_index': event_index_all, 'distance': distance_all, 'asymétrie': skewness_all, 
                                    'aire': aire_all, 'energie_max_trace': energie_trace, 'Impulsion': pente_all, 'starttime': starttime_all, 'endtime': endtime_all, 
                                    'durée': differences, 'numero2': ESEC_avalanches["numero"][event_index]})
    dataframe_event.to_csv(f'features/2_energie/data/dataframe_event_{event_index}.csv', index=False)




