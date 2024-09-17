"""Analyse seismic waveforms with ObsPy."""

import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from datetime import timedelta
import matplotlib.dates as mdates
import numpy as np
import obspy
import pandas as pd
import computations as cp
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import welch
import numpy as np
import computations as cp
import figures


def next2pow(x):
    return 1 if x == 0 else 2 ** (x - 1).bit_length()


def filtering(trace, freq_HP=9, freq_LP=0.5, max_percentage = 0.3):
    """
    DONNE UNE TRACE OU UN STREAM
    REND UNE TRACE OU UN STREAM FILTRe ET DETREND
    """
    #trace = trace.select(component="Z")
    trace = trace.detrend()
    trace = trace.filter("highpass", freq=freq_LP)  # filter pass haut
    trace = trace.filter("lowpass", freq=freq_HP)  # filter pass bas
    trace.filter("bandpass", freqmin=freq_LP, freqmax=freq_HP)
    trace = trace.taper(max_percentage=max_percentage, type="hann")

    return trace


def filter_stream(ESEC_avalanches, event_index, trace_index, freq_HP=9, freq_LP=0.5):
    event = ESEC_avalanches.loc[ESEC_avalanches["numero"] == event_index]
    stream = obspy.read(f"sismogrammes/{event_index:03d}.pickle")
    stream = stream.sort(keys=["distance"]) #trier le stream en fonction de la distance
    stream = stream.copy().select(component="Z")
    stream.remove_response(output="DISP", pre_filt=(0.01, 0.02, 20, 25), water_level=60)
    filtering(stream, freq_HP=freq_HP, freq_LP=freq_LP)

    trace = stream[trace_index]

    return event, stream, trace


def seuils(data_start_detection):
    data_rms = np.sqrt(np.mean(data_start_detection**2)) / 1.608571429 #2.608571429
    data_abs_median = 2 * np.median(np.abs(data_start_detection))

    # Threshold 1 : the noise
    threshold_1 = 3.399404762 * data_abs_median #
    data_above_threshold_1 = data_start_detection > threshold_1

    # Threshold 2 : the signal
    data_above_threshold_2 = data_above_threshold_1 & (np.abs(data_start_detection) > data_rms)

    data_rms = data_rms + 1.1 * threshold_1

    return threshold_1, data_rms, data_above_threshold_1, data_above_threshold_2


def methode_de_detection(trace, start_time_catalogue, end_time_catalogue, detection_yes_or_no, index):
    """
    DONNE UNE TRACE PRETE A LA DETECTION découpée par les temps de début et de fin du catalogue
    REND LA PARTIE DE LA TRACE DETECTE AVEC LES SEUILS
    """

    # Variables initiales
    time_raw = trace.times()
    distance = trace.stats.distance

    wavespeed = 6.5

    start_time = distance / wavespeed + start_time_catalogue
    end_time = distance / wavespeed + end_time_catalogue

    mask = (time_raw >= start_time) & (time_raw <= end_time)

    time_start_detection = time_raw[mask]
    data_start_detection = trace.data[mask]

    seuil_bas, seuil_haut, data_above_threshold_1, data_above_threshold_2 = seuils(data_start_detection)

    try:
        start_time = time_start_detection[data_above_threshold_1][0]
        end_time = time_start_detection[data_above_threshold_2][-1]

        mask = (time_start_detection >= start_time) & (time_start_detection <= end_time)

        trimmed_time = time_start_detection[mask]
        trimmed_data = data_start_detection[mask]

        detection_yes_or_no.append("TRUE")
        print("Detection on trace", index)

    except IndexError:
        print("No detection on trace", index)
        start_time = np.nan
        end_time = np.nan
        trimmed_time = np.nan
        trimmed_data = np.nan

        detection_yes_or_no.append("FALSE")

    duration = end_time - start_time

    if duration < 10:
        print("Durée de l'événement trop courte " + str(duration) + " secondes. Mauvaise détection on trace", index)
        start_time = np.nan
        end_time = np.nan
        trimmed_time = np.nan
        trimmed_data = np.nan

        detection_yes_or_no[-1] = "FALSE"

    if seuil_bas > seuil_haut:
        print("Seuil bruit trop élevé - no detection on trace", index)
        trimmed_time = np.nan 
        trimmed_data = np.nan 
        seuil_haut = np.nan 
        seuil_bas = np.nan 
        start_time = np.nan 
        end_time = np.nan 
        detection_yes_or_no[-1] = "FALSE"

    return trimmed_time, trimmed_data, seuil_haut, seuil_bas, start_time, end_time, detection_yes_or_no

#methode de detection totale
def methode_de_detection_totale(trace, ESEC_avalanches, event_index, add_start_time, add_end_time, detection_yes_or_no, index):
    #donner une trace
    #sort la partie de la trace pour la détection

    time_raw = trace.times()
    data_raw = trace.data

    distance = trace.stats.distance
    print("The distance of the trace is " + str(distance))

    start_time_catalogue = cp.conversion_du_temps_du_catalogue(trace, ESEC_avalanches["starttime"][event_index], add_start_time)
    end_time_catalogue = cp.conversion_du_temps_du_catalogue(trace, ESEC_avalanches["endtime"][event_index], add_end_time)

    wavespeed = 6.5

    start_time = distance / wavespeed + start_time_catalogue
    end_time = distance / wavespeed + end_time_catalogue

    mask = (time_raw >= start_time) & (time_raw <= end_time)

    time_start_detection = time_raw[mask]
    data_start_detection = trace.data[mask]

    #methode de detection
    trimmed_time, trimmed_data, seuil_haut, seuil_bas, start_time, end_time, detection_yes_or_no = methode_de_detection(trace, start_time_catalogue, end_time_catalogue, detection_yes_or_no, index)

    return time_start_detection, data_start_detection, trimmed_time, trimmed_data, time_raw, data_raw, seuil_haut, seuil_bas, detection_yes_or_no


def methode_de_detection_pour_un_evenement(stream, ESEC_avalanches, event_index, label1, label2):
    _, ax = plt.subplots(figsize=(10, 15), dpi=100)

    trimmed_time_starttime = []
    trimmed_time_endtime = []
    distance_all_trace = []
    detection_yes_or_no = []
    detection_stopped = False

    for index, trace in enumerate(stream):
        #Plot des traces brutes
        print('Work on trace', index)

        add_trace_color = False

        time_real, waveform = cp.plot_waveform_with_distance2(ESEC_avalanches, index, add_trace_color, trace.data, event_index, trace.times(), trace.stats.distance, ax)

        if not detection_stopped:
            #Detection method
            _, _, trimmed_time, trimmed_data, _, _, _, _, detection_yes_or_no = methode_de_detection_totale(trace, ESEC_avalanches, event_index, -30, 10, detection_yes_or_no, index)
        
            try:
                if len(detection_yes_or_no) >= 5 and all(d == "FALSE" for d in detection_yes_or_no[-5:]):
                    print("Signal non visible : 5 détections FALSE consécutives")
                    print("ARRET DE LA METHODE DE DETECTION !! LE SIGNAL N'EST PLUS VISIBLE")
                    detection_stopped = True
                    #stoppe la méthode de détection
                    trimmed_time_starttime.append(np.nan)
                    trimmed_time_endtime.append(np.nan)
                    distance_all_trace.append(trace.stats.distance)

                if not detection_stopped:
                    #Plot des traces détectées
                    add_trace_color = True
                    _, _, label1, label2 = cp.plot_waveform_with_distance(ESEC_avalanches, index, add_trace_color, trimmed_data, event_index, trimmed_time, trace.stats.distance, label1, label2, ax, color="blue", lw=0.6, rasterized=True)

                    trimmed_time_starttime.append(trimmed_time[0])
                    trimmed_time_endtime.append(trimmed_time[-1])
                    distance_all_trace.append(trace.stats.distance)

            except TypeError:
                print("Pas de détection !")
                trimmed_time_starttime.append(np.nan)
                trimmed_time_endtime.append(np.nan)
                distance_all_trace.append(trace.stats.distance)

        else:
            ax.plot(time_real, waveform + trace.stats.distance, color="k", lw=0.2, rasterized=True)
            trimmed_time_starttime.append(np.nan)
            trimmed_time_endtime.append(np.nan)
            distance_all_trace.append(trace.stats.distance)
        
        print("")

    ax.set_xlabel("Temps [s]", fontsize=14)
    ax.set_ylabel("Distance [km]", fontsize=14)
    # ax.set_xlim(np.min(time_real) + timedelta(seconds=410), np.max(time_real) - timedelta(seconds=410))
    ax.set_xlim(np.min(time_real) + timedelta(seconds=0), np.max(time_real) - timedelta(seconds=0))
    ax.set_ylim(0, 600)
    ax.margins(y=0)
    ax.tick_params(axis='both', which='major', labelsize=12)
    
    ax.legend(loc="upper right", fontsize=14, handlelength=4, handleheight=2, markerscale=1.5)

    formatter = mdates.DateFormatter('%H:%M:%S')
    ax.xaxis.set_major_formatter(formatter)
    
    figures.save(f"figures/event_detection/event_{event_index}.pdf")

    plt.show()

    return trimmed_time_starttime, trimmed_time_endtime, distance_all_trace, label1, label2


def methode_de_detection_pour_tous_les_evenements(ESEC_avalanches, stream, event_index, events_to_keep, label1, label2):
        trimmed_time_starttime = []
        trimmed_time_endtime = []
        distance_all_trace = []
        detection_yes_or_no = []

        detection_stopped = False

        for index, trace in enumerate(stream):
            print("Work on trace", index)

            if not detection_stopped:
                #Detection method
                _, _, trimmed_time, _, _, _, _, _, detection_yes_or_no = methode_de_detection_totale(trace, ESEC_avalanches, event_index, -30, 10, detection_yes_or_no, index)
            
                try:
                    if len(detection_yes_or_no) >= 5 and all(d == "FALSE" for d in detection_yes_or_no[-5:]):
                        print("Signal non visible : 5 détections FALSE consécutives")
                        print("ARRET DE LA METHODE DE DETECTION !! LE SIGNAL N'EST PLUS VISIBLE")
                        detection_stopped = True
                        #stoppe la méthode de détection
                        trimmed_time_starttime.append(np.nan)
                        trimmed_time_endtime.append(np.nan)
                        distance_all_trace.append(trace.stats.distance)

                    if not detection_stopped:
                        #Plot des traces détectées
                        trimmed_time_starttime.append(trimmed_time[0])
                        trimmed_time_endtime.append(trimmed_time[-1])
                        distance_all_trace.append(trace.stats.distance)

                except TypeError:
                    print("Pas de détection !")
                    trimmed_time_starttime.append(np.nan)
                    trimmed_time_endtime.append(np.nan)
                    distance_all_trace.append(trace.stats.distance)

            else:
                trimmed_time_starttime.append(np.nan)
                trimmed_time_endtime.append(np.nan)
                distance_all_trace.append(trace.stats.distance)
            
            print("")


        df = pd.DataFrame({'start_time': trimmed_time_starttime, 'end_time': trimmed_time_endtime, 'distance': distance_all_trace}) 
        df = df.sort_values('distance')
        df = df.reset_index(drop=True)
        df['detection'] = df['start_time'].apply(lambda x: False if np.isnan(x) else True)
        df['duration'] = df['end_time'] - df['start_time']

        if df['detection'].any():
            events_to_keep.append(event_index)  

        volume = [ESEC_avalanches["volume"][event_index]]

        for i in range(len(df["detection"])):
            if df["detection"][i] == True:
                if label1:
                    plt.scatter(volume, df["distance"][i], c='blue', marker='o', s=30, label="Detection")
                    label1 = False
                else:
                    plt.scatter(volume, df["distance"][i], c='blue', marker='o', s=30)
            else:
                if label2:
                    plt.scatter(volume, df["distance"][i], c='gray', marker='x', s=40, label="No detection")
                    label2 = False
                else:
                    plt.scatter(volume, df["distance"][i], c='gray', marker='x', s=40)
        
        # return label1, label2


def fit_line(frequencies, values):
    """ Ajuste une droite sur les données """
    model = LinearRegression() #tester d'autres modèles
    f = frequencies.reshape(-1, 1) 
    model.fit(f, values)
    return model.coef_[0], model.intercept_


def find_split_frequency(frequencies, values, min_freq=2.0, max_freq=10.0):
    """ Trouve la fréquence de séparation en détectant les discontinuités dans une plage spécifique """
    freq_mask = (frequencies >= min_freq) & (frequencies <= max_freq)
    filtered_frequencies = frequencies[freq_mask]
    filtered_values = values[freq_mask]
    
    if len(filtered_frequencies) < 2:
        return None

    derivative = np.diff(filtered_values) / np.diff(filtered_frequencies)
    split_index = np.argmax(np.abs(np.diff(derivative)))
    return filtered_frequencies[split_index + 1] if split_index < len(filtered_frequencies) - 1 else filtered_frequencies[-1]


def ajustement_de_segment(low_mask, frequencies_signal, psd_signal, ax, color='green', label="Ajustement bas", pltplot = True):
    if np.any(low_mask):
        low_freq = frequencies_signal[low_mask]
        low_psd = psd_signal[low_mask]
        low_slope, low_intercept = fit_line(np.log(low_freq), np.log(low_psd))

        if pltplot == True:
            plt.loglog(low_freq, np.exp(low_slope * np.log(low_freq) + low_intercept), color=color, label=label)
        else:
            ax.loglog(low_freq, np.exp(low_slope * np.log(low_freq) + low_intercept), color=color, label=label)

    return low_freq, low_slope, low_intercept, low_psd


def plot_spectre(trace, ESEC_avalanches, trimmed_data, trace_index, event_index, conserv_result=False):

    ax = plt
    curve_params = []

    # Welch
    segment_duration = 20
    noverlap = 12
    nperseg = int(segment_duration * trace.stats.sampling_rate)

    frequencies_signal, psd_signal = welch(trimmed_data, window='hamming', fs=trace.stats.sampling_rate, nperseg=nperseg, noverlap=noverlap)

    #Couper le spectre entre 0.5 et 9 Hz pour le filtrage
    mask = (frequencies_signal > 0.5) & (frequencies_signal < 9)
    frequencies_signal = frequencies_signal[mask]
    psd_signal = psd_signal[mask]

    split_freq = find_split_frequency(frequencies_signal, psd_signal, min_freq=1, max_freq=10)

    low_mask = frequencies_signal <= split_freq
    high_mask = frequencies_signal > split_freq

    #Ajustement des modeles
    low_freq, low_slope, low_intercept, low_psd = ajustement_de_segment(low_mask, frequencies_signal, psd_signal, ax, color='green', label="Modèle bas", pltplot = True)

    high_freq, high_slope, high_intercept, high_psd = ajustement_de_segment(high_mask, frequencies_signal, psd_signal, ax, color='blue', label="Modèle haut", pltplot = True)

    if conserv_result == True:
        curve_params.append({
                'Event Index': event_index,
                'Fréquence coin': split_freq,
                'Slope basse frequence': low_slope,
                'Intercept basse frequence': low_intercept,
                'First PSD basse frequence': low_psd[0],
                'PSD requence coin': low_psd[-1],
                'Slope haute frequence': high_slope,
                'Intercept haute frequence': high_intercept,
                'Last PSD haute frequence': high_psd[-1], 
                'numero1': ESEC_avalanches["numero"][event_index]
            })
        
    plt.loglog(frequencies_signal, psd_signal, color="C1", label="Spectre du signal sismique détecté")
    plt.legend()
    plt.margins(x=0)
    plt.xscale("log")
    plt.xlabel('Fréquences (Hz)')
    plt.ylabel(r'Densité Spectrale de Puissance du déplacement ($\mathrm{\frac{m^{2}}{Hz}}$)')

    if conserv_result == True:
        df = pd.DataFrame(curve_params)
        df.to_csv(f'features/1_fitting/data/curve_parameters_{event_index}.csv', index=False)

    figures.save(f"features/1_fitting/pictures/fitting_on_trace_{trace_index}_in_event_{event_index}.pdf")
    plt.show()

    return curve_params


def plot_detected_trace(time_raw, data_raw, time_start_detection, data_start_detection, trimmed_time, trimmed_data):
    _, ax = plt.subplots(figsize=(10, 5))
    ax.plot(time_raw, data_raw, label = "Signal")
    ax.plot(time_start_detection, data_start_detection, "C1", label="Signal for detection")
    ax.plot(trimmed_time, trimmed_data, "C2", label="Signal détecté")

    ax.set_xlabel("Temps [s]", fontsize = 14)
    ax.set_ylabel("Déplacement [m]", fontsize = 14)
    ax.set_xlim(time_raw[0], time_raw[-1])
    ax.legend(fontsize = 14)
    plt.show()
