"""Analyse seismic waveforms with ObsPy."""

from datetime import timedelta
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import obspy
import pandas as pd
from scipy.signal import welch
from sklearn.linear_model import LinearRegression

import computations as cp
import figures


def filter_stream(ESEC_avalanches, event_index, trace_index, freq_HP=9, freq_LP=0.5, max_percentage=0.3):
    """
    Filter a stream.

    Parameters
    ----------
    ESEC_avalanches : pandas.DataFrame
        The ESEC
    event_index : int
        The index of the event to filter.
    trace_index : int
        The index of the trace to retrieve and filter.
    freq_HP : float
        The high-pass filter frequency.
    freq_LP : float
        The low-pass filter frequency.
    max_percentage : float
        The maximum percentage for tapering. 

    Returns
    -------
    event : pandas.DataFrame
        The event data corresponding to the event_index.
    stream : obspy.Stream
        The filtered seismic stream.
    trace : obspy.Trace
        The specific trace filtered from the stream.
    """

    ## Retrieve the event data
    event = ESEC_avalanches.loc[ESEC_avalanches["numero"] == event_index]

    ## Read and prepare the stream
    stream = obspy.read(f"sismogrammes/{event_index:03d}.pickle")
    stream = stream.sort(keys=["distance"])
    
    ## Remove response, detrend and filter the stream
    stream = stream.select(component="Z")
    stream.remove_response(output="DISP", pre_filt=(0.01, 0.02, 20, 25), water_level=60)
    stream = stream.detrend()
    stream = stream.filter("highpass", freq=freq_LP)                     # High-pass filter
    stream = stream.filter("lowpass", freq=freq_HP)                      # Low-pass filter
    stream = stream.filter("bandpass", freqmin=freq_LP, freqmax=freq_HP) # Band-pass filter
    stream = stream.taper(max_percentage=max_percentage, type="hann")    # Taper

    ## Select one trace
    trace = stream[trace_index]

    return event, stream, trace


def thresholds(data_start_detection):
    """
    Calculate noise and signal thresholds based on the seismic signal data.

    Parameters
    ----------
    data_start_detection : numpy.ndarray
        Seismic signal data.

    Returns
    -------
    threshold_1 : float
        The calculated noise threshold.
    adjusted_rms : float
        The adjusted RMS value after applying the noise threshold.
    data_above_threshold_1 : numpy.ndarray
        Boolean array indicating which data points exceed the noise threshold.
    data_above_threshold_2 : numpy.ndarray
        Boolean array indicating which data points exceed both the noise and signal thresholds.
    """

    ## Calculate the RMS of the data
    data_rms = np.sqrt(np.mean(data_start_detection**2)) / 1.608571429

    ## Calculate the median of the absolute values of the data
    data_abs_median = 2 * np.median(np.abs(data_start_detection))

    ## Threshold 1 : the noise
    threshold_1 = 3.399404762 * data_abs_median
    data_above_threshold_1 = data_start_detection > threshold_1

    ## Threshold 2 : the seismic signal
    data_above_threshold_2 = data_above_threshold_1 & (np.abs(data_start_detection) > data_rms)

    ## Adjusted RMS value
    data_rms = data_rms + 1.1 * threshold_1

    return threshold_1, data_rms, data_above_threshold_1, data_above_threshold_2


def detection_on_one_trace(trace, ESEC_avalanches, event_index, trace_index, add_start_time, add_end_time, detection_yes_or_no):
    """
    Prepares a trace for detection by trimming it according to the start and end times from the catalog.
    Returns the detected part of the trace with thresholds.

    Parameters
    ----------
    trace : obspy.Trace
        The trace to process for detection.
    ESEC_avalanches : pandas.DataFrame
        The ESEC.
    event_index : int
        The index of the event in the DataFrame.
    trace_index : int
        The index of the trace being processed.
    add_start_time : float
        Additional time to add to the start time.
    add_end_time : float
        Additional time to add to the end time.
    detection_yes_or_no : 
    
    Returns
    -------
    time_start_detection : numpy.ndarray
        Time for the detection method.
    data_start_detection : numpy.ndarray
        Data for the detection method.
    trimmed_time : numpy.ndarray
        Detected time.
    trimmed_data : numpy.ndarray
        Detected data from the seismic signal.
    time_raw : numpy.ndarray
        Original time of the trace.
    data_raw : numpy.ndarray
        Original data of the trace.
    upper_threshold : float
        Upper threshold for detection.
    lower_threshold : float
        Lower threshold for detection.
    """

    ### Step 1 : Trim the trace using the ESEC start and end times to refine the detection area

    ## Extract information from the trace
    time_raw = trace.times()
    data_raw = trace.data
    distance = trace.stats.distance
    print("The distance of the trace is " + str(distance))

    ## Extract the start and end times from ESEC
    start_time_catalogue = cp.conversion_du_temps_du_catalogue(trace, ESEC_avalanches["starttime"][event_index], add_start_time)
    end_time_catalogue = cp.conversion_du_temps_du_catalogue(trace, ESEC_avalanches["endtime"][event_index], add_end_time)

    ## Calculate adjusted times based on wave speed
    wavespeed = 6.5
    start_time = distance / wavespeed + start_time_catalogue
    end_time = distance / wavespeed + end_time_catalogue

    ## Create a mask to trim the trace
    mask = (time_raw >= start_time) & (time_raw <= end_time)
    time_start_detection = time_raw[mask]
    data_start_detection = trace.data[mask]


    ### Step 2 : The detection with the thresholds

    ## Compute thresholds
    lower_threshold, upper_threshold, data_above_threshold_1, data_above_threshold_2 = thresholds(data_start_detection)

    try:
        ## Compute start and end times of the trace
        start_time = time_start_detection[data_above_threshold_1][0]
        end_time = time_start_detection[data_above_threshold_2][-1]

        ## Trim the trace
        mask = (time_start_detection >= start_time) & (time_start_detection <= end_time)
        trimmed_time = time_start_detection[mask]
        trimmed_data = data_start_detection[mask]
        detection_yes_or_no.append("TRUE")

        print("Detection on trace", trace_index)

    except IndexError:
        ## If the detection is not possible, an error occurs.
        print("No detection on trace", trace_index)
        start_time = end_time = trimmed_time = trimmed_data = np.nan
        detection_yes_or_no.append("FALSE")


    ### Step 3 : Add conditions to improve the detection method

    ## Compute the duration of the avalanche
    duration = end_time - start_time

    ## Add condition if the duration is too small
    if duration < 10:
        print(f"Event duration too short: {str(duration)} seconds. Bad detection on trace {trace_index}.")
        start_time = end_time = trimmed_time = trimmed_data = np.nan
        detection_yes_or_no[-1] = "FALSE"

    ## Add condition if the lower threshold is upper the upper threshold
    if lower_threshold > upper_threshold:
        print(f"Noise threshold too high - no detection on trace {trace_index}.")
        trimmed_time = trimmed_data = upper_threshold = lower_threshold = start_time = end_time = np.nan
        detection_yes_or_no[-1] = "FALSE"

    return time_start_detection, data_start_detection, trimmed_time, trimmed_data, time_raw, data_raw, upper_threshold, lower_threshold, detection_yes_or_no


def detected_method_for_one_event(stream, ESEC_avalanches, event_index, label1):
    """
    Perform detection method on seismic traces for one event.

    Parameters
    ----------
    stream : obspy.Stream
        The stream containing traces to be analyzed.
    ESEC_avalanches : pandas.DataFrame
        The ESEC.
    event_index : int
        Index of the event to be analyzed.
    label1 : str
        The label. If true, the legend will appears.

    Returns
    -------
    trimmed_time_starttime : list
        List of start times of detected signals.
    trimmed_time_endtime : list
        List of end times of detected signals.
    distance_all_trace : list
        List of distances for each trace.
    label1 : str
        The updated label.
    """

    ## Initialize a figure for plotting and lists to store detection results
    _, ax = plt.subplots(figsize=(10, 15), dpi=100)

    trimmed_time_starttime = []
    trimmed_time_endtime = []
    distance_all_trace = []
    detection_yes_or_no = []
    detection_stopped = False

    ## Loop over each trace in the stream
    for index, trace in enumerate(stream):
        print('Work on trace', index)

        ## Plot the raw waveform
        time_real, waveform = figures.plot_waveform_with_distance_signal(ESEC_avalanches, trace.data, event_index, trace.times(), trace.stats.distance, ax)

        if not detection_stopped:
            ## The detection method
            _, _, trimmed_time, trimmed_data, _, _, _, _, detection_yes_or_no = detection_on_one_trace(trace, ESEC_avalanches, event_index, index, -30, 10, detection_yes_or_no)

            try:
                ## Check if five consecutive detections are "FALSE". If true, the detection has stopped
                if len(detection_yes_or_no) >= 5 and all(d == "FALSE" for d in detection_yes_or_no[-5:]):
                    print("Detection stopped: 5 consecutive FALSE detections !")
                    detection_stopped = True

                    ## Mark the trace as undetected
                    trimmed_time_starttime.append(np.nan)
                    trimmed_time_endtime.append(np.nan)
                    distance_all_trace.append(trace.stats.distance)

                ## Plot the detected signal if detection has not stopped
                if not detection_stopped:
                    _, _, label1 = figures.plot_waveform_with_distance_detected_signal(ESEC_avalanches, trimmed_data, event_index, trimmed_time, trace.stats.distance, label1, ax, color="blue", lw=0.6)

                    ## Store results
                    trimmed_time_starttime.append(trimmed_time[0])
                    trimmed_time_endtime.append(trimmed_time[-1])
                    distance_all_trace.append(trace.stats.distance)


            ## If it is not possible to trace the detected seismic signal, it means that the method failed to detect a seismic signal so an error will be generated.
            except TypeError:
                print("No detection !")
                trimmed_time_starttime.append(np.nan)
                trimmed_time_endtime.append(np.nan)
                distance_all_trace.append(trace.stats.distance)

        ## If detection has stopped, plot only the raw waveform without detection
        else:
            ax.plot(time_real, waveform + trace.stats.distance, color="k", lw=0.2)
            trimmed_time_starttime.append(np.nan)
            trimmed_time_endtime.append(np.nan)
            distance_all_trace.append(trace.stats.distance)
        
        print("")


    ax.set_xlabel("Time [s]", fontsize=14)
    ax.set_ylabel("Distance [km]", fontsize=14)
    ax.set_xlim(np.min(time_real) + timedelta(seconds=0), np.max(time_real) - timedelta(seconds=0))
    ax.set_ylim(0, 600)
    ax.margins(y=0)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.legend(loc="upper right", fontsize=14, handlelength=4, handleheight=2, markerscale=1.5)

    ## Format X-axis to show time in 'H:M:S' format
    formatter = mdates.DateFormatter('%H:%M:%S')
    ax.xaxis.set_major_formatter(formatter)
    
    figures.save(f"figures/event_detection/event_{event_index}.pdf")

    plt.show()

    return trimmed_time_starttime, trimmed_time_endtime, distance_all_trace, label1


def detection_method_for_all_events(ESEC_avalanches, stream, event_index, events_to_keep, label1, label2):
        """
        Perform detection method on seismic traces for one event.

        Parameters
        ----------
        stream : obspy.Stream
            The stream containing traces to be analyzed.
        ESEC_avalanches : pandas.DataFrame
            The ESEC.
        event_index : int
            Index of the event to be analyzed.
        events_to_keep : list
            The events to keep because the detection method can work.
        label1 : str
            The label. If true, the legend will appears.
        label2 : str
            The label. If true, the legend will appears.
        """

        ## Initialize lists to store detection results
        trimmed_time_starttime = []
        trimmed_time_endtime = []
        distance_all_trace = []
        detection_yes_or_no = []
        detection_stopped = False

        ## Loop over each trace in the stream
        for index, trace in enumerate(stream):
            print("Work on trace", index)

            if not detection_stopped:
                ## The detection method
                _, _, trimmed_time, _, _, _, _, _, detection_yes_or_no = detection_on_one_trace(trace, ESEC_avalanches, event_index, index, -30, 10, detection_yes_or_no)
            
                try:
                    ## Check if five consecutive detections are "FALSE". If true, the detection has stopped
                    if len(detection_yes_or_no) >= 5 and all(d == "FALSE" for d in detection_yes_or_no[-5:]):
                        print("Detection stopped: 5 consecutive FALSE detections !")
                        detection_stopped = True

                        ## Mark the trace as undetected
                        trimmed_time_starttime.append(np.nan)
                        trimmed_time_endtime.append(np.nan)
                        distance_all_trace.append(trace.stats.distance)

                    ## Plot the detected signal if detection has not stopped
                    if not detection_stopped:

                        ## Store results
                        trimmed_time_starttime.append(trimmed_time[0])
                        trimmed_time_endtime.append(trimmed_time[-1])
                        distance_all_trace.append(trace.stats.distance)

                ## If it is not possible to trace the detected seismic signal, it means that the method failed to detect a seismic signal so an error will be generated.
                except TypeError:
                    print("No detection !")
                    trimmed_time_starttime.append(np.nan)
                    trimmed_time_endtime.append(np.nan)
                    distance_all_trace.append(trace.stats.distance)

            ## If detection has stopped, plot only the raw waveform without detection
            else:
                trimmed_time_starttime.append(np.nan)
                trimmed_time_endtime.append(np.nan)
                distance_all_trace.append(trace.stats.distance)
            
            print("")


        ## Create DataFrame with the start and end time of the detected method and the distance of the stations
        df = pd.DataFrame({'start_time': trimmed_time_starttime, 'end_time': trimmed_time_endtime, 'distance': distance_all_trace})
        df = df.sort_values('distance') ## Sort by distance
        df = df.reset_index(drop=True) ## Reset index
        df['detection'] = df['start_time'].apply(lambda x: False if np.isnan(x) else True) ## In a new column named "detection", add True if detection is possible. Add False if detection is not possible
        df['duration'] = df['end_time'] - df['start_time'] ## Compute the duration of the event using the detected method
        
        ## Extract the volume of the event
        volume = [ESEC_avalanches["volume"][event_index]]

        ## Keep only event with a detection in the first trace
        if df['detection'].any():
            events_to_keep.append(event_index)

        ## Plot the results
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


def fit_line(frequencies, values):
    """
    Fits a linear model in the data.

    Parameters:
    ------------
    frequencies : numpy.ndarray
        The frequency values for the PSD.
    values : numpy.ndarray
        The PSD values.

    Returns:
    ---------
    model.coef_[0] : float
        The slope of the fitted line.
    model.intercept_ : float
        The intercept of the fitted line.
    """

    model = LinearRegression()     ## The linear model
    f = frequencies.reshape(-1, 1) ## Reshapes the frequencies array into a 2D array
    model.fit(f, values)           ## Adjust the model

    return model.coef_[0], model.intercept_


def find_split_frequency(frequencies, values, min_freq=2.0, max_freq=10.0):
    """
    Finds the split frequency by detecting discontinuities in a specific frequency range.
    
    Parameters:
    ------------
    frequencies : numpy.ndarray
        The array of frequencies.
    values : numpy.ndarray
        The array of PSD values corresponding to the frequencies.
    min_freq : float, optional
        The minimum frequency to search for the split.
    max_freq : float, optional
        The maximum frequency to search for the split.

    Returns:
    ---------
    float
        The frequency where a significant discontinuity is detected, which will be used for splitting low and high-frequency ranges.
    """
    
    ## Create a mask to filter frequencies within the range defined by min_freq and max_freq
    freq_mask = (frequencies >= min_freq) & (frequencies <= max_freq)
    filtered_frequencies = frequencies[freq_mask]
    filtered_values = values[freq_mask]

    ## Calculate the derivative of the PSD values
    derivative = np.diff(filtered_values) / np.diff(filtered_frequencies)

    ## Find the index of the split frequency
    split_index = np.argmax(np.abs(np.diff(derivative)))

    return filtered_frequencies[split_index + 1] if split_index < len(filtered_frequencies) - 1 else filtered_frequencies[-1]


def ajustement_de_segment(mask, frequencies_signal, psd_signal, ax, color='green', label="Ajustement bas", pltplot = True):
    """
    Fits a linear model to a segment of the PSD data and plots the result.

    Parameters:
    ------------
    mask : numpy.ndarray
        Boolean mask array to filter.
    frequencies_signal : numpy.ndarray
        The full array of frequency values.
    psd_signal : numpy.ndarray
        The full array of PSD values corresponding to the frequencies.
    ax : matplotlib axis
        The axis object on which to plot the adjusted line.
    color : str
        The color of the fitted line.
    label : str
        The label for the fitted line.
    pltplot : bool
        If True, uses plt to plot directly, otherwise uses the provided axis.

    Returns:
    ---------
    freq : numpy.ndarray
        The filtered values.
    slope : float
        The slope of the fitted line.
    intercept : float
        The intercept of the fitted line
    psd : numpy.ndarray
        The PSD values
    """

    ## Check if there are any valid points in the mask 
    if np.any(mask):

        ## Extract the frequencies and PSD values
        freq = frequencies_signal[mask]
        psd = psd_signal[mask]

        ## Fit the model
        slope, intercept = fit_line(np.log(freq), np.log(psd))

        ## Plot the result
        if pltplot == True:
            plt.loglog(freq, np.exp(slope * np.log(freq) + intercept), color=color, label=label)
        else:
            ax.loglog(freq, np.exp(slope * np.log(freq) + intercept), color=color, label=label)

    return freq, slope, intercept, psd


def plot_spectre(trace, ESEC_avalanches, trimmed_data, trace_index, event_index, conserv_result=False):
    """
    Plots the power spectral density (PSD) of a seismic signal with Welch method and fits two linear models.

    Parameters:
    ------------
    trace : obspy.Trace
        The seismic trace containing the signal data to be analyzed.
    ESEC_avalanches : pandas.DataFrame
        The ESEC.
    trimmed_data : numpy.ndarray
        The detected signal data.
    trace_index : int
        Index of the seismic trace.
    event_index : int
        Index of the avalanche event.
    conserv_result : bool, optional
        If True, saves the fitting parameters in a CSV file and the plots.
    """

    ## Initialize list to store results
    curve_params = []

    ## Welch parameters
    segment_duration = 20
    noverlap = 12
    nperseg = int(segment_duration * trace.stats.sampling_rate)

    ## Welch method
    frequencies_signal, psd_signal = welch(trimmed_data, window='hamming', fs=trace.stats.sampling_rate, nperseg=nperseg, noverlap=noverlap)

    ## Cut the spectrum between filter frequencies
    mask = (frequencies_signal > 0.5) & (frequencies_signal < 9)
    frequencies_signal = frequencies_signal[mask]
    psd_signal = psd_signal[mask]

    ## Find the split frequency
    split_freq = find_split_frequency(frequencies_signal, psd_signal, min_freq=1, max_freq=10)
    low_mask = frequencies_signal <= split_freq
    high_mask = frequencies_signal > split_freq

    ## Adjusting two models
    _, low_slope, low_intercept, low_psd = ajustement_de_segment(low_mask, frequencies_signal, psd_signal, plt, color='green', label="Modèle bas", pltplot = True)
    _, high_slope, high_intercept, high_psd = ajustement_de_segment(high_mask, frequencies_signal, psd_signal, plt, color='blue', label="Modèle haut", pltplot = True)

    ## Store results
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
        
    ## Save results in a dataframe
    if conserv_result == True:
        df = pd.DataFrame(curve_params)
        df.to_csv(f'features/1_fitting/data/curve_parameters_{event_index}.csv', index=False)
        
    plt.loglog(frequencies_signal, psd_signal, color="C1", label="Spectre du signal sismique détecté")
    plt.legend()
    plt.margins(x=0)
    plt.xscale("log")
    plt.xlabel('Fréquences (Hz)')
    plt.ylabel(r'Densité Spectrale de Puissance du déplacement ($\mathrm{\frac{m^{2}}{Hz}}$)')

    figures.save(f"features/1_fitting/pictures/fitting_on_trace_{trace_index}_in_event_{event_index}.pdf")

    plt.show()