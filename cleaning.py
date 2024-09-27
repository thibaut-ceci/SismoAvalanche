"""
Cleaning data in ESEC

This module contains functions to clean data from ESEC.
"""

import os
from tqdm.notebook import tqdm

import obspy
from obspy import ObsPyException

tqdm.pandas()


def events_without_digital_data(ESEC_avalanches):
    """
    Remove the events with no digital data.

    Parameters:
    ------------
    ESEC_avalanches : pandas.DataFrame
        The ESEC.

    Returns:
    ---------
    ESEC_avalanches : pandas.DataFrame
        The new ESEC without the events with no digital data.
    """

    ## The number of events before filtering
    initial_length = len(ESEC_avalanches)

    ## Remove events associated with "no digital data"
    ESEC_avalanches = ESEC_avalanches[ESEC_avalanches["datalocation"] != "No digital data"]

    ## The number of events after filtering
    final_length = len(ESEC_avalanches)

    ## Number of events removed
    removed_events = initial_length - final_length

    ## Print the number of events removed
    print(f"Number of events removed : {removed_events}")

    return ESEC_avalanches


def events_without_pickle_file(ESEC_avalanches):
    """
    Function to remove the events without pickle file.

    Parameters:
    ------------
    ESEC_avalanches : pandas.DataFrame
        The ESEC.

    Returns:
    ---------
    ESEC_avalanches : pandas.DataFrame
        The new ESEC without the events without pickle file.
    """

    ## Loop over all events
    for numero_event in tqdm(ESEC_avalanches["numero"], total=len(ESEC_avalanches)):

        ## Check if events have pickle file
        try:
            stream = obspy.read(f"sismogrammes/{numero_event:03d}.pickle")

        ## If an event has not pickle file, it is removed
        except FileNotFoundError:
            print("No pickle file in stream " + str(numero_event))

            ESEC_avalanches = ESEC_avalanches.drop(numero_event)

            print("line " + str(numero_event) + " removed")
            print("")

            continue

    print("In this catalog, there are now", len(ESEC_avalanches), "avalanches.")

    return ESEC_avalanches


def trace_with_a_frequency_below_of(a_frequency_below_of, ESEC_avalanches):
    """
    Function to remove traces with a frequency below of 20 Hz (for example).

    Parameters:
    ------------
    a_frequency_below_of : int
        The frequency to remove the trace
    ESEC_avalanches : pandas.DataFrame
        The ESEC.

    Returns:
    ---------
    ESEC_avalanches : pandas.DataFrame
        The ESEC without removed events.
    """

    ## List to store removed events
    events_to_remove = []

    ## Check if all traces in an event have a frequency below of 20 Hz
    for numero_event in tqdm(ESEC_avalanches["numero"], total=len(ESEC_avalanches)):

        ## Load the event
        stream = obspy.read(f"sismogrammes/{numero_event:03d}.pickle", format="PICKLE")
        print("Select stream " + str(numero_event))

        ## Store the trace index if it has a frequency below of 20 Hz
        traces_to_remove = []

        for numero_trace, trace in enumerate(stream):
            fs = trace.stats.sampling_rate

            if fs <= a_frequency_below_of:
                print("Trace " + str(numero_trace) + " has freq = " + str(fs))
                traces_to_remove.append(numero_trace)

        print(f"Number of traces with a frequency below of {a_frequency_below_of} Hz : " + str(len(traces_to_remove)))

        ## Remove traces marked for removal (in "traces_to_remove")
        for index in sorted(traces_to_remove, reverse=True):
            stream.pop(index)
            print("Trace " + str(index) + " deleted")
        
        ## Write the new stream without the traces with a frequency below of 20 Hz.
        try:
            stream.write(f"sismogrammes/{numero_event:03d}.pickle", format="PICKLE")

        ## An error occurs if the stream is empty. So, the event is removed
        except ObsPyException:
            print("All traces were deleted in stream", str(numero_event))
            events_to_remove.append(numero_event)
            os.remove(f"sismogrammes/{numero_event:03d}.pickle")
            continue
        
        print("Stream", str(numero_event), "has", len(stream), "traces")
        print("End of stream", str(numero_event))
        print("")

    ## Remove events from ESEC after the loop
    ESEC_avalanches = ESEC_avalanches[~ESEC_avalanches["numero"].isin(events_to_remove)]

    return ESEC_avalanches
    

def trace_without_metadata(ESEC_avalanches):
    """
    Function to remove the traces without instrumental response

    Parameters:
    ------------
    ESEC_avalanches : pandas.DataFrame
        The ESEC.
    """

    ## Check if all the traces in an event have an instrumental response
    for numero_event in tqdm(ESEC_avalanches["numero"], total=len(ESEC_avalanches)):

        ## Load the event
        stream = obspy.read(f"sismogrammes/{numero_event:03d}.pickle", format="PICKLE")
        print("Select stream " + str(numero_event))

        ## Check if the trace has an instrumental response and store the trace index if not
        traces_to_remove = []

        for numero_trace, trace in enumerate(stream):
            try:
                trace.remove_response()
            except TypeError:
                print("Trace " + str(numero_trace) + " has no instrumental response")
                traces_to_remove.append(numero_trace)

        print("Number of traces having no instrumental response : " + str(len(traces_to_remove)))

        ## Reload the stream because all the traces are in displacement now.
        stream = obspy.read(f"sismogrammes/{numero_event:03d}.pickle", format="PICKLE")

        ## Remove traces marked for removal
        for index in sorted(traces_to_remove, reverse=True):
            stream.pop(index)
            
            print("Trace " + str(index) + " deleted")

        ## Save the new stream without the traces without instrumental response
        stream.write(f"sismogrammes/{numero_event:03d}.pickle", format="PICKLE")
        print("End of stream " + str(numero_event))
        print("")
