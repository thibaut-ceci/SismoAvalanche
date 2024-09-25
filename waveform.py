"""
Manage seismic waveforms with ObsPy.

Write by Léonard Seydoux (seydoux@ipgp.fr) and Thibaut Céci (thi.ceci@gmail.com)
"""

from obspy import UTCDateTime, Stream
from obspy.clients.fdsn import Client
from obspy.geodetics import locations2degrees
from tqdm import tqdm

## Connect to the FSDN server, or the IRIS datacenter. Because this variable is defined outside any function, it is visible by all functions of this module.
client = Client("IRIS")


def download_inventory(event, maxradius=1.0, retries=10):
    """
    Downloads the station inventory for a given seismic event.
    
    Parameters:
    -----------
    event : Event
        Seismic event object containing attributes like latitude, longitude, starttime, and endtime.
    maxradius : float
        Maximum radius around the event's location to search for stations (in deg).
    retries : int
        Number of retry attempts in case of failure.
        
    Returns:
    --------
    Inventory
        Inventory object if successful
    """
    
    ## Start time and end time of the event
    start, end = UTCDateTime(event.starttime), UTCDateTime(event.endtime)

    ## Download station inventory from client
    for attempt in range(retries):
        try:
            return client.get_stations(
                latitude=event.latitude,
                longitude=event.longitude,
                startbefore=start,
                endafter=end,
                maxradius=maxradius,
                matchtimeseries=True,
                channel="BH*,HH*",
            )
        
        ## If an error occurs, the code restarts
        except Exception as e:
            print(f"Error for download inventory. Attempt {attempt + 1} of {retries}. Error: {e}")
            #time.sleep(delay)
    return []


def download_stream(event, time_margins=150, print_error=False, retries=3):
    """
    Download stream of waveforms for a given event.

    Parameters
    ----------
    event : pd.Series
        A series containing the event information.
    time_margins : float
        The number of seconds to add to the start and end times of the event.
    print_error : bool
        Whether to print error messages.
    retries : int
        Number of times to retry the request in case of an error.

    Returns
    -------
    stream : obspy.Stream
        A stream of waveforms for the event.
    """

    ## Start time and end time of the event
    start, end = UTCDateTime(event.starttime), UTCDateTime(event.endtime)
    start, end = start - time_margins, end + time_margins

    ## The stream for adding traces
    stream = Stream()

    ## Download and add traces in stream
    for network in tqdm(event.inventory):
        if any(char.isdigit() for char in network.code):
            continue

        for station in network:
            for attempt in range(retries):
                try:
                    traces = client.get_waveforms(
                        network.code,
                        station.code,
                        "*",
                        "BH*,HH*",
                        start,
                        end,
                        attach_response=True,
                    )

                    traces.merge(method=1, fill_value="interpolate")

                    distance = locations2degrees(
                        event.latitude,
                        event.longitude,
                        station.latitude,
                        station.longitude,
                    )

                    for trace in traces:
                        trace.stats.distance = distance * 111.19

                    stream += traces
                    break  ## Break the retry loop if successful
                
                ## If an error occurs, the code restarts
                except Exception as e:
                    if print_error:
                        print(f"Error with station {station}. Attempt {attempt + 1} of {retries}. Error: {e}")
                    continue  ## Retry

    return stream