"""Manage seismic waveforms with ObsPy."""

from obspy import UTCDateTime, Stream
from obspy.clients.fdsn import Client
from obspy.geodetics import locations2degrees
from tqdm import tqdm
import time

# Connect to the FSDN server, or the IRIS datacenter. Because this variable is
# defined outside any function, it is visible by all functions of this module.
client = Client("IRIS")


def download_inventory(event, maxradius=1.0, retries=10, delay=0):
    start, end = UTCDateTime(event.starttime), UTCDateTime(event.endtime)

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
        except Exception as e:
            print(f"ERROR FOR DOWNLOAD INVENTORY. Attempt {attempt + 1} of {retries}. Error: {e}")
            time.sleep(delay)
    return []




def download_stream(event, time_margins=150, print_error=False, retries=3, delay=5):
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
    delay : int
        Delay between retries in seconds.

    Returns
    -------
    stream : obspy.Stream
        A stream of waveforms for the event.
    """
    start, end = UTCDateTime(event.starttime), UTCDateTime(event.endtime)
    start, end = start - time_margins, end + time_margins

    stream = Stream()

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
                    break  # Break the retry loop if successful

                except Exception as e:
                    if print_error:
                        print(f"Error with station {station}. Attempt {attempt + 1} of {retries}. Error: {e}")
                    time.sleep(delay)
                    continue  # Retry

    return stream