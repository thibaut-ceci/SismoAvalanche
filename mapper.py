"""
Manage creation of maps.

Write by Léonard Seydoux (seydoux@ipgp.fr) and Thibaut Céci (thi.ceci@gmail.com)
"""

import warnings

from cartopy import crs as ccrs
from cartopy import feature as cfeature
import matplotlib.pyplot as plt


def select_extent():
    """
    Returns the geographical extent (latitude and longitude) for the selected area.
    """
    World = [-180, 180, -90, 90]
    Europe = [-5, 28, 37, 52]
    Japon = [123, 150, 4, 46]
    Groenland = [-75, -10, 55, 85]
    France = [-5, 15, 42, 51]
    Amerique = [-170, -50, 10, 60]
    Amerique_du_nord = [-170, -50, 10, 70]
    Iles_Mariana = [144, 147, 12.5, 20]
    Suisse = [5.8, 15, 45.25, 47.75]
    Chine = [73, 135, 18, 54]
    World_without_Antarctica = [-180, 180, -60, 90]
    Canada = [-155, -132, 62, 58]
    USA = [-127, -105, 36, 52.5]

    return World, Europe, Japon, Groenland, France, Amerique, Amerique_du_nord, Iles_Mariana, Suisse, Chine, World_without_Antarctica, Canada, USA


def add_equator_and_tropics(ax):
    """
    Adds the Equator and Tropics (Cancer and Capricorn) lines to a map.

    Parameters:
    ------------
    ax : matplotlib.axes._subplots.AxesSubplot
        The axis to plot on.

    Returns:
    --------
    ax : matplotlib.axes._subplots.AxesSubplot
        The axis with the equator and tropics added.
    """
    ## Add Equator and Tropics
    ax.plot([-180, 180], [0, 0], color="black", transform=ccrs.PlateCarree(), alpha=0.5, lw = 0.5)
    ax.plot([-180, 180], [23.4368, 23.4368], color="black", transform=ccrs.PlateCarree(), alpha=0.5, lw = 0.5)
    ax.plot([-180, 180], [-23.4368, -23.4368], color="black", transform=ccrs.PlateCarree(), alpha=0.5, lw = 0.5)

    ## Add labels for Equator and Tropics
    ax.text(-179, 1, 'Equator', va='bottom', ha='left', color="black", alpha=0.5, fontsize=8, transform=ccrs.PlateCarree())
    ax.text(-179, 23.4368, 'Tropic of Cancer', va='bottom', ha='left', color="black", alpha=0.5, fontsize=8, transform=ccrs.PlateCarree())
    ax.text(-179, -23.4368, 'Tropic of Capricorn', va='bottom', ha='left', color="black", alpha=0.5, fontsize=8, transform=ccrs.PlateCarree())

    return ax
    

def coastlines_features(ax):
    """
    Adds coastlines and geographical features (land, ocean, rivers, lakes, borders) to a map.

    Parameters:
    ------------
    ax : matplotlib.axes._subplots.AxesSubplot
        The axis to plot on.

    Returns:
    --------
    ax : matplotlib.axes._subplots.AxesSubplot
        The axis with features added.
    """

    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle="-", color="0.7")
    ax.add_feature(cfeature.OCEAN, color="lightblue")
    ax.add_feature(cfeature.LAND, color="0.95")
    ax.add_feature(cfeature.RIVERS, edgecolor="lightblue")
    ax.add_feature(cfeature.LAKES, edgecolor="lightblue")
    ax.add_feature(cfeature.COASTLINE, edgecolor="#9B9B9B", alpha=0.95)
    ax.add_feature(cfeature.STATES, edgecolor="#9B9B9B", alpha=0.3)
    ax.add_feature(cfeature.NaturalEarthFeature("cultural", "populated_places", "10m"))

    return ax


def gridlines(ax):
    """
    Adds gridlines in a map.

    Parameters:
    ------------
    ax : matplotlib.axes._subplots.AxesSubplot
        The axis to add gridlines on.

    Returns:
    --------
    ax : matplotlib.axes._subplots.AxesSubplot
        The axis with gridlines added.
    """

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.2, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    return ax


def show_select_area(avalanches, extent, add_equator_and_tropics_on_map = False, add_legend = True):
    """
    Plot a map with avalanche events

    Parameters:
    ------------
    avalanches : pandas.DataFrame
        The ESEC.
    extent : list
        Extent of the map.
    add_equator_and_tropics_on_map : bool
        If True, adds the equator and tropics lines to the map.
    add_legend : bool
        If True, adds a legend to the map.
    """

    ## Ignore specific warnings
    warnings.filterwarnings("ignore", category=UserWarning)
        
    _ = plt.figure(figsize=(15,15))
    ax = plt.axes(projection=ccrs.Mercator())

    ax.plot(avalanches.longitude, avalanches.latitude, "C3*", transform=ccrs.PlateCarree(), label="Avalanche", markeredgecolor='black', markeredgewidth=0.50, markersize = 25)
    
    ## Add features
    coastlines_features(ax)

    ## Optionally add the equator and tropics
    if add_equator_and_tropics_on_map == True:
        add_equator_and_tropics(ax)

    ## Set the extent of the map
    ax.set_extent(extent)

    ## Add gridlines with labels
    gridlines(ax)

    ## Optionally add a legend
    if add_legend == True:
        plt.legend()


def show(event, max_extent = 10):
    """
    Show a map with the event and available stations

    Parameters
    ----------
    event : obspy.core.event.Event
        The event to show on the map. Also contains the inventory of stations.
    max_extent : float
        The maximum extent (in degrees) around the avalanche.
    """
    
    ## Ignore specific warnings
    warnings.filterwarnings("ignore", category=UserWarning, message=".*Approximating coordinate system.*")

    ## Plot stations
    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.Mercator()}, figsize = (8, 12))

    station_longitudes = []
    station_latitudes = []
    for net in event.inventory:
        for sta in net:
            station_longitudes.append(sta.longitude)
            station_latitudes.append(sta.latitude)
    ax.plot(station_longitudes, station_latitudes, "C0v", markersize=14, transform=ccrs.PlateCarree(), label="Stations", markerfacecolor='C0', markeredgecolor='black', markeredgewidth=0.40)
    
    ## Plot the event
    ax.plot(event.longitude,event.latitude,"C3*",markersize=16,transform=ccrs.PlateCarree(), label="Avalanche",  markerfacecolor="red", markeredgecolor='black', markeredgewidth=0.40)

    ## Set the extent of the map
    ax.set_extent([event.longitude - max_extent, event.longitude + max_extent, event.latitude - max_extent/1.80, event.latitude + max_extent/1.80], crs=ccrs.PlateCarree())

    ## Plot a circle
    ax.tissot(rad_km=555.55, lons=event.longitude, lats=event.latitude, n_samples=100, facecolor="C0", alpha=0.2)

    ## Labels and features
    plt.legend(loc="lower left", fontsize=12)
    coastlines_features(ax)
    gridlines(ax)

    ## Add global location map at the top right corner
    ax_pos = ax.get_position()
    projection = ccrs.Orthographic(central_latitude=event.latitude, central_longitude=event.longitude)
    axins = fig.add_axes([ax_pos.x1 - 0.095, ax_pos.y1 - 0.085, 0.19, 0.19], projection=projection)
    axins.set_global()
    axins.coastlines(linewidth=0.5)
    axins.background_img()

    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)


def show_select_area_with_NoS(avalanches, extent, column):
    """
    Show a map with selected areas based on the number of stations.

    Parameters
    ----------
    avalanches : pandas.DataFrame
        The ESEC.
    extent : list
        The extent to display on the map.
    column : str
        The column "number of stations" in the ESEC.
    """
    
    ## Ignore specific warnings
    warnings.filterwarnings("ignore", category=UserWarning)

    ## Split data for plotting
    avalanches_filtre_NOS0 = avalanches.loc[(avalanches[column] >= 3) & (avalanches[column] <= 10)]
    avalanches_filtre_NOS1 = avalanches.loc[(avalanches[column] >= 10) & (avalanches[column] <= 20)]
    avalanches_filtre_NOS2 = avalanches.loc[(avalanches[column] >= 20) & (avalanches[column] <= 50)]
    avalanches_filtre_NOS3 = avalanches.loc[(avalanches[column] >= 50) & (avalanches[column] <= 100)]

    ## Plot
    fig = plt.figure(figsize=(10,10))
    ax = plt.axes(projection=ccrs.Mercator())

    ax.plot(avalanches_filtre_NOS3.longitude, avalanches_filtre_NOS3.latitude, "v", color="purple", transform=ccrs.PlateCarree(), label="50 - 100 stations", markeredgecolor='black', markeredgewidth=0.50, markersize = 20)    
    ax.plot(avalanches_filtre_NOS2.longitude, avalanches_filtre_NOS2.latitude, "v", color="red", transform=ccrs.PlateCarree(), label="20 - 50 stations", markeredgecolor='black', markeredgewidth=0.50, markersize = 16)
    ax.plot(avalanches_filtre_NOS1.longitude, avalanches_filtre_NOS1.latitude, "v", color="orange", transform=ccrs.PlateCarree(), label="10 - 20 stations", markeredgecolor='black', markeredgewidth=0.50, markersize = 12)
    ax.plot(avalanches_filtre_NOS0.longitude, avalanches_filtre_NOS0.latitude, "v", color="yellow", transform=ccrs.PlateCarree(), label="3 - 10 stations", markeredgecolor='black', markeredgewidth=0.50, markersize = 8)
    
    ## Add coastlines, gridlines, set extent and show legend
    coastlines_features(ax)
    gridlines(ax)
    ax.set_extent(extent)
    plt.legend()


def show_select_area_with_volume(avalanches, extent, column):
    """
    Show a map with selected areas based on the volume of the avalanches.

    Parameters
    ----------
    avalanches : pandas.DataFrame
        The ESEC.
    extent : list
        The extent to display on the map.
    column : str
        The column "volume" in the ESEC.
    """

    ## Ignore specific warnings
    warnings.filterwarnings("ignore", category=UserWarning)

    ## Split data for plotting
    avalanches_filtre_VOL0 = avalanches.loc[(avalanches[column] >= 0) & (avalanches[column] <= 100000)]
    avalanches_filtre_VOL3 = avalanches.loc[(avalanches[column] >= 100000) & (avalanches[column] <= 1000000)]
    avalanches_filtre_VOL4 = avalanches.loc[(avalanches[column] >= 1000000) & (avalanches[column] <= 10000000)]
    avalanches_filtre_VOL5 = avalanches.loc[(avalanches[column] >= 10000000)]

    ## Plot
    _ = plt.figure(figsize=(10,10))
    ax = plt.axes(projection=ccrs.Mercator())

    ax.plot(avalanches_filtre_VOL5.longitude, avalanches_filtre_VOL5.latitude, "v", color="purple", transform=ccrs.PlateCarree(), label=r"10 000 000 - 72 000 000 $\mathrm{m}^3$", markeredgecolor='black', markeredgewidth=0.50, markersize = 20)
    ax.plot(avalanches_filtre_VOL4.longitude, avalanches_filtre_VOL4.latitude, "v", color="red",  transform=ccrs.PlateCarree(), label=r"1 000 000 - 10 000 000 $\mathrm{m}^3$", markeredgecolor='black', markeredgewidth=0.50, markersize = 16)
    ax.plot(avalanches_filtre_VOL3.longitude, avalanches_filtre_VOL3.latitude, "v", color="orange", transform=ccrs.PlateCarree(), label=r"100 000 - 1 000 000 $\mathrm{m}^3$", markeredgecolor='black', markeredgewidth=0.50, markersize = 12)
    ax.plot(avalanches_filtre_VOL0.longitude, avalanches_filtre_VOL0.latitude, "v", color="yellow", transform=ccrs.PlateCarree(), label=r"6 000 - 100 000 $\mathrm{m}^3$", markeredgecolor='black', markeredgewidth=0.50, markersize = 8)

    ## Add coastlines, gridlines, set extent and show legend
    coastlines_features(ax)
    gridlines(ax)
    ax.set_extent(extent)
    plt.legend()


def show_select_area_with_subtype(avalanches, extent, column):
    """
    Show a map with selected areas based on the subtype.

    Parameters
    ----------
    avalanches : pandas.DataFrame
        The ESEC.
    extent : list
        The extent to display on the map.
    column : str
        The column "subtype" in the ESEC.
    """
    
    ## Ignore specific warnings
    warnings.filterwarnings("ignore", category=UserWarning)

    ## Split data for plotting
    avalanches_filtre_subtype0 = avalanches.loc[avalanches[column] == 'Rock/ice/debris avalanches and slides']
    avalanches_filtre_subtype1 = avalanches.loc[avalanches[column] == 'Snow avalanches']

    ## Plot
    fig = plt.figure(figsize=(10,10))
    ax = plt.axes(projection=ccrs.Mercator())

    ax.plot(avalanches_filtre_subtype0.longitude, avalanches_filtre_subtype0.latitude, "C1*", transform=ccrs.PlateCarree(), label="Rock/ice/debris avalanches and slides", markersize = 12, markeredgecolor='black', markeredgewidth=0.40)
    ax.plot(avalanches_filtre_subtype1.longitude, avalanches_filtre_subtype1.latitude, "C0*", transform=ccrs.PlateCarree(), label="Snow avalanches", markersize = 12, markeredgecolor='black', markeredgewidth=0.40)

    ## Add coastlines, gridlines, set extent and show legend
    coastlines_features(ax)
    gridlines(ax)
    ax.set_extent(extent)
    plt.legend()


def show_select_area_with_type(avalanches, extent, column):
    """
    Show a map with selected areas based on the type.

    Parameters
    ----------
    avalanches : pandas.DataFrame
        The ESEC.
    extent : list
        The extent to display on the map.
    column : str
        The column "type" in the ESEC.
    """

    ## Ignore specific warnings
    warnings.filterwarnings("ignore", category=UserWarning)

    ## Split data for plotting
    value = avalanches[column].value_counts()
    value = value.index.tolist()

    ## Plot
    fig = plt.figure(figsize=(10, 10))
    ax = plt.axes(projection=ccrs.Mercator())

    colors = ["blue", "green", "red", "cyan", "magenta", "yellow", "black", "orange", "purple", "brown", "pink", "gray", "olive", "navy", "teal", "maroon", "gold", "indigo", "coral", "lime", "orchid", "peru"]

    for type, color in zip(value, colors):
        avalanches_filtre_subtype = avalanches.loc[avalanches[column] == type]
        ax.plot(avalanches_filtre_subtype.longitude, avalanches_filtre_subtype.latitude, "*", c=color, transform=ccrs.PlateCarree(), label=type, markersize = 15, markeredgecolor='black', markeredgewidth=0.40)

    ## Add coastlines, gridlines, set extent and show legend
    coastlines_features(ax)
    gridlines(ax)
    ax.set_extent(extent)
    plt.legend(bbox_to_anchor=(1.0, 1), loc='upper left')


def show_select_area_with_length(avalanches, extent, column):
    """
    Show a map with selected areas based on the length of the avalanches.

    Parameters
    ----------
    avalanches : pandas.DataFrame
        The ESEC.
    extent : list
        The extent to display on the map.
    column : str
        The column "length" in the ESEC.
    """

    ## Ignore specific warnings
    warnings.filterwarnings("ignore", category=UserWarning)

    ## Split data for plotting
    avalanches_filtre_L0 = avalanches.loc[(avalanches[column] >= 0) & (avalanches[column] <= 1000)]
    avalanches_filtre_L1 = avalanches.loc[(avalanches[column] >= 1000) & (avalanches[column] <= 5000)]
    avalanches_filtre_L2 = avalanches.loc[(avalanches[column] >= 5000) & (avalanches[column] <= 10000)]
    avalanches_filtre_L3 = avalanches.loc[(avalanches[column] >= 10000) & (avalanches[column] <= 20000)]
    avalanches_filtre_L4 = avalanches.loc[(avalanches[column] >= 20000)]

    ## Plot
    fig = plt.figure(figsize=(10,10))
    ax = plt.axes(projection=ccrs.Mercator())

    ax.plot(avalanches_filtre_L0.longitude, avalanches_filtre_L0.latitude, "C0*", transform=ccrs.PlateCarree(), label="0 - 1000 m", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)
    ax.plot(avalanches_filtre_L1.longitude, avalanches_filtre_L1.latitude, "C1*", transform=ccrs.PlateCarree(), label="1000 - 5000 m", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)
    ax.plot(avalanches_filtre_L2.longitude, avalanches_filtre_L2.latitude, "C2*", transform=ccrs.PlateCarree(), label="5000 - 10000 m", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)
    ax.plot(avalanches_filtre_L3.longitude, avalanches_filtre_L3.latitude, "C3*", transform=ccrs.PlateCarree(), label="10000 - 20000 m", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)
    ax.plot(avalanches_filtre_L4.longitude, avalanches_filtre_L4.latitude, "C4*", transform=ccrs.PlateCarree(), label="+ de 20000 m", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)

    ## Add coastlines, gridlines, set extent and show legend
    coastlines_features(ax)
    gridlines(ax)
    ax.set_extent(extent)
    plt.legend()


def show_select_area_with_height(avalanches, extent, column):
    """
    Show a map with selected areas based on the height of the avalanches.

    Parameters
    ----------
    avalanches : pandas.DataFrame
        The ESEC.
    extent : list
        The extent to display on the map.
    column : str
        The column "height" in the ESEC.
    """
   
    ## Ignore specific warnings
    warnings.filterwarnings("ignore", category=UserWarning)

    ## Split data for plotting
    avalanches_filtre_H0 = avalanches.loc[(avalanches[column] >= 0) & (avalanches[column] <= 500)]
    avalanches_filtre_H1 = avalanches.loc[(avalanches[column] >= 500) & (avalanches[column] <= 1000)]
    avalanches_filtre_H2 = avalanches.loc[(avalanches[column] >= 1000) & (avalanches[column] <= 1500)]
    avalanches_filtre_H3 = avalanches.loc[(avalanches[column] >= 1500) & (avalanches[column] <= 2000)]
    avalanches_filtre_H4 = avalanches.loc[(avalanches[column] >= 2000)]

    ## Plot
    fig = plt.figure(figsize=(10,10))
    ax = plt.axes(projection=ccrs.Mercator())

    ax.plot(avalanches_filtre_H0.longitude, avalanches_filtre_H0.latitude, "C0*", transform=ccrs.PlateCarree(), label="0 - 500 m", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)
    ax.plot(avalanches_filtre_H1.longitude, avalanches_filtre_H1.latitude, "C1*", transform=ccrs.PlateCarree(), label="500 - 1000 m", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)
    ax.plot(avalanches_filtre_H2.longitude, avalanches_filtre_H2.latitude, "C2*", transform=ccrs.PlateCarree(), label="1000 - 1500 m", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)
    ax.plot(avalanches_filtre_H3.longitude, avalanches_filtre_H3.latitude, "C3*", transform=ccrs.PlateCarree(), label="1500 - 2000 m", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)
    ax.plot(avalanches_filtre_H4.longitude, avalanches_filtre_H4.latitude, "C4*", transform=ccrs.PlateCarree(), label="+ de 2000 m", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)

    ## Add coastlines, gridlines, set extent and show legend
    coastlines_features(ax)
    gridlines(ax)
    ax.set_extent(extent)
    plt.legend()


def show_select_area_with_date(avalanches, extent, column):
    """
    Show a map with selected areas based on the starttime of the avalanches.

    Parameters
    ----------
    avalanches : pandas.DataFrame
        The ESEC.
    extent : list
        The extent to display on the map.
    column : str
        The column "starttime" in the ESEC.
    """

    ## Ignore specific warnings
    warnings.filterwarnings("ignore", category=UserWarning)
    
    ## Split data for plotting
    avalanches[column][0]

    avalanches_filtre_ARS0 = avalanches.loc[(avalanches[column] >= str(1900)) & (avalanches[column] <= str(1980))]
    avalanches_filtre_ARS1 = avalanches.loc[(avalanches[column] >= str(1980)) & (avalanches[column] <= str(1990))]
    avalanches_filtre_ARS2 = avalanches.loc[(avalanches[column] >= str(1990)) & (avalanches[column] <= str(2000))]
    avalanches_filtre_ARS3 = avalanches.loc[(avalanches[column] >= str(2000)) & (avalanches[column] <= str(2010))]
    avalanches_filtre_ARS4 = avalanches.loc[(avalanches[column] >= str(2010)) & (avalanches[column] <= str(2020))]
    avalanches_filtre_ARS5 = avalanches.loc[(avalanches[column] >= str(2020))]

    ## Plot
    fig = plt.figure(figsize=(10,10))
    ax = plt.axes(projection=ccrs.Mercator())

    ax.plot(avalanches_filtre_ARS0.longitude, avalanches_filtre_ARS0.latitude, "C0*", transform=ccrs.PlateCarree(), label="avant 1980", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)
    ax.plot(avalanches_filtre_ARS1.longitude, avalanches_filtre_ARS1.latitude, "C1*", transform=ccrs.PlateCarree(), label="1980 - 1990", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)
    ax.plot(avalanches_filtre_ARS2.longitude, avalanches_filtre_ARS2.latitude, "C2*", transform=ccrs.PlateCarree(), label="1990 - 2000", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)
    ax.plot(avalanches_filtre_ARS3.longitude, avalanches_filtre_ARS3.latitude, "C3*", transform=ccrs.PlateCarree(), label="2000 - 2010", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)
    ax.plot(avalanches_filtre_ARS4.longitude, avalanches_filtre_ARS4.latitude, "C4*", transform=ccrs.PlateCarree(), label="2010 - 2020", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)
    ax.plot(avalanches_filtre_ARS5.longitude, avalanches_filtre_ARS5.latitude, "C5*", transform=ccrs.PlateCarree(), label="après 2020", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)

    ## Add coastlines, gridlines, set extent and show legend
    coastlines_features(ax)
    gridlines(ax)
    ax.set_extent(extent)
    plt.legend()