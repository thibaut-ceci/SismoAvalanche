"""Manage creation of maps."""

import warnings

from cartopy import crs as ccrs
from cartopy import feature as cfeature
from matplotlib import pyplot as plt


def select_extent():
    # Pour sélectionner des zones pour le extent
    World = [-180, 180, -90, 90]
    Europe = [-5, 28, 37, 52]
    Japon = [123, 150, 4, 46]
    Groenland = [-75, -10, 55, 85]
    France = [-5, 15, 42, 51]
    Special_Amerique = [-170, -50, 10, 60]
    Iles_Mariana = [144, 147, 12.5, 20]
    Suisse = [5.8, 15, 45.25, 47.75]
    Chine = [73, 135, 18, 54]
    World_without_Antarctica = [-180, 180, -60, 90]
    manuel = [-170, -100, 30, 60]
    manuel2 = [5, 15, 43, 48]

    #Amerique_du_nord = [-170, -50, 10, 70] #Amérique du Nord
    #Europe2 = [-12, 40, 35, 70] #Europe
    #Groenland2 = [-56, -50, 70, 73] #Groenland
    #Special_Amerique2 = [-160, -100, 35, 60]
    #Canada = [-155, -132, 62, 58] #Canada
    #USA = [-127, -105, 36, 52.5] #USA

    return World, Europe, Japon, Groenland, France, Special_Amerique, Iles_Mariana, Suisse, Chine, World_without_Antarctica, manuel, manuel2


def add_equator_and_tropics(ax):
    #Ajout de l'équateur et des tropiques
    ax.plot([-180, 180], [0, 0], color="black", transform=ccrs.PlateCarree(), alpha=0.5, lw = 0.5)
    ax.plot([-180, 180], [23.4368, 23.4368], color="black", transform=ccrs.PlateCarree(), alpha=0.5, lw = 0.5)
    ax.plot([-180, 180], [-23.4368, -23.4368], color="black", transform=ccrs.PlateCarree(), alpha=0.5, lw = 0.5)

    #Add labels for the equator and tropics
    ax.text(-179, 1, 'Equateur', va='bottom', ha='left', color="black", alpha=0.5, fontsize=8, transform=ccrs.PlateCarree())
    ax.text(-179, 23.4368, 'Tropique du Cancer', va='bottom', ha='left', color="black", alpha=0.5, fontsize=8, transform=ccrs.PlateCarree())
    ax.text(-179, -23.4368, 'Tropique du Capricorne', va='bottom', ha='left', color="black", alpha=0.5, fontsize=8, transform=ccrs.PlateCarree())

    return ax
    

def coastlines_features(ax):
    #Add coastlines
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle="-", color="0.7")
    ax.add_feature(cfeature.OCEAN, color="lightblue")
    ax.add_feature(cfeature.LAND, color="0.95")
    ax.add_feature(cfeature.RIVERS, edgecolor="lightblue")
    ax.add_feature(cfeature.LAKES, edgecolor="lightblue")
    ax.add_feature(cfeature.COASTLINE, edgecolor="gray")
    ax.add_feature(cfeature.STATES, edgecolor="gray")
    ax.add_feature(cfeature.NaturalEarthFeature("cultural", "populated_places", "10m"))

    return ax


def gridlines(ax):
    # Add gridlines with labels
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.2, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    return ax


def show_select_area(avalanches, extent, add_equator_and_tropics_on_map = False, add_legend = True):

    warnings.filterwarnings("ignore", category=UserWarning)
        
    _ = plt.figure(figsize=(15,15))
    ax = plt.axes(projection=ccrs.Mercator())

    #Plot the top of the events
    ax.plot(avalanches.longitude, avalanches.latitude, "C3*", transform=ccrs.PlateCarree(), label="Evénement", markeredgecolor='black', markeredgewidth=0.50, markersize = 25)
    
    coastlines_features(ax)

    if add_equator_and_tropics_on_map == True:
        add_equator_and_tropics(ax)

    #Définition de l'étendue
    ax.set_extent(extent)

    #Add gridlines with labels
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.2, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    # gl.xlabel_style = {'size': 12, 'color': 'gray'}
    # gl.ylabel_style = {'size': 12, 'color': 'gray'}

    if add_legend == True:
        plt.legend()


def show(event, max_extent=1):
    """Show a map with the event and available stations,

    Parameters
    ----------
    event : obspy.core.event.Event
        The event to show on the map. Also contains the inventory of stations.
    """
    
    warnings.filterwarnings("ignore", category=UserWarning, message=".*Approximating coordinate system.*")
        
    # Create figure
    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.Mercator()}, figsize = (8, 12))#figsize = (10, 15), dpi=300)

    # Plot stations
    station_longitudes = []
    station_latitudes = []
    for net in event.inventory:
        for sta in net:
            station_longitudes.append(sta.longitude)
            station_latitudes.append(sta.latitude)
    ax.plot(station_longitudes, station_latitudes, "C0v", markersize=14, transform=ccrs.PlateCarree(), label="Stations", markerfacecolor='C0', markeredgecolor='black', markeredgewidth=0.40)
    
    # Plot event
    ax.plot(event.longitude,event.latitude,"C3*",markersize=16,transform=ccrs.PlateCarree(), label="Événement", markeredgecolor='black', markeredgewidth=0.40)

    ax.set_extent([event.longitude - max_extent, event.longitude + max_extent, event.latitude - max_extent/1.80, event.latitude + max_extent/1.80], crs=ccrs.PlateCarree())

    # Plot a circle with the radius
    # Filter out the warning about the projection issues
    ax.tissot(rad_km=555.55, lons=event.longitude, lats=event.latitude, n_samples=100, facecolor="C0", alpha=0.2)

    # Calculate the distance from each station to the event
    # distances = np.sqrt((np.array(station_longitudes) - event.longitude)**2 + (np.array(station_latitudes) - event.latitude)**2)

    # # Get the indices that would sort the distances
    # indices = np.argsort(distances)

    # # # Sort the longitudes, latitudes, and distances using the indices
    # station_longitudes = np.array(station_longitudes)[indices]
    # station_latitudes = np.array(station_latitudes)[indices]
    # distances = distances[indices]

    # # Plot the stations with numbers
    # for i, (lon, lat) in enumerate(zip(station_longitudes, station_latitudes)):
    #     ax.text(lon, lat, str(i), transform=ccrs.PlateCarree())
    
    # Plot the stations with numbers
    # for i, (lon, lat) in enumerate(zip(station_longitudes, station_latitudes)):
    #     color = 'orange' if i in [14, 95] else 'C0'
    #     ax.plot(lon, lat, marker='v', color=color, markersize=14, transform=ccrs.PlateCarree(), markeredgecolor='black', markeredgewidth=0.40)
        # ax.text(lon, lat, str(i), transform=ccrs.PlateCarree())

    # add_zebra_frame(ax, lw=5, crs=ccrs.Mercator(), zorder=3)
    # Labels
    plt.legend(loc="lower left", fontsize=12)

    # Title
    # title = title or event.title
    # plt.title(title, loc="left")

    # Labels
    # ax.coastlines()
    # ax.add_feature(cfeature.BORDERS, linestyle="-", color="0.7")
    ax.add_feature(cfeature.OCEAN, color="lightblue")
    # ax.add_feature(cfeature.LAND, color="0.95")
    ax.add_feature(cfeature.RIVERS, edgecolor="lightblue")
    ax.add_feature(cfeature.LAKES, edgecolor="lightblue")
    ax.add_feature(cfeature.COASTLINE, edgecolor="#9B9B9B", alpha=0.95)
    ax.add_feature(cfeature.STATES, edgecolor="#9B9B9B", alpha=0.3)
    # ax.add_feature(cfeature.NaturalEarthFeature("cultural", "populated_places", "10m"))
    
    gridlines = ax.gridlines(draw_labels=True, linestyle="--", color="#9B9B9B", alpha=0, linewidth=0.0001)

    # # Existing code
    gridlines.top_labels = False
    gridlines.right_labels = False

    # Add global location map at the top right corner
    ax_pos = ax.get_position()
    projection = ccrs.Orthographic(
        central_latitude=event.latitude, central_longitude=event.longitude
    )
    axins = fig.add_axes([ax_pos.x1 - 0.095, ax_pos.y1 - 0.085, 0.19, 0.19], projection=projection)
    axins.set_global()
    axins.coastlines(linewidth=0.5)
    axins.background_img()
    # axins.gridlines(draw_labels=False)
    axins.plot(event.longitude, event.latitude, "C3*", markersize=12, transform=ccrs.PlateCarree(), label="Événement", markerfacecolor="red", markeredgecolor='black', markeredgewidth=0.40)

    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)


def show_select_area_with_NoS(avalanches, extent, column):

    warnings.filterwarnings("ignore", category=UserWarning)

    avalanches_filtre_NOS0 = avalanches.loc[(avalanches[column] >= 3) & (avalanches[column] <= 10)]
    avalanches_filtre_NOS1 = avalanches.loc[(avalanches[column] >= 10) & (avalanches[column] <= 20)]
    avalanches_filtre_NOS2 = avalanches.loc[(avalanches[column] >= 20) & (avalanches[column] <= 50)]
    avalanches_filtre_NOS3 = avalanches.loc[(avalanches[column] >= 50) & (avalanches[column] <= 100)]

    fig = plt.figure(figsize=(10,10))
    ax = plt.axes(projection=ccrs.Mercator())

    ax.plot(avalanches_filtre_NOS3.longitude, avalanches_filtre_NOS3.latitude, "v", color="purple", transform=ccrs.PlateCarree(), label="50 - 100 stations", markeredgecolor='black', markeredgewidth=0.50, markersize = 20)    
    ax.plot(avalanches_filtre_NOS2.longitude, avalanches_filtre_NOS2.latitude, "v", color="red", transform=ccrs.PlateCarree(), label="20 - 50 stations", markeredgecolor='black', markeredgewidth=0.50, markersize = 16)
    ax.plot(avalanches_filtre_NOS1.longitude, avalanches_filtre_NOS1.latitude, "v", color="orange", transform=ccrs.PlateCarree(), label="10 - 20 stations", markeredgecolor='black', markeredgewidth=0.50, markersize = 12)
    ax.plot(avalanches_filtre_NOS0.longitude, avalanches_filtre_NOS0.latitude, "v", color="yellow", transform=ccrs.PlateCarree(), label="3 - 10 stations", markeredgecolor='black', markeredgewidth=0.50, markersize = 8)
    
    coastlines_features(ax)
    gridlines(ax)
    ax.set_extent(extent)
    plt.legend()


def show_select_area_with_volume(avalanches, extent, column):

    warnings.filterwarnings("ignore", category=UserWarning)

    avalanches_filtre_VOL0 = avalanches.loc[(avalanches[column] >= 0) & (avalanches[column] <= 100000)]
    avalanches_filtre_VOL3 = avalanches.loc[(avalanches[column] >= 100000) & (avalanches[column] <= 1000000)]
    avalanches_filtre_VOL4 = avalanches.loc[(avalanches[column] >= 1000000) & (avalanches[column] <= 10000000)]
    avalanches_filtre_VOL5 = avalanches.loc[(avalanches[column] >= 10000000)]

    _ = plt.figure(figsize=(10,10))
    ax = plt.axes(projection=ccrs.Mercator())

    ax.plot(avalanches_filtre_VOL5.longitude, avalanches_filtre_VOL5.latitude, "v", color="purple", transform=ccrs.PlateCarree(), label=r"10 000 000 - 72 000 000 $\mathrm{m}^3$", markeredgecolor='black', markeredgewidth=0.50, markersize = 20)
    ax.plot(avalanches_filtre_VOL4.longitude, avalanches_filtre_VOL4.latitude, "v", color="red",  transform=ccrs.PlateCarree(), label=r"1 000 000 - 10 000 000 $\mathrm{m}^3$", markeredgecolor='black', markeredgewidth=0.50, markersize = 16)
    ax.plot(avalanches_filtre_VOL3.longitude, avalanches_filtre_VOL3.latitude, "v", color="orange", transform=ccrs.PlateCarree(), label=r"100 000 - 1 000 000 $\mathrm{m}^3$", markeredgecolor='black', markeredgewidth=0.50, markersize = 12)
    ax.plot(avalanches_filtre_VOL0.longitude, avalanches_filtre_VOL0.latitude, "v", color="yellow", transform=ccrs.PlateCarree(), label=r"6 000 - 100 000 $\mathrm{m}^3$", markeredgecolor='black', markeredgewidth=0.50, markersize = 8)

    coastlines_features(ax)
    gridlines(ax)
    ax.set_extent(extent)
    plt.legend()


def show_select_area_with_subtype(avalanches, extent, column):

    warnings.filterwarnings("ignore", category=UserWarning)

    avalanches_filtre_subtype0 = avalanches.loc[avalanches[column] == 'Rock/ice/debris avalanches and slides']
    avalanches_filtre_subtype1 = avalanches.loc[avalanches[column] == 'Snow avalanches']

    fig = plt.figure(figsize=(10,10))
    ax = plt.axes(projection=ccrs.Mercator())

    ax.plot(avalanches_filtre_subtype0.longitude, avalanches_filtre_subtype0.latitude, "C1*", transform=ccrs.PlateCarree(), label="Rock/ice/debris avalanches and slides", markersize = 12, markeredgecolor='black', markeredgewidth=0.40)
    ax.plot(avalanches_filtre_subtype1.longitude, avalanches_filtre_subtype1.latitude, "C0*", transform=ccrs.PlateCarree(), label="Snow avalanches", markersize = 12, markeredgecolor='black', markeredgewidth=0.40)

    coastlines_features(ax)
    gridlines(ax)
    ax.set_extent(extent)
    plt.legend()


def show_select_area_with_type(avalanches, extent, column):

    warnings.filterwarnings("ignore", category=UserWarning)

    x = avalanches[column].value_counts()

    x = x.index.tolist()

    fig = plt.figure(figsize=(10, 10))
    ax = plt.axes(projection=ccrs.Mercator())

    colors = ["blue", "green", "red", "cyan", "magenta", "yellow", "black", "orange", "purple",
            "brown", "pink", "gray", "olive", "navy", "teal", "maroon", "gold", "indigo",
            "coral", "lime", "orchid", "peru"]

    for type, color in zip(x, colors):
        avalanches_filtre_subtype = avalanches.loc[avalanches[column] == type]
        ax.plot(avalanches_filtre_subtype.longitude, avalanches_filtre_subtype.latitude, "*", c=color, transform=ccrs.PlateCarree(), label=type, markersize = 15, markeredgecolor='black', markeredgewidth=0.40)

    coastlines_features(ax)
    gridlines(ax)
    ax.set_extent(extent)
    plt.legend(bbox_to_anchor=(1.0, 1), loc='upper left')


def show_select_area_with_length(avalanches, extent, column):

    warnings.filterwarnings("ignore", category=UserWarning)

    avalanches_filtre_L0 = avalanches.loc[(avalanches[column] >= 0) & (avalanches[column] <= 1000)]
    avalanches_filtre_L1 = avalanches.loc[(avalanches[column] >= 1000) & (avalanches[column] <= 5000)]
    avalanches_filtre_L2 = avalanches.loc[(avalanches[column] >= 5000) & (avalanches[column] <= 10000)]
    avalanches_filtre_L3 = avalanches.loc[(avalanches[column] >= 10000) & (avalanches[column] <= 20000)]
    avalanches_filtre_L4 = avalanches.loc[(avalanches[column] >= 20000)]

    fig = plt.figure(figsize=(10,10))
    ax = plt.axes(projection=ccrs.Mercator())


    ax.plot(avalanches_filtre_L0.longitude, avalanches_filtre_L0.latitude, "C0*", transform=ccrs.PlateCarree(), label="0 - 1000 m", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)
    ax.plot(avalanches_filtre_L1.longitude, avalanches_filtre_L1.latitude, "C1*", transform=ccrs.PlateCarree(), label="1000 - 5000 m", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)
    ax.plot(avalanches_filtre_L2.longitude, avalanches_filtre_L2.latitude, "C2*", transform=ccrs.PlateCarree(), label="5000 - 10000 m", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)
    ax.plot(avalanches_filtre_L3.longitude, avalanches_filtre_L3.latitude, "C3*", transform=ccrs.PlateCarree(), label="10000 - 20000 m", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)
    ax.plot(avalanches_filtre_L4.longitude, avalanches_filtre_L4.latitude, "C4*", transform=ccrs.PlateCarree(), label="+ de 20000 m", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)

    coastlines_features(ax)
    gridlines(ax)
    ax.set_extent(extent)
    plt.legend()


def show_select_area_with_height(avalanches, extent, column):

    warnings.filterwarnings("ignore", category=UserWarning)

    avalanches_filtre_H0 = avalanches.loc[(avalanches[column] >= 0) & (avalanches[column] <= 500)]
    avalanches_filtre_H1 = avalanches.loc[(avalanches[column] >= 500) & (avalanches[column] <= 1000)]
    avalanches_filtre_H2 = avalanches.loc[(avalanches[column] >= 1000) & (avalanches[column] <= 1500)]
    avalanches_filtre_H3 = avalanches.loc[(avalanches[column] >= 1500) & (avalanches[column] <= 2000)]
    avalanches_filtre_H4 = avalanches.loc[(avalanches[column] >= 2000)]

    fig = plt.figure(figsize=(10,10))
    ax = plt.axes(projection=ccrs.Mercator())

    ax.plot(avalanches_filtre_H0.longitude, avalanches_filtre_H0.latitude, "C0*", transform=ccrs.PlateCarree(), label="0 - 500 m", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)
    ax.plot(avalanches_filtre_H1.longitude, avalanches_filtre_H1.latitude, "C1*", transform=ccrs.PlateCarree(), label="500 - 1000 m", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)
    ax.plot(avalanches_filtre_H2.longitude, avalanches_filtre_H2.latitude, "C2*", transform=ccrs.PlateCarree(), label="1000 - 1500 m", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)
    ax.plot(avalanches_filtre_H3.longitude, avalanches_filtre_H3.latitude, "C3*", transform=ccrs.PlateCarree(), label="1500 - 2000 m", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)
    ax.plot(avalanches_filtre_H4.longitude, avalanches_filtre_H4.latitude, "C4*", transform=ccrs.PlateCarree(), label="+ de 2000 m", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)

    coastlines_features(ax)
    gridlines(ax)
    ax.set_extent(extent)
    plt.legend()


def show_select_area_with_date(avalanches, extent, column):

    warnings.filterwarnings("ignore", category=UserWarning)
    
    avalanches[column][0]

    avalanches_filtre_ARS0 = avalanches.loc[(avalanches[column] >= str(1900)) & (avalanches[column] <= str(1980))]
    avalanches_filtre_ARS1 = avalanches.loc[(avalanches[column] >= str(1980)) & (avalanches[column] <= str(1990))]
    avalanches_filtre_ARS2 = avalanches.loc[(avalanches[column] >= str(1990)) & (avalanches[column] <= str(2000))]
    avalanches_filtre_ARS3 = avalanches.loc[(avalanches[column] >= str(2000)) & (avalanches[column] <= str(2010))]
    avalanches_filtre_ARS4 = avalanches.loc[(avalanches[column] >= str(2010)) & (avalanches[column] <= str(2020))]
    avalanches_filtre_ARS5 = avalanches.loc[(avalanches[column] >= str(2020))]

    fig = plt.figure(figsize=(10,10))
    ax = plt.axes(projection=ccrs.Mercator())

    ax.plot(avalanches_filtre_ARS0.longitude, avalanches_filtre_ARS0.latitude, "C0*", transform=ccrs.PlateCarree(), label="avant 1980", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)
    ax.plot(avalanches_filtre_ARS1.longitude, avalanches_filtre_ARS1.latitude, "C1*", transform=ccrs.PlateCarree(), label="1980 - 1990", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)
    ax.plot(avalanches_filtre_ARS2.longitude, avalanches_filtre_ARS2.latitude, "C2*", transform=ccrs.PlateCarree(), label="1990 - 2000", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)
    ax.plot(avalanches_filtre_ARS3.longitude, avalanches_filtre_ARS3.latitude, "C3*", transform=ccrs.PlateCarree(), label="2000 - 2010", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)
    ax.plot(avalanches_filtre_ARS4.longitude, avalanches_filtre_ARS4.latitude, "C4*", transform=ccrs.PlateCarree(), label="2010 - 2020", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)
    ax.plot(avalanches_filtre_ARS5.longitude, avalanches_filtre_ARS5.latitude, "C5*", transform=ccrs.PlateCarree(), label="après 2020", markersize = 15, markeredgecolor='black', markeredgewidth=0.40)

    coastlines_features(ax)
    gridlines(ax)
    ax.set_extent(extent)
    plt.legend()