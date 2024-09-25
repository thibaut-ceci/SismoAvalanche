import covseisnet as csn
import matplotlib.pyplot as plt
import obspy
import figures
import aside.analysis2 as analysis2
import computations as cp

import numpy as np

from tqdm.notebook import tqdm

import aside.analysis2 as analysis2

tqdm.pandas()


def filter_stream_with_covseisnet(event_index, trim):
    stream = obspy.read(f"sismogrammes/{event_index:03d}.pickle")
    stream = stream.sort(keys=["distance"])
    stream = stream.select(component="Z")

    stream = csn.arraystream.ArrayStream([trace for trace in stream if trace.stats.distance < 100])
    stream.trim(starttime=stream[0].stats.starttime + trim[0], endtime=stream[0].stats.starttime + trim[1])
    
    min_starttime = max(tr.stats.starttime for tr in stream)
    max_endtime = min(tr.stats.endtime for tr in stream)
    stream.trim(min_starttime, max_endtime)
    stream.resample(40)

    stream.detrend("demean")
    stream.preprocess(window_duration_sec=300, epsilon=1e-10)
    stream.taper(max_percentage=0.01)

    min_starttime = max(tr.stats.starttime for tr in stream)
    max_endtime = min(tr.stats.endtime for tr in stream)
    stream.trim(min_starttime, max_endtime)

    # for tr in stream:
    #     if tr.stats.npts > 24000:
    #         print("NPTS :", tr.stats.npts)
    #         stream.remove(tr)


    print(stream)

    print("Nombre de traces dans le stream", len(stream))

    return stream


def get_entropy(stream, window_size=50, average=4):
    """ Explique ce que chaque ligne fait et peut-être, explique-le en anglais """
    times, frequencies, covariances = csn.covariancematrix.calculate(stream, window_size, average)
    entropy = covariances.coherence(kind="entropy")
    return times, frequencies, entropy


# def compute_neguentropie(stream, av0):
#     norm = min(len(stream), av0)
#     shanon_index = np.exp(entropy0)
#     entropy0 = (np.exp(entropy0)-1)/norm
#     neguentropie = 1 - entropy0

#     return neguentropie, shanon_index


def compute_entropy(stream, ws0, av0, ax):
    times0, frequencies0, entropy0 = get_entropy(stream, window_size=ws0, average=av0 - 1)
    
    f = frequencies0[1:]
    frequency_mask = (f > 0.1) & (f < 18) #0.5
    f = f[frequency_mask]
    shanon_index = np.exp(entropy0)
    shanon_index = shanon_index[:, frequency_mask].max(axis=0)
    shanon_index_normalised = shanon_index/len(stream)

    ax[2].plot(f, shanon_index_normalised, lw=0.75)

    return f, shanon_index_normalised


def entropy(ESEC_avalanches, stream_csn, ws0, av0, ax, curve_params, event_index):
    f, shanon_index = compute_entropy(stream_csn, ws0, av0, ax)
    
    #Moyenne glissante sur f et shanon_index
    cp.moyenne_glissante(f, shanon_index, ax, window_size = 40)

    split_freq = analysis2.find_split_frequency(f, shanon_index, min_freq=4.0, max_freq=10.0)

    print("Valeur de la fréquence coin : ", split_freq)

    low_mask = f <= split_freq
    high_mask = f > split_freq

    low_freq, low_slope, low_intercept, low_shanon = analysis2.ajustement_de_segment(low_mask, f, shanon_index, ax[2], color='green', label="Ajustement bas", pltplot = False)

    high_freq, high_slope, high_intercept, high_shanon = analysis2.ajustement_de_segment(high_mask, f, shanon_index, ax[2], color='blue', label="Ajustement haut", pltplot = False)
        
    curve_params.append({'Event Index': event_index,'FC_ent': split_freq, 'Slope_BF_ent': low_slope,'Intercept_BF_ent': low_intercept,
                         'First Value_BF_ent': low_shanon[0],'FC_value_ent': low_shanon[-1],'Slope_HF_ent': high_slope,'Intercept_HF_ent': high_intercept,
                         'Last Value_HF_ent': high_shanon[-1],'type': ESEC_avalanches["type"][event_index], 'numero3': ESEC_avalanches["numero"][event_index], 'volume': ESEC_avalanches["volume"][event_index],
                         'length': ESEC_avalanches["length"][event_index], 'height': ESEC_avalanches["height"][event_index]})


def entropy_total(ESEC_avalanches, trim, ws0, av0, curve_params):
    for event_index in tqdm(ESEC_avalanches["numero"], total=len(ESEC_avalanches)):
        try:
            print("-------------------------")
            print("Event numéro", event_index)
            print("-------------------------")

            stream = obspy.read(f"sismogrammes/{event_index:03d}.pickle")
            stream = stream.sort(keys=["distance"])
            stream = stream.select(component="Z")
            stream.trim(starttime=stream[0].stats.starttime + trim[0], endtime=stream[0].stats.starttime + trim[1])
            analysis2.filtering(stream, freq_HP=9, freq_LP=0.5)
            trace = stream[0]

            stream_csn = filter_stream_with_covseisnet(event_index, trim)

            if len(stream_csn) == 1:
                raise Exception("Il n'y a qu'une seule trace dans le stream. Calculer l'entropie est impossible")
            

            fig, ax = plt.subplots(3, 1, figsize=(10, 8), constrained_layout=True, gridspec_kw=dict(height_ratios=[1, 1, 2]))

            ax[0].plot(trace.times(), trace.data)
            ax[0].set_ylabel("Nombres numériques")
            ax[0].set_xmargin(0)
            ax[0].set_xticklabels([])

            ax[1].plot(stream_csn[0].times(), stream_csn[0].data)
            ax[1].set_xmargin(0)
            ax[1].set_xlabel("Temps [s]")
            ax[1].set_ylabel("Nombres numériques")

            entropy(ESEC_avalanches, stream_csn, ws0, av0, ax, curve_params, event_index)

            ax[2].legend()
            ax[2].set_xscale("log")
            ax[2].set_ylabel("Entropie")
            ax[2].set_xlabel("Fréquences [Hz]")
            ax[2].set_xmargin(0)

            # plt.title(ESEC_avalanches["type"][event_index])
            figures.save(f"features/3_entropie/pictures/entropie_{event_index}.pdf", tight_layout=False)
            plt.show()

        except Exception as e:
            print('Pas de trace dans le stream', e)

            ESEC_avalanches = ESEC_avalanches.drop(event_index)

    return curve_params, ESEC_avalanches














"""

#===================================================================================================
# Function to read a pickel file and return a stream
# input:
# pickelFile: the pickel file to read
# wd: the window duration used to preprocess the stream
# trim: the time window to trim the stream
# mpc: the maximum percentage of the taper
# output:
# stream: the stream read from the pickel file
def get_trace(pickelFile,wd=250,trim=[300,800],mpc=0.1):
    stream = csn.read(pickelFile)
    stream.preprocess(window_duration_sec=wd)
    stream.taper(max_percentage=mpc)
    stream.trim(starttime=stream[0].stats.starttime + trim[0], endtime=stream[0].stats.starttime + trim[1]);
    return stream

#===================================================================================================
# Function to calculate the entropy of a stream
# input:
# stream: the stream to calculate the entropy
# window_size: the window size used to calculate the entropy
# average: the number of windows used to average the entropy
# output:
# times: the time axis of the entropy
# frequencies: the frequency axis of the entropy
# entropy: the entropy
def get_entropy(stream, average_step, window_size=50, average=4):
    times, frequencies, covariances = csn.covariancematrix.calculate(stream, window_size, average, average_step)
    entropy = covariances.coherence(kind="entropy")
    return times, frequencies, entropy


#===================================================================================================
# Function to plot the entropy
#
# The plot is divided in two subplots: the first one shows the trace and the second one the entropy
# The entropy is plotted as a pcolormesh plot
# input:
# stream: the stream to plot
# times: the time axis of the entropy
# frequencies: the frequency axis of the entropy
# entropy: the entropy to plot
# window_size: the window size used to calculate the entropy
# average: the number of windows used to average the entropy
def plot_entropy(stream,times, frequencies, entropy, window_size=50, average=4):
    # Plot the entropy

    # #user Helvetica font for plots
    # plt.rcParams["font.family"] = "Helvetica"
    # plt.rcParams["font.size"] = 10
    # plt.rcParams["mathtext.fontset"] = "custom"
    # plt.rcParams["mathtext.rm"] = "Helvetica"
    # plt.rcParams["mathtext.it"] = "Helvetica:italic"
    # plt.rcParams["mathtext.bf"] = "Helvetica:bold"
    
    fig, ax = plt.subplots(2, 1, sharex=True, constrained_layout=True,dpi=150)
    # Plot the trace and entropogram
    ax[0].plot(stream[0].times(), stream[0].data, color="k") # Plot only the first trace
    mappable = ax[1].pcolormesh(times, frequencies, entropy.T, cmap="RdYlBu", shading="flat", rasterized=True)
    # Label the plot
    ax[0].set_ylabel("DN [*]")
    ax[0].grid()
    ax[1].set_ylim(frequencies[1], frequencies.max() / 2)
    ax[1].set_yscale("log")
    ax[1].set_ylabel("Frequency (Hz)")
    ax[1].set_xlabel("Time (seconds)")
    ax[1].title.set_text("Entropy (window size: {}s, average: {})".format(window_size, average))
    plt.colorbar(mappable, ax=ax[1], label="Entropy");

"""