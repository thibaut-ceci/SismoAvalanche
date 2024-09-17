import matplotlib.pyplot as plt
import numpy as np


def correlation_matrix(data, label, print_data=False):
    # Compute correlation matrix
    corr_matrix = data.corr()

    plt.imshow(corr_matrix, cmap='hot', interpolation='nearest')
    plt.colorbar(label="Corrélation")

    plt.xticks(range(data.shape[1]), label, rotation=45)
    plt.yticks(range(data.shape[1]), label)

    #Print the correlaton matrix for see the values
    if print_data == True:
        print(corr_matrix)


def remove_fig(fig, axs):
    num_features = axs.shape[0]

    for i in range(num_features):
        for j in range(i + 1, num_features):
            fig.delaxes(axs[i, j])


def remove_ticks(X, Y, axs):
    for i in range(X):
        for j in range(Y):
            if i == j:
                axs[i, j].xaxis.set_ticklabels([])
                axs[i, j].yaxis.set_ticklabels([])
                axs[i, j].set_xticks([])
                axs[i, j].set_yticks([])

            if j == 0 and i < 8:
                if not (i == X - 1 and j == 0):  # Keep the x-ticks for the bottom-left figure
                    axs[i, j].xaxis.set_ticklabels([])
                    axs[i, j].set_xticks([])

            if j > 0 and (i != j or j == 8):
                if i != X - 1:  # Skip the last row
                    axs[i, j].xaxis.set_ticklabels([])
                    axs[i, j].set_xticks([])
                axs[i, j].yaxis.set_ticklabels([])
                axs[i, j].set_yticks([])


def label_size(X, Y, axs, labelsize = 14):
    for i in range(X):
        axs[i, 0].tick_params(axis='both', which='major', labelsize=labelsize)

    for j in range(1, Y):
        axs[Y-1, j].tick_params(axis='both', which='major', labelsize=labelsize)


def labels_scatter(x_labels, y_labels, X, Y, axs, fontsize=25):   
    for i, label in enumerate(y_labels):
        axs[i + 1, 0].set_ylabel(label, fontsize=fontsize)

    for j, label in enumerate(x_labels):
        axs[Y - 1, j].set_xlabel(label, fontsize=fontsize)


def plot_results(columns, axs):
    num_features = len(columns)  # The length of columns

    for i in range(num_features):
        for j in range(num_features):
            if i == j:
                continue
            if i > j:
                axs[i, j].plot(columns[j], columns[i], 'o', color="C0")
            else:
                axs[j, i].plot(columns[i], columns[j], 'o', color="C0")


def histogram(features, axs, bin=10, fontsize=6):
    bin_edges = np.logspace(np.log10(min(features)), np.log10(max(features)), num=bin)
    axs.hist(features, bin_edges, color='C0')

    ax1 = axs.twinx()

    ylim0 = axs.get_ylim()
    ax1.set_ylim(ylim0)
    ax1.set_ylabel("Nombre d'événements", rotation=270, va='bottom', fontsize=fontsize)


def histogram_normal(features, axs, bin=10, fontsize=6):
    axs.hist(features, 20, color='C0')

    ax1 = axs.twinx()

    ylim0 = axs.get_ylim()
    ax1.set_ylim(ylim0)
    ax1.set_ylabel("Nombre d'événements", rotation=270, va='bottom', fontsize=fontsize)


def histogram_flip(features, axs, bin=10, fontsize=6):
    #bin_edges = np.logspace(np.log10(min(features["hl"])), np.log10(max(features["hl"])), num=20)
    axs.hist(features, bin, orientation='horizontal')

    ax5 = axs.twiny()

    ylim = axs.get_xlim()
    ax5.set_xlim(ylim)
    ax5.set_xlabel("Nombre d'événements", rotation=0, va='top', fontsize=fontsize)
    ax5.xaxis.tick_bottom()
    ax5.xaxis.set_label_position('bottom')
    ax5.xaxis.set_label_coords(0.5,-0.15)


def set_axe_log(axs, size):
    for i in range(0, size):
        axs[2, i].set_yscale('log')
        axs[3, i].set_yscale('log')
        axs[5, i].set_yscale('log')
        axs[6, i].set_yscale('log')
        axs[7, i].set_yscale('log')
        axs[8, i].set_yscale('log')
        axs[9, i].set_yscale('log')
        axs[11, i].set_yscale('log')
        axs[14, i].set_yscale('log')
        axs[16, i].set_yscale('log')

    for i in range(0, size):
        axs[i, 0].set_xscale('log')
        axs[i, 2].set_xscale('log')
        axs[i, 3].set_xscale('log')
        axs[i, 5].set_xscale('log')
        axs[i, 6].set_xscale('log')
        axs[i, 7].set_xscale('log')
        axs[i, 8].set_xscale('log')
        axs[i, 9].set_xscale('log')
        axs[i, 11].set_xscale('log')
        axs[i, 14].set_xscale('log')
        axs[i, 16].set_xscale('log')
