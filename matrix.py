"""
ESEC scatter matrix

Library to create a scatter matrix.
"""

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def correlation_matrix(data, labels, annotate=True):
    """
    Computes and plots the correlation matrix of multiple columns in ESEC.

    Parameters:
    ------------
    data : pandas.DataFrame
        The columns from the ESEC.
    labels : list of str
        List of column labels to use for x and y axes.
    annotate : bool
        If True, annotates the heatmap with the correlation coefficients.
    """
    ## Compute the correlation matrix
    corr_matrix = data.corr()

    ## Display it
    plt.figure(figsize=(10, 8))
    sns.heatmap(corr_matrix, annot=annotate, cmap='coolwarm', xticklabels=labels, yticklabels=labels, cbar_kws={'label': 'Correlation'}, square=True)

    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()

    plt.show()


def remove_fig(fig, axs):
    """
    Removes the upper triangular subplots from the scatter matrix.

    Parameters:
    -----------
    fig : matplotlib.figure.Figure
        The figure object from which subplots will be removed.
    
    axs : numpy.ndarray
        A 2D array of Axes objects in the figure.
    """
    
    ## Determine the number of features (subplots along one axis)
    num_features = axs.shape[0]

    ## Loop through each subplot in the upper triangular part of the grid
    for i in range(num_features):
        for j in range(i + 1, num_features):

            ## Remove the subplot
            fig.delaxes(axs[i, j])


def remove_ticks(X, Y, axs):
    """
    Removes tick labels and ticks from certain subplots in a grid.

    Parameters:
    -----------
    X : int
        The number of rows in the scatter matrix.
    Y : int
        The number of columns in the scatter matrix.
    axs : numpy.ndarray
        A 2D array of Axes objects representing the subplots in the figure.
    """

    ## Loop through each subplot in the grid
    for i in range(X):
        for j in range(Y):

            ## Remove both x and y tick labels and ticks for the diagonal subplots
            if i == j:
                axs[i, j].xaxis.set_ticklabels([])  ## Remove x-axis tick labels
                axs[i, j].yaxis.set_ticklabels([])  ## Remove y-axis tick labels
                axs[i, j].set_xticks([])            ## Remove x-axis ticks
                axs[i, j].set_yticks([])            ## Remove y-axis ticks

            ## For subplots in the first column (j == 0) remove x-ticks except the bottom-left one
            if j == 0 and i < 8:
                if not (i == X - 1 and j == 0):
                    axs[i, j].xaxis.set_ticklabels([])
                    axs[i, j].set_xticks([])

            ## For all subplots except those on the diagonal (i != j)
            if j > 0 and (i != j or j == 8):
                if i != X - 1:
                    axs[i, j].xaxis.set_ticklabels([])
                    axs[i, j].set_xticks([])
                axs[i, j].yaxis.set_ticklabels([])
                axs[i, j].set_yticks([])


def label_size(X, Y, axs, labelsize=14):
    """
    Adjusts the size of tick labels for the scatter matrix.

    Parameters:
    -----------
    X : int
        The number of rows in the scatter matrix.
    Y : int
        The number of columns in the scatter matrix.
    axs : numpy.ndarray
        A 2D array of Axes objects representing the subplots in the figure.
    labelsize : int
        The size of the labels.
    """

    ## Loop through each row of subplots and set label size for the first column
    for i in range(X):
        axs[i, 0].tick_params(axis='both', which='major', labelsize=labelsize)

    ## Loop through the last row of subplots and set label size for the rest of the columns
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
        axs[5, i].set_yscale('log')
        axs[6, i].set_yscale('log')
        axs[12, i].set_yscale('log')

    for i in range(0, size):
        axs[i, 0].set_xscale('log')
        axs[i, 2].set_xscale('log')
        axs[i, 5].set_xscale('log')
        axs[i, 6].set_xscale('log')
        axs[i, 12].set_xscale('log')

