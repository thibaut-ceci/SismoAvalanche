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
    """
    Adds labels to the scatter plot grid's x and y axes.

    Parameters:
    -----------
    x_labels : list of str
        A list of labels for the x-axis for each subplot column.
    y_labels : list of str
        A list of labels for the y-axis for each subplot row.
    X : int
        The number of rows in the scatter matrix.
    Y : int
        The number of columns in the scatter matrix.
    axs : np.ndarray
        A 2D array of Axes objects representing the grid of subplots.
    fontsize : int
        The font size to use for the labels.
    """

    ## Loop through y_labels and assign them as the y-axis labels for the first column of the grid
    for i, label in enumerate(y_labels):
        axs[i + 1, 0].set_ylabel(label, fontsize=fontsize)

    ## Loop through x_labels and assign them as the x-axis labels for the last row of the grid
    for j, label in enumerate(x_labels):
        axs[Y - 1, j].set_xlabel(label, fontsize=fontsize)


def plot_results(columns, axs):
    """
    Plots scatter plots for each pair of features in a grid of subplots.

    Parameters:
    -----------
    columns : list or pandas.DataFrame
        A collection of data columns (features) to be plotted against each other.

    axs : numpy.ndarray
        A 2D array of Axes objects representing the grid of subplots.
    """

    ## Compute the number of columns or features
    num_features = len(columns)  

    ## Loop through every possible combination of features
    for i in range(num_features):
        for j in range(num_features):

            ## Skip the diagonal
            if i == j:
                continue 

            ## Plot the scatter plot for lower triangle
            if i > j:
                axs[i, j].plot(columns[j], columns[i], 'o', color="C0")

            ## Plot the scatter plot for upper triangle
            else:
                axs[j, i].plot(columns[i], columns[j], 'o', color="C0")


def histogram(features, axs, bin=10, fontsize=6):
    """
    Plots a histogram of the given features (in log).

    Parameters:
    -----------
    features : np.array
        The features.
    axs : matplotlib.axes.Axes
        The Axes object on which the histogram will be plotted.
    bin : int
        The number of bins to use for the histogram.
    fontsize : int
        Font size for the y-axis label.
    """
    ## Define logarithmically spaced bin edges 
    bin_edges = np.logspace(np.log10(min(features)), np.log10(max(features)), num=bin)

    ## Plot the histogram in log
    axs.hist(features, bin_edges, color='C0')

    ## Create an axis to the right
    ax = axs.twinx()
    ylim0 = axs.get_ylim()

    ax.set_ylim(ylim0)
    ax.set_ylabel("Number of events", rotation=270, va='bottom', fontsize=fontsize)


def histogram_normal(features, axs, fontsize=6):
    """
    Plots a histogram of the given features.

    Parameters:
    -----------
    features : np.array
        The features.
    axs : matplotlib.axes.Axes
        The Axes object on which the histogram will be plotted.
    fontsize : int
        Font size for the y-axis label.
    """
    ## Plot the histogram
    axs.hist(features, 20, color='C0')

    ## Create an axis to the right
    ax = axs.twinx()
    ylim0 = axs.get_ylim()
    
    ax.set_ylim(ylim0)
    ax.set_ylabel("Number of events", rotation=270, va='bottom', fontsize=fontsize)


def histogram_flip(features, axs, fontsize=6):
    """
    Plots a histogram of the given features and flip it.

    Parameters:
    -----------
    features : np.array
        The features.
    axs : matplotlib.axes.Axes
        The Axes object on which the histogram will be plotted.
    fontsize : int
        Font size for the y-axis label.
    """

    ## Plot the histogram
    axs.hist(features, orientation='horizontal')

    ## Create an axis to the right
    ax = axs.twiny()
    ylim = axs.get_xlim()
    ax.set_xlim(ylim)
    ax.set_xlabel("Number of events", rotation=0, va='top', fontsize=fontsize)
    ax.xaxis.tick_bottom()

    ## Flip the histogram
    ax.xaxis.set_label_position('bottom')
    ax.xaxis.set_label_coords(0.5,-0.15)


def set_axe_log(axs, size, y_log_rows, x_log_cols):
    """
    Sets specific axes to logarithmic scale based on the row or column index.

    Parameters:
    -----------
    axs : 2D array of matplotlib.axes.Axes
        The axes on which to set the log scale.
    size : int
        The total number of axes in the grid.
    y_log_rows : list of int
        The row indices for which the y-axis should be set to log scale.
    x_log_cols : list of int
        The column indices for which the x-axis should be set to log scale.
    """

    ## Apply log scale on the y-axis
    for i in range(size):
        if i in y_log_rows:
            for j in range(size):
                axs[i, j].set_yscale('log')

                ## Keep y ticks on the left-most column
                if j != 0:
                    axs[i, j].tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

    ## Apply log scale on the x-axis
    for j in range(size):
        if j in x_log_cols:
            for i in range(size):
                axs[i, j].set_xscale('log')

                ## Keep x ticks on the bottom-most row
                if i != size - 1:
                    axs[i, j].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

