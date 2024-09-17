"""ESEC catalog management.

This module contains functions to read and filter the catalog obtained by the
ESEC project. A REFAIIIIIIIIIIIIIIIIIIIIIIIRE
"""

import pickle

import pandas as pd

def display_params_for_catalog():
    # Display all columns and rows of the DataFrame. For use it, just import pandas as pd and call this function.
    pd.set_option("display.max_rows",None)
    pd.set_option("display.max_columns",None)
    pd.set_option("display.width",None) 
    pd.set_option("display.max_colwidth",None) 


def load(path):
    return pickle.load(open(path, "rb"))


def get_the_most_or_least_instrumented_avalanches(most_instrumented_avalanches, number, ascending=False):
    #Sort the column "Number of stations" in the descending ordre
    most_instrumented_avalanches.sort_values(by="Number of stations", ascending=ascending, inplace=True)
    most_instrumented_avalanches = most_instrumented_avalanches.head(number)

    return most_instrumented_avalanches