import pandas as pd
import numpy as np
import os
import pyodbc
import math
from scipy.interpolate import interp1d
import datetime
import time


def process_delay_column(text):
    '''
    The delay column in the target distribution is written as an interval
    Process this function to it calculates the mid_point of the interval

    Parameters
    ----------
    text : str
        text in form (a, b] to split
    
    Returns
    -------
    mid_point : float
        mid point of a and b
    '''
    
    end_points = [int(float(x)) for x in text.strip("(]").split(', ')]
    
    mid_point = (end_points[0] + end_points[1])/2.
    
    return mid_point

def process_TT_e2e_file(input_dir, file_name):
    '''
    Process the contact tracing distribution file
    Ensure they are normalised and that the distributions
    is defined an interget hourly intervals

    Parameters
    ----------
    input_dir : str
        name of directory in which TT distribution is stored
    file_name : str
        name of file containing TT distribution
    
    Returns
    -------
    e2e_df : DataFrame
        dataframe with frequency column and delay column. Delay column
        specifies the hour (at every hour) and frequency is the pdf
    '''

    e2e_df = pd.read_csv(os.path.join(input_dir, file_name), index_col = 0).transpose().reset_index()
    e2e_df.columns = ['delay', 'frequency']
    e2e_df['delay'] = e2e_df['delay'].apply(process_delay_column)

    # interpolate so the frequency is defined at hourly intervals
    x = e2e_df['delay']
    x_new = np.arange(math.ceil(min(x)), math.floor(max(x)), 1)
    f = interp1d(x, e2e_df['frequency'], kind = 'cubic')
    freq_new = f(x_new)

    # normalise frequency
    freq_new /= freq_new.sum()

    e2e_df = pd.DataFrame({'delay' : x_new, 'frequency' : freq_new})
    e2e_df['frequency'] = e2e_df['frequency']/e2e_df['frequency'].sum()

    return e2e_df

def pad_with_zeros(df_list, min_delay_bound):
    '''
    Pads every dataframe in df_list with zeros to ensure that they are of the same length
    and the delay is symmetric about zero

    Parameters
    ----------
        df_list: list
            list of dataframes
        min_delay_bound: integer
            minimum (half-)length of delay. This should be 40*24 as the epi parameters
            are defined across the interval (-40*24,40*24), where 40 represents the days
    
    Returns
    -------
        new_df_list: list
            list of padded dataframes in same order as df_list
    '''

    # standardise the delay array across all dataframes
    # the delay must be symmetric about zero for convolution
    min_delay = math.floor(min([df['delay'].min() for df in df_list]))
    max_delay = math.ceil(max([df['delay'].max() for df in df_list]))
    delay_bound = max(max_delay, -min_delay, min_delay_bound)
    delay = pd.DataFrame({'delay' : np.arange(-delay_bound, delay_bound + 1, 1)})

    # pad the dataframe with the expanded delay column
    new_df_list = []
    for df in df_list:
        df = delay.merge(df, on = 'delay', how = 'left').fillna(0)
        new_df_list.append(df)

    return new_df_list