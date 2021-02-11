import numpy as np
import pandas as pd
import os
import sys

from scipy.interpolate import interp1d
import math
from scipy.stats import gamma

import rum_model.utils as utils

def dgamma(x,shape,rate=1):
    """
    Calculates the density/point estimate of the Gamma-distribution
    It is formatted in the same way as the dgamma function in R

    Parameters
    ----------
    x : float array
        array to return the point density estimates of
    shape : float
    rate : float

    Returns
    -------
    result: float array
        probability density function defined at points specified in x

    """

    result=rate*gamma.pdf(x=rate*x,a=shape,loc=0,scale=1)

    return result

def get_incubation_period(delay_bound):
    '''
    Produces a distribution for the time from infection to symptom onset.
    The distribution is assumed to follow a gamma distribution with shape = 4.23 and rate = 0.81
    The distribution is defined at hourly intervals and is padded with zeros to match the 
    shape of other distributions in the model.

    Paramters
    ---------
    delay_bound : float
        maximum (or absolute value of minimum) delay time of the test and trace distributions in hours

    Returns
    -------
    incubation_period : dataframe
        dataframe with columns delay, defining the delay time in hourly intervals, 
        and frequency, defining the point probability density.
    '''

    max_coord = 40 # the maximum day to calculate the gamma pdf for
    delay_days = np.arange(0, max_coord+0.001, 1/24) # array of the delay in days
    freq = dgamma(delay_days, shape = 4.23, rate = 0.81)
    freq /= freq.sum() # normalise the pdf
    delay = np.arange(0,max_coord*24 + 1, 1) # transform the delay in days to hours
    incubation_period = pd.DataFrame({'delay' : delay, 'frequency': freq})

    # pad the dataframe out with zeros to match the length of the test and trace distributions
    delay_padded = pd.DataFrame({'delay' : np.arange(-delay_bound, delay_bound + 1, 1)})
    incubation_period = delay_padded.merge(incubation_period, how = 'left', on = 'delay').fillna(0)

    return incubation_period

def get_serial_interval(delay_bound):
    '''
    Produces a distribution for the time from symptom onset in the primary case 
    to symptom onset in the secondary case
    The distribution is assumed to follow a gamma distribution with shape = 5.18 and rate = 0.96
    The distribution is defined at hourly intervals and is padded with zeros to match the 
    shape of other distributions in the model.

    Paramters
    ---------
    delay_bound : float
        maximum (or absolute value of minimum) delay time of the test and trace distributions in hours

    Returns
    -------
    serial_interval_df : dataframe
        dataframe with columns delay, defining the delay time in hourly intervals, 
        and frequency, defining the point probability density.
    '''

    max_coord = 40 # the maximum day to calculate the gamma pdf for
    delay_days = np.arange(0, max_coord+0.001, 1/24) # array of the delay in days
    freq = dgamma(delay_days, shape = 5.18, rate = 0.96)
    freq /= freq.sum() # normalise the pdf
    delay = np.arange(0,max_coord*24 + 1, 1) # transform the delay in days to hours
    serial_interval = pd.DataFrame({'delay' : delay, 'frequency': freq})

    # pad the dataframe out with zeros to match the length of the test and trace distributions
    delay_padded = pd.DataFrame({'delay' : np.arange(-delay_bound, delay_bound + 1, 1)})
    serial_interval_df = delay_padded.merge(serial_interval, how = 'left', on = 'delay').fillna(0)

    return serial_interval_df

def get_symptom_to_onward_vector(delay_bound, option = 'he'):
    '''
    Produces a distribution for the time between symptom onset in the primary case to 
    onward infection of the secondary case.
    It is calculated by taking the difference between the incubation period, and the serial
    interval vector which is defined as the distribution of time between infection of
    the primary case and infection of the secondary case.

    Where option = 'he':
    The incubation period is assumed to be gamma distribution with shape 4.23 and rate 0.81
    The serial interval vector is assumed to be a gamma distribution with shape 5.18 and rate 0.96

    We also assume that infection is not possible more than 2 days before symptom onset.

    The difference between the two distributions is calculate by taking the convolution

    Where option = 'ashcroft':
    The distribution is taken from the aschroft paper. It is defined daily, so we interpolate
    to obtain the distribution defined at hourly intervals

    Parameters
    ----------
    delay_bound : float
        maximum (or absolute value of minimum) delay time of the test and trace distributions in hours
    option : str
        'he' or 'ashcroft'

    Returns
    -------
    symptom_to_onward_data : dataframe
        dataframe with columns delay, defining the delay time in hourly intervals, 
        and frequency, defining the point probability density.    
    '''

    if option == 'he':

        serial_interval_vector = get_serial_interval(delay_bound)
        incubation_period_vector = get_incubation_period(delay_bound)

        # take difference between the serial interval and incubation period to calculate the symptom to onward distribution
        frequency = np.convolve(serial_interval_vector['frequency'], incubation_period_vector['frequency'][::-1], 'same')
        frequency /= frequency.sum()

        symptom_to_onward_vector = {'delay' : serial_interval_vector['delay'], 'frequency' : frequency}
        symptom_to_onward_data = pd.DataFrame(symptom_to_onward_vector)

        # filter frequency to ensure that transmission doesn't occur more the two days
        # before symptom onset
        symptom_to_onward_data.loc[symptom_to_onward_data['delay'] < -2*24, 'frequency'] = 0
        symptom_to_onward_data['frequency'] /= symptom_to_onward_data['frequency'].sum() # ensure probabilities add up to 1

    elif option == 'ashcroft':

        # infectiousness curve defined delay between -5 and 10 days
        coarse_infectiousness = [
            '0.0223781937259083',
            '0.0353117463128404',
            '0.0662651166344249',
            '0.103030125971041',
            '0.134813059932273',
            '0.150505975534097',
            '0.145114140517949',
            '0.122150668663839',
            '0.0906372826285407',
            '0.0598005731289022',
            '0.0353572882950919',
            '0.0188663868762238',
            '0.00914344380961272',
            '0.00404823605805759',
            '0.00164610895845332',
            '0.000617722128227433'
        ]

        # interpolate to get hourly data
        coarse_delay = [d for d in range(-5*24,11*24,24)]
        fine_delay = np.arange(-5*24,10*24 + 1, 1)

        f = interp1d(coarse_delay, coarse_infectiousness, kind = 'cubic')
        symptom_to_onward = f(fine_delay)

        # ensure normalised
        symptom_to_onward_data = pd.DataFrame({'delay': fine_delay, 'frequency' : symptom_to_onward})
        symptom_to_onward_data['frequency'] /= symptom_to_onward_data['frequency'].sum()

        # pad with zeros
        delay_padded = pd.DataFrame({'delay' : np.arange(-delay_bound, delay_bound + 1, 1)})
        symptom_to_onward_data = delay_padded.merge(symptom_to_onward_data, how = 'left', on = 'delay').fillna(0)

    return symptom_to_onward_data

def calculate_time_distributions(input_dir, config_settings):
    '''
    Function that loads and calculates the necessary contact tracing and epi delay distributions
    depending on the parameters defined in the config file

    Parameters
    ----------
    config_settings : dict
        contains the model run settings contained in the config file.
    
    Returns
    -------
    time_distributions_dict : dict
        contains the epidemiological temporal distributions, and contact tracing distributions
    '''

    # determine whether to use he or ashcroft infectiousness curve
    infectiousness_option = config_settings['infectiousness_option']

    time_distributions_dict = {}
    time_distributions_dict['infectiousness_option'] = infectiousness_option

    e2e_file = config_settings["total_time_to_contact_file"]
    target_e2e_file = config_settings['target_total_time_to_contact_file']

    # process the e2e contact tracing files
    e2e_df = utils.process_TT_e2e_file(input_dir, e2e_file)
    target_e2e_df = utils.process_TT_e2e_file(input_dir, target_e2e_file)

    # pad the distributions with zeros so they are the same size 
    [target_e2e_df, e2e_df] = utils.pad_with_zeros(
        [target_e2e_df, e2e_df],
        40*24
        )

    delay_bound = e2e_df['delay'].max() # this is to ensure that the distribution dataframes are padded with zeros to be the same size

    time_distributions_dict['total_time_to_contact'] = e2e_df.copy()
    time_distributions_dict['target_total_time_to_contact'] = target_e2e_df.copy()

    # Calculate epidemiological parameters #
    ########################################

    # get time from symptoms to onward infection
    symptoms_to_secondary = get_symptom_to_onward_vector(delay_bound, option = infectiousness_option)
    time_distributions_dict['secondary_symptoms_to_tertiary'] = symptoms_to_secondary.copy()
    # get serial interval, which is the time from symptom onset in index case, to symptom onset in the infected case
    symptoms_to_symptoms = get_serial_interval(delay_bound)
    time_distributions_dict['serial_interval'] = symptoms_to_symptoms.copy()
    # time from symptom onset in primary case to infection of tertiary case
    time_to_tertiary_infection = get_time_to_tertiary_infection(
        secondary_symptoms_to_tertiary = symptoms_to_secondary,
        symptoms_to_symptoms = symptoms_to_symptoms
        )
    time_distributions_dict['time_to_tertiary_infection'] = time_to_tertiary_infection.copy()

    return time_distributions_dict


def get_time_to_tertiary_infection(
    secondary_symptoms_to_tertiary,
    symptoms_to_symptoms
    ):
    '''
    Calculates the time from symptom onset in the primary case to
    infection of the tertiary case.

    It is calculated as the sum of:
         the time from symptom onset in the primary case to symptom onset in the secondary case
         and the time from symptom onset in the secondary case to infection of the tertiary case

    Parameters
    ----------
    secondary_symptoms_to_tertiary : dataframe
        Dataframe containining distribution from symptom onset in the secondary case to 
        infection of the tertiary case
    symptoms_to_symptoms : dataframe
        Equivalent to serial interval

    Returns
    ------
    time_to_tertiary_infection : dataframe
        Dataframe with columns : 'delay' defining the time delay in hourly intervals and
        matches the delay arrays of all other test and trace dataframe, and 'frequency' defining
        the point probability density of the time between symptom onset of the primary
        case and infection of the tertiary case.
    '''

    time_to_tertiary = np.convolve(secondary_symptoms_to_tertiary['frequency'], symptoms_to_symptoms['frequency'], 'same')
    time_to_tertiary_infection = pd.DataFrame({'delay' : secondary_symptoms_to_tertiary['delay'], 'frequency' : time_to_tertiary})


    return time_to_tertiary_infection


def get_contact_isolation_impact(
    time_to_tertiary_infection, 
    total_time_to_contact
    ):
    '''
    Returns the probability that the time between symptom onset in the primary case and 
    infection of the tertiary case is more than to time between symptom onset and the 
    secondary case being contacted.

    Parameters
    ----------
    time_to_tertiary_infection : dataframe
        Output of get_time_to_tertiary_infection
    total_time_to_contact : dataframe
        Output of get_toal_time_to_contact
    
    Returns
    -------
    prob : float
        Probability time_to_tertiary_symptoms > total_time_to_contact
    '''

    # subtract total_time_to_contact distribution from time_to_tertiary infection
    # using convolutions
    diff_freq = np.convolve(time_to_tertiary_infection['frequency'], total_time_to_contact['frequency'][::-1], 'same')
    diff_df = pd.DataFrame({'delay' : time_to_tertiary_infection['delay'], 'frequency' : diff_freq})
    prob = diff_df[diff_df['delay'] > 0]['frequency'].sum()

    return prob


def get_symptom_and_contact_success(
    symptom_isolation_success,
    symptomatic_rate,
    symptomatic_ascertainment_rate,
    compliance_with_contact_isolation,
    percentage_notified,
    secondary_symptoms_to_tertiary,
    serial_interval,
    total_time_to_contact
):
    ''' 
    It assumes symptom isolation success occurs in the secondary generation

    Returns the probability that isolation occurred on symptom onset and occurred before onward infection
    AND isolation occurred on contact and occurred before onward infection

    When calculating the probability that contact isolation was impactful
    we must assume that symptom isolation was impactful

    Parameters
    ----------

    Returns
    -------
    test_and_contact_success : float
        Probability
    '''

    # enforce symptom_to_secondary to be positive - as this is the condition for symptomatic
    # isolation to be successful
    secondary_symptoms_to_tertiary.loc[secondary_symptoms_to_tertiary['delay'] < 0 , 'frequency'] = 0
    # re-normalise the frequency
    secondary_symptoms_to_tertiary['frequency'] /= secondary_symptoms_to_tertiary['frequency'].sum()

    time_to_tertiary_infection = get_time_to_tertiary_infection(
        secondary_symptoms_to_tertiary = secondary_symptoms_to_tertiary,
        symptoms_to_symptoms = serial_interval
        )
    conditional_contact_impact = get_contact_isolation_impact(
        time_to_tertiary_infection, total_time_to_contact)
    symptom_and_contact =  symptom_isolation_success * \
                           symptomatic_rate * \
                           symptomatic_ascertainment_rate * \
                           percentage_notified * \
                           compliance_with_contact_isolation * \
                           conditional_contact_impact
                           

    return symptom_and_contact


