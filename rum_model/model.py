import numpy as np
import sys
import pandas as pd

import rum_model.parameters as parameters
import rum_model.utils as utils

def TT_model_rum(
    symptomatic_ascertainment_rate, 
    symptomatic_rate,
    percentage_notified, 
    compliance_with_symptom_isolation_test,
    compliance_with_symptom_isolation_no_test,
    compliance_with_contact_isolation,
    serial_interval,
    secondary_symptoms_to_tertiary,
    time_to_tertiary_infection,
    total_time_to_contact
    ):
    '''
    Runs the Rum model as defined in the technical document

    Parameters
    ----------
    symptomatic_ascertainment_rate : float (between 0 and 1)
        The probability a symptomatic person is tested
    symptomatic_rate : float (between 0 and 1)
        The probability an infected person is symptomaitc
    percentage_notified : float (between 0 and 1)
        The proportion of infected contacts that are successfully contacted
    compliance_with_symptom_isolation_test : float (between 0 and 1)
        Probability an individual isolates on symptoms given that they are tested
    compliance_with_symptom_isolation_no_test : float (between 0 and 1)
        Probability an individual isolates on symptoms given that they are not tested
    compliance_with_contact_isolation : float (between 0 and 1)
        Probability an individual isolates on contact
    serial_interval : DataFrame
        Distribution of the time from sympom onset in one case to symptom onset in the individual they infect
    secondary_symptoms_to_tertiary : DataFrame
        Distribution of the time from symptom onset in the secondary case to infection of the tertiary case
    time_to_tertiary_infection : DataFrame
        Time from sympom onset in the primary case to infection of the tertiary case
    total_time_to_contact : DataFrame
        Time from sympom onset in the primary case to contact of the secondary case

    Returns
    -------
    contributions_dict : dictionary
        A dictionary containing the probability of the various transmission events
    '''


    ## SYMPTOM METRICS ##

    # Probability case isolates on symptoms - taking into account isolating on test or no test
    adherence_symptom_isolation = symptomatic_rate * \
                                          symptomatic_ascertainment_rate * \
                                          compliance_with_symptom_isolation_test + \
                                          (1 - symptomatic_ascertainment_rate) * \
                                          symptomatic_rate * \
                                          compliance_with_symptom_isolation_no_test
    # probability that if the case isolated, it occured before the secondary case was infected
    adherence_symptom_isolation_impact = secondary_symptoms_to_tertiary[secondary_symptoms_to_tertiary['delay'] >= 0]['frequency'].sum()
    # probability that the case isolated and did no before onward infection
    symptom_isolation_success = adherence_symptom_isolation_impact*adherence_symptom_isolation

    # TEST METRICS ## 

    # Probability the primary case was tested (assuming only symptomatic cases are test)
    primary_tested = symptomatic_rate * symptomatic_ascertainment_rate

    ## CONTACT METRICS ##

    # Probability that the secondary case was contacted
    proportion_contacts_reached = primary_tested * percentage_notified
    # Probability that the secondary case isolated on contact
    proportion_contacts_reached_compliant = proportion_contacts_reached * compliance_with_contact_isolation
    # Probability that if the secondary case isolated, that they did so before onward transmission occurred
    # i.e. the probability that the time between symptom onset in the primary case and 
    # infection of the tertiary case is more than to time between symptom onset and the 
    # secondary case being contacted.
    adherence_contact_isolation_impact = parameters.get_contact_isolation_impact(time_to_tertiary_infection.copy(), total_time_to_contact.copy())

    proportion_transmission_pre_contact = (1 - adherence_contact_isolation_impact)
    transmission_occuring_pre_symptom = (1 - adherence_symptom_isolation_impact)


    # probability that isolation on contact was successful
    contact_isolation_success = proportion_contacts_reached_compliant*adherence_contact_isolation_impact

    # the intersection of symptom of contact success
    symptom_and_contact_success = parameters.get_symptom_and_contact_success(
        symptom_isolation_success,
        symptomatic_rate,
        symptomatic_ascertainment_rate,
        compliance_with_contact_isolation,
        percentage_notified,
        secondary_symptoms_to_tertiary.copy(),
        serial_interval.copy(),
        total_time_to_contact.copy()
    )


    # the overall success is the probability
    # that either symptom isolation is successful or contact isolation is successful
    # we must deduct the intersection of symptom and contact success
    # as the two events are not disjoint. that is:
    # P(A OR B) = P(A) + P(B) - P(A AND B)
    overall_success = symptom_isolation_success + \
                                contact_isolation_success - \
                                symptom_and_contact_success

    # transmission averted in the percetnage of chains that are broken
    # and is equivalent to the overall success

    # the marginal impact of conatct tracing
    marginal_impact = overall_success - symptom_isolation_success

    # write outputs to dictionary
    # the main outputs of interest are transmission averted and marginal impact
    # but the other intermediate outputs are included for reference
    # For example, it is useful to know how much transmission ocurred pre-symptom
    # and so what is the 'maximum' amount of tranmission that can be prevented
    # by symptom isolation
    contributions_dict = {
        'primary_tested' : primary_tested,
        'adherence_symptom_isolation' : adherence_symptom_isolation,
        'adherence_symptom_isolation_impact' : adherence_symptom_isolation_impact,
        'symptom_isolation_success' : symptom_isolation_success,
        'proportion_contacts_reached' : proportion_contacts_reached,
        'proportion_contacts_reached_compliant' : proportion_contacts_reached_compliant,
        'adherence_contact_isolation_impact' : adherence_contact_isolation_impact,
        'contact_isolation_success' : contact_isolation_success,
        'transmission_pre_contact' : proportion_transmission_pre_contact,
        'transmission_occuring_pre_symptom' : transmission_occuring_pre_symptom,
        'transmission_averted' : overall_success,
        'marginal_impact' : marginal_impact,
        'symptom_and_contact_success' : symptom_and_contact_success
        }


    return contributions_dict

