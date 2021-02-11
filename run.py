import pandas as pd
import os
import sys
import json
import rum_model.parameters as parameters
import rum_model.model as model


def main(model_settings, directory_settings):
    '''
    Wrapper function that iterates through the scenarios, runs the rum model for each
    and saves the output to file.

    Parameters
    ----
    model_settings : dict
        A dictionary containing the settings for the model run as defined in the config file
        The dictionary must contain:
            "total_time_to_contact_file"
            "target_total_time_to_contact_file"
            "scenarios_file"
            "output_file"
            "infectiousness_option"
    directory_settings: dict
        A dictionary that states where files are stored or should be written
        The dictionary must contain:
            "input_dir"
            "output_dir"
            "scenarios_dir"

    Returns
    -------
    output_df : DataFrame
        A dataframe containing the parameters and output metrics for each scenario or row in the
        scenarios file
    '''

    # read directory names
    input_dir = directory_settings['input_dir']
    scenarios_dir = directory_settings['scenarios_dir']
    output_dir = directory_settings['output_dir']

    # check directories exist.
    # if output dir does not exist, then create the directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for directory in [scenarios_dir, input_dir]:
        if not os.path.exists(directory):
            print('Path: ' + directory + ' does not exist')
            sys.exit()


    # Load file settings
    output_file = model_settings['output_file']
    scenarios_file = model_settings['scenarios_file']

    # Load model settings #
    ####################### 

    ## Load scenarios data to run
    scenarios_df = pd.read_csv(os.path.join(scenarios_dir,scenarios_file))

    # read end to end contact tracing time distributions
    # and calculate the epidemiological delay distributions
    time_distributions_dict = parameters.calculate_time_distributions(input_dir, model_settings)

    # initialise output dataframe
    output_df = pd.DataFrame()

    # RUN RUM MODEL #
    ## Iterate through scenarios in the scnearios dataframe, where each row is a
    ## different T&T scenario
    for index, row in scenarios_df.iterrows():

        # define input parameters
        symptomatic_rate = row['symptomatic_rate']
        symptomatic_ascertainment_rate = row['symptomatic_ascertainment_rate']
        percentage_notified = row['percentage_notified']
        compliance_with_symptom_isolation_test = row["compliance_with_symptom_isolation_test"]
        compliance_with_symptom_isolation_no_test = row['compliance_with_symptom_isolation_no_test']
        compliance_with_contact_isolation = row["contact_compliance_with_isolation"]
        distribution = row["Distribution"]
        scenario_name = row["Scenario"]

        # load the epi delay distributions
        secondary_symptoms_to_tertiary = time_distributions_dict['secondary_symptoms_to_tertiary']
        time_to_tertiary_infection = time_distributions_dict['time_to_tertiary_infection']
        serial_interval = time_distributions_dict['serial_interval']

        # define the scenario e2e contact time options
        # this depends on whether we wish the use the real October distribution
        # or the target distribution as is specified by the scenario
        if distribution == "October Distribution":
            total_time_to_contact = time_distributions_dict['total_time_to_contact']
        else:
            total_time_to_contact = time_distributions_dict['target_total_time_to_contact']

        # define parameter dictionary to keep record of input parameters
        parameter_dict = {
            'Scenario' : scenario_name,
            'symptomatic_ascertainment_rate': symptomatic_ascertainment_rate,
            'symptomatic_rate' : symptomatic_rate,
            'percentage_notified' : percentage_notified,
            'compliance_with_symptom_isolation_test': compliance_with_symptom_isolation_test,
            'compliance_with_symptom_isolation_no_test' : compliance_with_symptom_isolation_no_test,
            'compliance_with_contact_isolation' : compliance_with_contact_isolation,
            'distribution' : distribution
        }

        # run TT model
        # output is a dictionary of output metrics, included the proportion of
        # transmission averted and marginal impact of contact tracing
        output = model.TT_model_rum(
            symptomatic_ascertainment_rate,
            symptomatic_rate,
            percentage_notified,
            compliance_with_symptom_isolation_test,
            compliance_with_symptom_isolation_no_test,
            compliance_with_contact_isolation,
            serial_interval.copy(),
            secondary_symptoms_to_tertiary.copy(),
            time_to_tertiary_infection.copy(),
            total_time_to_contact.copy(),
        )

        # print output
        print(scenario_name)
        print(output)
        print('-------------------')

        # combine parameters and outputs and add a new row to the output dataframe
        output_df = output_df.append({**parameter_dict, **output}, ignore_index = True)

    # save to file
    output_df.to_csv(os.path.join(output_dir, output_file))

    return output_df

if __name__ == "__main__":

    # Check that the user has supplied name of config file
    n_args = len(sys.argv)

    if n_args != 2:
        print("Provide the name of a config file : python run.py config_file_name")
        sys.exit()

    config_name = sys.argv[1]

    ## Load config file and define working directory
    config_dir = 'configs'
    config_file = os.path.join(config_dir, config_name)
    with open(config_file) as json_file:
        config_settings = json.load(json_file)

    ## Load directory settings    
    directory_settings = config_settings['directory_settings']

    # Iterate through different model runs as defined in the config file.
    # If using our specifed config file, we are just iterating over
    # the two infectiousness methodologies
    for model_name in config_settings["model_runs"]:
        print("RUNNING MODEL: ", model_name)
        main(config_settings["model_runs"][model_name], directory_settings)

    
    


