# The R&ugrave;m Model

## Model Description

The R&ugrave;m Model estimates the impact of test, trace and isolation on transmission reduction.

We assume an infected population, whom we refer to as tertiary cases. We then back-populate to find the secondary case, and primary case that led to infection of the tertiary cases. We assume the tertiary cases were infected in a world without any testing, tracing or self-isolating. As a result, there is no dependence on the R number.

For each step in the transmission chain, we estimate the probability of success of isolation, so that onward transmission was prevented. The percentage of transmission chains stopped, is the overall success or transmission reduction.

As the transmission chains are constructed assuming no self-isolation, the transmission reduction is relative to the counter factual where nobody self-isolates. This must be considered when interpreting the results.

## Instructions

### Installation

Clone this repository, and enter the repo directory.

Requirements for this model are quite basic. The package versions on which this model was built are in the [requirements.txt](requirements.txt) file.

To install the requiremnets:
* pip install -r requirements.txt

### Config files

All config files must be stored in the [configs/](configs/) directory.

We suggest you use the [config file](configs/config.json) provided. If you wish to create a new model run, create a new config file with a different name.

* config.json is an example of a config file that would produce the outputs in the technical document, given the right scenario settings and contact tracing distributions

In order to run the R&ugrave;m model, the necessary model run config inputs are:

| Key Name | Description |
| -------- | -----------|
| infectiousness_option | Either 'ashcroft' or 'he'. Determines which symptom to onward vector to use. Read the technical document for further details.|
| total_time_to_contact_file | Name of the file containing the distribution of end to end contact time.|
| target_total_time_to_contact_file |  Same as above but for the target scenario. |
| scenarios_file | Name of csv file containing scenarios. **Set to scenarios.csv as default**|
| output_file | name of csv file to store the model outputs |


### Run the R&ugrave;m model

* python run.py model_config_file_name
* For example, to run using the provided config: python run.py config.json

The config file must be located in the [configs](configs/) directory. The [config.json](configs/config.json) reproduces outputs presented in the technical document.

## Data Input

The contact tracing delay distributions and scenario settings defining the T&T parameters are required as input

### Delay distributions

You must provide a contact tracing distribution, and a desired (or 'target') contact tracing distribution. We store these in the [data/input/TT]([data/input/TT]) directory, but you can specify another directory in the config file.

Each of the distributions must be for the time from symptom onset in the primary case, to contact of the tertiary case. More details of how this is calculated in practice are provided in the technical document.

We include the distribution files used in the technical document. This are defined over hourly bins.

| File Name | Description |
| --------- | ----------- |
| Combined_E2E_from_symptoms.csv | Actual distribution for October from symptom onset to contact made |
| Combined_Target_E2E.csv | Target distribution from symptom onset to contact made |


### Scenario files - T&T parameter settings

We provide an example scenario file: [data/input/scenarios/scenarios.csv](data/input/scenarios/scenarios.csv). **Do not change this scenario file**. This file is the approved range of scenarios and reproduces the output in the technical annex.

If you wish to create a new scenarios file, create a new file with a different name. It must be a csv and include the following columns with names as prescribed:

| Column Name | Description |
| ------------| ------------ |
| Scenario | Name of the scenario to be modelled - is an identifier |
| symptomatic_rate | Percentage of infected individuals who are symptomatic at any point |
| symptomatic_ascertainment_rate | Percentage of symptomatic individuals that are tested |
| percentage_notified | Percentage of people who an individual infects that are contacted by test and trace |
| compliance_with_symptom_isolation_test | Probability of self-isolating on symptoms given that they are tested |
| compliance_with_symptom_isolation_no_test | Probability of self-isolating on symptoms without ordering a test |
| compliance_with_test_result_isolation | Probability of self-isolating on a test result, without isolating on symptoms |
| compliance_with_contact_isolation | Probability of self-isolating on contact |
| Distribution | Whether to use the actual October test and trace delay distributions, or the hypothetical test and trace delay dsitrbutions. Should take the value of 'October Distribution' or 'Target Distribution' 

Each row describes a different scenario.

### Output

Outputs are stored in the directory you specify in the config file. The important R&ugrave;m outputs are as follows:

| Output Name | Description |
| ----------- | ----------- |
| transmission_averted | The percentage of transmission chains stopped by isolating on symptoms or contact |
| marginal_impact | The marginal impact of contact tracing |
| symptom_isolation_success | The proportion of chains stopped by isolation on symptoms |
| contact_isolation_success | The proportion of transmission chains stopped by contact tracing |


## Project structure

* [rum_model/model.py](rum_model/model.py) contains the model function
* [rum_model/parameters.py](rum_model/parameters.py) contains functions to calculate necessary contact tracing or epidemiological parameters, and event probabilties.
* [rum_model/utils.py](rum_model/utils.py) provides utility function for processing distributions for example
* [run.py](run.py) is the main script to run through all scenarios and model settings.

## LICENSE

This publication is licensed under the terms of the Open Government Licence v3.0 except where otherwise stated. To view this licence, visit <https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3>
