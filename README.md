# SINEFIA
A wrapper to run and store CLOUDY models, and compare with observations through a Bayesian approach.

## Introduction
CLOUDY is a commonly used radiative transfer code (https://trac.nublado.org/) that simulates the physical conditions and emission spectra of a cloud of gas and dust under various astrophysical environments. The tools we present here are designed to facilitate running CLOUDY models for starburst and AGN environments, and comparing the results with observations. We also aim to procedurally build a database of line emission from the most common emission lines as a function of different CLOUDY parameters. Users can immediately utilize results from the available lines without needing to run the often time-consuming and computationally expensive CLOUDY models. 

As the architecture of the models and the parameters involved differ depending on whether an AGN or a starburst (SB) is considered as the main source of ionization, the two are treated separately. However, the structure and functions of the codes remain practically identical. 

## Brief description
When running the Python scripts, a series of parameters require user input, all of which are indicated and described in the self-explanatory codes. A binary switch (cloudy_switch) triggers the running of CLOUDY models: users define values for the cloud and ionizing source parameters (e.g., density, metallicity, ionization parameter), with the option to run a grid of models for some parameters. A second switch (extractor_switch) takes the outputs from the models and stores the line emissions in dictionaries, which are created and updated to function as a line database. A third switch (bestmodel_switch) allows for comparison with observations, identifying the available CLOUDY model that best fits the user-defined emission lines. Future updates will include the possibility of creating 1-D and 2-D posterior distributions to visualize the parameter space when comparing with observations. 
