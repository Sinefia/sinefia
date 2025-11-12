#######################################################
#######################         #######################
####################### SINEFIA #######################
#######################         #######################
#######################################################


#######################################################
#######################################################


import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import subprocess
import shutil
import itertools
import json
from scipy.stats import chi2


#######################################################
#######################################################
#######################################################

Ncpus = 8 # Number of CPUs used for the parallelization in grids
cloudy_executable = "/Users/roman/Documents/PhD/cloudy/c17.03/source/cloudy.exe" # Location of the cloudy.exe file


### Input AGN variables ###

# For non-grid parameters set step to 0

aox_init = -1.4 # X-ray to UV-optical ratio
aox_end = -1.4
aox_step = 0 #

U_init = -3.5 # Ionization parameter, start of the grid (in log)
U_end = -2.0 # Ionization parameter, end of the grid
U_step = 0.5 # Ionization parameter, step of the grid

### Input cloud variables ###

density_law = 'constant density' #diff model, not param

nh_init = 1.0 # Hydrogen density at the iluminated face of the cloud, start of the grid (in log)
nh_end = 6.0 # Hydrogen density at the iluminated face of the cloud, end of the grid
nh_step = 1.0 # Hydrogen density, step of the grid


Z_init = 1. # Solar metallicity of both metals and grains. Abundance also should vary for low metallicities! And will be as density_law, not a param
#O/H
Z_end = 1. 
Z_step = 0

covfac_init = 0.3 # Covering factor
covfac_end = 0.3
covfac_step = 0

Nh_init = 23.0 # Column density, the stopping criteria (in log)
Nh_end = 23.0
Nh_step = 0

### Observational values to compare with ###

obs_lines_file = "ObsLineList.txt"

#Things missing:
    # wrapper for other params (aox, Z, cov fac, Nh)
    # cornerplots

### SWITCHES ###

cloudy_switch=0 # runs CLOUDY

extractor_switch=0 # takes from each folder the lines (.lineoutput), the grid of params (grid_U_hden.txt), and merges.

bestmodel_switch=1 # takes the parameters written above and estimate likelihoods, comparing with observations

#cornerplot_switch=0 #takes the likelihoods and makes a corner plot, COMBINE

#######################################################
#######################################################
#######################################################

output_folder = 'models'

# Ensure the output folder exists
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    
# Create the folder to run cloudy. If it exists it will automatically overwrite the files on it.
folder_name = f"AGNmodel_aox{aox_init}_U{U_init}to{U_end}inc{U_step}_nH{nh_init}to{nh_end}inc{nh_step}_Z{Z_init}_covfac{covfac_init}_NH{Nh_init}_ctedens" 
if not os.path.exists(folder_name):
    os.makedirs(folder_name)
    
# Arrays
if aox_step==0:
    aox = np.arange(aox_init, aox_init+1, 1)
else:
    aox = np.arange(aox_init, aox_end+aox_step, aox_step)

if U_step==0:
    U = np.arange(U_init, U_init+1, 1)
else:
    U = np.arange(U_init, U_end+U_step, U_step)
    
if nh_step==0:
    nh = np.arange(nh_init, nh_init+1, 1)
else:
    nh = np.arange(nh_init, nh_end+nh_step, nh_step)
    
if Z_step==0:
    Z = np.arange(Z_init, Z_init+1, 1)
else:
    Z = np.arange(Z_init, Z_end+Z_step, Z_step)

if covfac_step==0:
    cov = np.arange(covfac_init, covfac_init+1, 1)
else:
    cov = np.arange(covfac_init, covfac_end+covfac_step, covfac_step)

if Nh_step==0:
    NH = np.arange(Nh_init, Nh_init+1, 1)
else:
    NH = np.arange(Nh_init, Nh_end+Nh_step, Nh_step)
    
#print(aox,U,nh,Z,cov,NH)

combinations = list(itertools.product(aox, U, nh, Z, cov, NH))
#print(combinations)

list_lines_file = "LineList.dat"

# Read the linelist file to get the desired lines
with open(list_lines_file, 'r') as linelist_file:
    desired_lines = [line.strip().replace(' ', '_') for line in linelist_file if not line.startswith('#')]
    #print(desired_lines)
    
# Relative abundances for C, N and He, as in Nicholls+17 and Decarli+23
oh = 3.19*10**(-4) # Default O abundance in ISM
c_abundance = np.log10(10**(-0.8)+10**(np.log10(oh*Z_init)+2.72))+np.log10(oh*Z_init)
n_abundance = np.log10(10**(-1.732)+10**(np.log10(oh*Z_init)+2.19))+np.log10(oh*Z_init)
he_abundance = -1.0783+np.log10(1+0.1703*Z_init)

#######################################################
    
if cloudy_switch==1: 

    for item in desired_lines:
        file_name = os.path.join(output_folder, item+'.json')
        if os.path.exists(file_name):
            with open(file_name, 'r') as json_file:
                existing_data = json.load(json_file)
                existing_combinations = list(zip(existing_data['aox'], existing_data['U'], existing_data['nh'], existing_data['Z'], existing_data['cov'], existing_data['NH']))
            
                existing_matches = [combination for combination in combinations if combination in existing_combinations]
                if existing_matches:
                    print('_______________________') 
                    print("Some parameter combination (aox, log(U), log(nh), Z, cov, Nh) already exist in the database:")
                    print(existing_matches)
                    print('_______________________') 
                    print('You may want to change the parameters to avoid running that models again and save time.') 
                    while True:
                        user_input = input("Do you still want to run your input models? [y/n]: ").strip().lower()
                        if user_input == "y":
                            break
                        elif user_input == "n":
                            sys.exit(1)
                        else:
                            print("Invalid input. Please enter 'y' or 'n'")       
        else:
            pass
            
    print('_______________________') 
    print('As a reference, the grid of models will take approximately ',len(combinations)*4,' minutes to run with 8 CPUs')
    print('_______________________') 

    cloudy_input_file = "agnmodel.in"

    # Copy CLOUDY input files to new folder and change to it
    shutil.copy(cloudy_input_file, folder_name)
    shutil.copy(list_lines_file, folder_name)
    os.chdir(folder_name)
    
    # Read the CLOUDY input file
    with open(cloudy_input_file, "r") as file:
        init_lines = file.readlines()
        
    modified_lines = [line.replace("{aox_init}", str(aox_init)).replace("{U_init}", str(U_init)).replace("{U_end}", str(U_end)).replace("{U_step}", str(U_step)).replace("{nh_init}", str(nh_init)).replace("{nh_end}", str(nh_end)).replace("{nh_step}", str(nh_step)).replace("{Z_init}", str(Z_init)).replace("{covfac_init}", str(covfac_init)).replace("{Nh_init}", str(Nh_init)).replace("{c_abundance}", str(c_abundance)).replace("{n_abundance}", str(n_abundance)).replace("{he_abundance}", str(he_abundance)) for line in init_lines] #Simple version without grid 

    # Modify and write the lines in the CLOUDY input file with the dictionary variables
    with open(cloudy_input_file, "w") as file:
        file.writelines(modified_lines)


    # Run CLOUDY in the folder
    command = [cloudy_executable, cloudy_input_file]
    subprocess.run(command)

    os.chdir("..")


#######################################################
#######################################################

if extractor_switch==1:
    
    lineoutput_path = os.path.join(folder_name,  'agnmodellineoutput.txt')        
    try:
        lineoutput = open(lineoutput_path, 'r')
        #lineoutput.close()
    except FileNotFoundError:
        print('_______________________')
        print(f"File not found: {lineoutput_path}. Run the CLOUDY models for that parameter combination.")
        sys.exit(1)
        
    lines = []
    def it3_finder(file):
        for line in file:
            if 'iteration' in line:
                lines.append(line.split()[2:])
    
    it3_finder(lineoutput)

    # Loop for each dessired line      
    for i, item in zip(range(len(desired_lines)), desired_lines):
        file_name = os.path.join(output_folder, item+'.json')

        # Initialize lists to store parameter values and intensity values for each line
        param_names = ['aox', 'U', 'nh', 'Z', 'cov', 'NH', 'intensity']
        param_values = {param: [] for param in param_names}
        
        # Assign values to keys
        for j, (a, b, c, d, e, f) in enumerate(combinations):
            param_values['aox'].append(a)
            param_values['U'].append(b)
            param_values['nh'].append(c)
            param_values['Z'].append(d)
            param_values['cov'].append(e)
            param_values['NH'].append(f)
                
        for vals in lines:
            param_values['intensity'].append(float(vals[i]))
        
        print('_______________________')
        print(f'Computed models for {item} are:')
        print(param_values)
    
        # Check if the file already exists, and if so load the existing dictionary
        if os.path.exists(file_name):
            with open(file_name, 'r') as json_file:
                existing_data = json.load(json_file)
        else:
            existing_data = {'aox':[],'U':[],'nh':[],'Z':[],'cov':[],'NH':[],'intensity':[]}

        print('_______________________')
        
        # Define the function to merge with existing data
        def merge(existing_data, param_values):
            # Check if the combination already exists in existing data
            existing_combinations = list(zip(
                existing_data['aox'],
                existing_data['U'],
                existing_data['nh'],
                existing_data['Z'],
                existing_data['cov'],
                existing_data['NH']
            ))
            
            existing_matches = [combination for combination in combinations if combination in existing_combinations]
            #print(existing_matches)
            nonexisting_matches = [combination for combination in combinations if combination not in existing_combinations]
            #print(nonexisting_matches)
            
            if nonexisting_matches:      
                for combination in nonexisting_matches:
                    print(f"New (aox, log(U), log(nh), Z, cov, Nh) combination {combination} added to existing data")
                for key in param_values:
                    existing_data[key].extend(param_values[key])                
            if existing_matches:
                for combination in existing_matches:
                    print(f"(aox, log(U), log(nh), Z, cov, Nh) combination {combination} already exists")
                
            #for combination in combinations:
            #    if combination not in existing_combinations:
                    
        # Update existing data with new
        merge(existing_data, param_values)    

        with open(file_name, 'w') as json_file:
            json.dump(existing_data, json_file, indent=4)

        print('_______________________')
        print(f'Data saved to {file_name}')
        print('_______________________')
    
   # os.system('rm -rf '+folder_name)
   
   
#######################################################
#######################################################

if bestmodel_switch==1: 

    ### Read observations file   
    obs_linenames = []
    obs_linevalues = []
    obs_lineerrors = []
    
    with open(obs_lines_file, 'r') as obslines_file:
        for line in obslines_file:
            if line.startswith('#'):
                continue
            parts = line.split()
            obs_linenames.append(parts[0])
            obs_linevalues.append(float(parts[1]))
            obs_lineerrors.append(float(parts[2]))
    #print(obs_linenames,' ',obs_linevalues,' ',obs_lineerrors)

    # Load existing dictionaries of desired lines
    for item in obs_linenames:
        file_name = os.path.join(output_folder, item+'.json')
        if os.path.exists(file_name):
            with open(file_name, 'r') as json_file:
                globals()[item] = json.load(json_file) # Name the file like the line
        else:  
            print(f'The line {item} does not exist in the database')
            
     # Create a new dictionary to store the lines
    model_data = {}
        
    for item in obs_linenames:
        current_data = globals()[item]
        model_data[item] = {}
            
        for key in ["aox", "U", "nh", "Z", "cov", "NH"]:
            model_data[item][key] = current_data[key]
            model_data[item]["intensity"] = current_data["intensity"]
    
    # And now store only the intensities
    n_lines = len(obs_linenames)
    model_fluxes = np.zeros((len(combinations), n_lines))
    for line_idx, (line_name, line_data) in enumerate(model_data.items()):  # Iterate over lines
        for comb_idx in range(len(combinations)):  # Iterate over parameter combinations
            model_fluxes[comb_idx, line_idx] = line_data['intensity'][comb_idx]
            
    
    # Scaling and chi2 calculation
    numerator = np.nansum(obs_linevalues*model_fluxes/(np.array((obs_lineerrors))**2), axis=1)
    denominator = np.nansum(model_fluxes**2/np.array((obs_lineerrors))**2 , axis=1)
    scal = numerator/denominator #1 scaling per model
    
    residuals = obs_linevalues - scal[:, np.newaxis]*model_fluxes 
    #print('residuals',residuals)
    chi_sq = np.nansum((residuals**2) / np.array((obs_lineerrors))**2, axis=1)
    
    # Find the index of the smallest chi-squared value
    min_chi2_idx = np.argmin(chi_sq)
    min_chi2_value = chi_sq[min_chi2_idx]

    # Print corresponding aox and U for the model with the smallest chi-squared
    b_aox = model_data[list(model_data.keys())[0]]["aox"][min_chi2_idx]
    b_U = model_data[list(model_data.keys())[0]]["U"][min_chi2_idx]
    b_nh = model_data[list(model_data.keys())[0]]["nh"][min_chi2_idx]
    b_Z = model_data[list(model_data.keys())[0]]["Z"][min_chi2_idx]
    b_cov = model_data[list(model_data.keys())[0]]["cov"][min_chi2_idx]
    b_NH = model_data[list(model_data.keys())[0]]["NH"][min_chi2_idx]

    print(f"Smallest chi-squared value: {min_chi2_value}")
    print(f"Corresponding best parameters (aox, U, nh, Z, cov, NH): ({b_aox}, {b_U}, {b_nh}, {b_Z}, {b_cov}, {b_NH})")

    print('_______________________')
    
#######################################################
#######################################################

    
    
    

    
    
