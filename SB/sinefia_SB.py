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
import glob
import itertools
import json
from scipy.stats import chi2
import time
from galaxy_data_loader import load_galaxy_data  # Import the external loader function 



#######################################################
#######################################################


Ncpus = 56 # Number of CPUs used for the parallelization in grids
cloudy_executable = "/Users/cristina/Desktop/Work/cloudy/c22.02/source/cloudy.exe" # Location of the cloudy.exe file



# Set the paths relative to this script location (no matter where are located, it must be in the same place of the script)
base_path = os.path.dirname(os.path.realpath(__file__))  # Base path of the script 
galaxy_folder = os.path.join(base_path, "Spitzer_galaxies_for_CLOUDY_comparison")

list_lines_file = os.path.join(base_path, "linelist.dat")


##### Observational values to compare with #####

# Load the list of galaxy names
galaxies = np.genfromtxt('/Users/cristina/Desktop/Work/SPITZER/Scripts/example.dat', usecols=[2], dtype=str)

# Load the observational data
galaxies_data, obs_files = load_galaxy_data(galaxy_folder, galaxies, return_obs_files=True)

print("Galaxy data structure populated. Ready for comparison.")




# Example of processing and using the data
for file_key, obs_data in galaxies_data.items():
    print(f"Processing file: {file_key}")

    obs_linenames = [entry[0] for entry in obs_data]
    obs_linevalues = [entry[1] for entry in obs_data]
    obs_lineerrors = [entry[2] for entry in obs_data]

    print("Line names:", obs_linenames)
    print("Line values:", obs_linevalues)
    print("Line errors:", obs_lineerrors)
    print('--------------------------')



# Function to convert numpy types to Python types
def convert_numpy_types(data):
    if isinstance(data, np.ndarray):
        return data.tolist()  # Convert numpy arrays to lists
    elif isinstance(data, np.int64):
        return int(data)  # Convert numpy int64 to Python int
    elif isinstance(data, np.float64):
        return float(data)  # Convert numpy float64 to Python float
    elif isinstance(data, dict):
        return {key: convert_numpy_types(value) for key, value in data.items()}
    elif isinstance(data, list):
        return [convert_numpy_types(item) for item in data]
    else:
        return data
    
    
    

### Input SB variables ###

# For non-grid parameters set step to 0

age = [1000000, 10000000, 19952623.1496888]  # Corresponding to 1 Myr, 10 Myr, 20 Myr according to the stellar template considered


U_init = -4.5 # Ionization parameter, start of the grid (in log)
U_end = -3.5 # Ionization parameter, end of the grid
U_step = 0.5 # Ionization parameter, step of the grid

nh_init = 2 # Hydrogen density at the illuminated face of the cloud, start of the grid (in log)
nh_end = 4 # Hydrogen density at the illuminated face of the cloud, end of the grid
nh_step = 1.0 # Hydrogen density, step of the grid


Z_init = -1.69897000433602 # to vary based on the stellar template considered
Z_end = -1.69897000433602 
Z_step = 0

covfac_init = 0.3 # Covering factor
covfac_end = 0.3
covfac_step = 0


# Column density, the stopping criteria (in log)
Nh_init = 21.0 
Nh_end = 21.0
Nh_step = 0


# Relative abundances for C, N and He, as in Nicholls+17 and Decarli+23
oh = 3.19*10**(-4) # Default O abundance in ISM
#c_abundance = np.log10(10**(-0.8)+10**(np.log10(oh*Z_init)+2.72))+np.log10(oh*Z_init)
#n_abundance = np.log10(10**(-1.732)+10**(np.log10(oh*Z_init)+2.19))+np.log10(oh*Z_init)
#he_abundance = -1.0783+np.log10(1+0.1703*Z_init)


 
    
 
### SWITCHES ###

cloudy_switch=0 # runs CLOUDY

extractor_switch=1 # takes from each folder the lines (.lineoutput), the grid of params (grid_U_hden.txt), and merges.

bestmodel_switch=1 # takes the parameters written above and estimate likelihoods, comparing with observations



#######################################################
#######################################################



output_folder = os.path.join(base_path, 'models') 

# Ensure the output folder exists
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    
# Path to your PopStar SED file
ascii_file = "/Users/cristina/Desktop/Work/cloudy/c22.02/BPASS/BPASSv2_imf135_300_burst_single.ascii"  # Location of the stellar template ascii file
 
    
# Arrays
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
    
#print(age,U,nh,Z,cov,NH)



# Create all parameter combinations
combinations = list(itertools.product(age, U, nh, Z, cov, NH))

# Read the linelist file to get the desired lines
with open(list_lines_file, 'r') as linelist_file:
    desired_lines = [line.strip().replace('  ', '_') for line in linelist_file if not line.startswith('#')]
    #print(desired_lines)
    
 


#######################################################
#######################################################




if cloudy_switch == 1:

    
    # Time estimation
    total_combinations = len(combinations)
    average_time_per_model = 20  # Let's assume 20 minutes as an average per model
    estimated_total_time = total_combinations * average_time_per_model

    print('_______________________')
    print(f'Estimated time to complete the grid: {estimated_total_time} minutes (~{estimated_total_time // 60} hours)')
    print(f'Running {total_combinations} models with {Ncpus} CPUs...')
    print('_______________________')

    start_time = time.time()  # Record the start time


    # Iterate over each combination of parameters
    for i, (age_value, U_value, nh_value, Z_value, cov_value, NH_value) in enumerate(combinations):
     # Create the folder to run cloudy. If it exists it will automatically overwrite the files on it.
        folder_name = f"SBmodel_age{age_value:.2f}_U{U_init}to{U_end}inc{U_step}_nh{nh_init}to{nh_end}inc{nh_step}_Z{Z_init}_covfac{covfac_init}_NH{Nh_init}_ctedens" 
        full_folder_path = os.path.join(base_path, folder_name) 

        #if not os.path.exists(folder_name):
            #os.makedirs(folder_name)
            
        # Ensure the folder exists before creating the input file
        if not os.path.exists(full_folder_path):
            os.makedirs(full_folder_path)
            
        # Generate the CLOUDY input file path 
        #cloudy_input_file = os.path.join(folder_name, "SBmodel.in")
        cloudy_input_file = os.path.join(full_folder_path, "SBmodel.in") #NEW

            
        # Check for existing models to avoid re-running if necessary
        for item in desired_lines:
            file_name = os.path.join(output_folder, item + '.json')

            if os.path.exists(file_name):
                with open(file_name, 'r') as json_file:
                    existing_data = json.load(json_file)

                    # Use .get() to avoid KeyError if the key doesn't exist
                    age_data = existing_data.get('age', [])
                    U_data = existing_data.get('U', [])
                    nh_data = existing_data.get('nh', [])
                    Z_data = existing_data.get('Z', [])
                    cov_data = existing_data.get('cov', [])
                    NH_data = existing_data.get('NH', [])

                    # Create a list of existing combinations if all keys are found
                    if age_data and U_data and nh_data and Z_data and cov_data and NH_data:
                        existing_combinations = list(zip(age_data, U_data, nh_data, Z_data, cov_data, NH_data))

                        existing_matches = [combination for combination in combinations if combination in existing_combinations]

                        if existing_matches:
                            print('_______________________')
                            print("Some parameter combination (age, log(U), log(nh), Z, cov, Nh) already exist in the database:")
                            print(existing_matches)
                            print('_______________________')
                            print('You may want to change the parameters to avoid running those models again and save time.')

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
            
 
        # Write the CLOUDY input commands into SBmodel.in
        with open(cloudy_input_file, "w") as file:
            file.write(f"# title SB model\n")
            file.write(f"# SB -> U, NH, nH, Zgas, Zstar, Age_star\n")
            file.write(f"table star \"../BPASS/BPASSv2_imf135_300_burst_single.ascii\" age = {age_value:.5f} years Z = {Z_value}\n")  # Use the correct PopStar syntax
            file.write(f"ionization parameter {U_value} vary\n")
            file.write(f"grid range from {U_init} to {U_end} with {U_step} dex steps ncpus {Ncpus}\n")
            file.write(f"cosmic rays background\n")
            file.write(f"hden {nh_value} vary\n")
            file.write(f"grid range from {nh_init} to {nh_end} with {nh_step} dex steps ncpus {Ncpus}\n")
            file.write(f"abundance ISM no grains\n")
            file.write(f"grains ISM\n")
            file.write(f"grains PAH\n")
            file.write(f"metals and grains {Z_value}\n")
            #file.write(f"element abundance carbon {c_abundance}\n")
            #file.write(f"element abundance nitrogen {n_abundance}\n")
            #file.write(f"element abundance helium {he_abundance}\n")
            file.write(f"constant pressure\n")
            file.write(f"covering factor {cov_value}\n")
            file.write(f"stop H2 column density {NH_value}\n")
            file.write(f"save last line list emergent absolute \"lineoutput.txt\" \"linelist.dat\"\n")

        # Copy necessary files and run CLOUDY
        #shutil.copy(cloudy_input_file, folder_name)
        shutil.copy(list_lines_file, folder_name)
        shutil.copy(ascii_file, folder_name)  # Copy the ASCII SED file into the folder 
        
        #os.chdir(folder_name)
        
        # Change the working directory to the model folder #NEW
        current_working_directory = os.getcwd()
        os.chdir(full_folder_path)  # Change to the folder where the model input file is

        
    
        # Print progress
        print(f"Running model {i + 1} of {total_combinations}...")
        
        # Read the CLOUDY input file
        with open(cloudy_input_file, "r") as file:
            init_lines = file.readlines()
            
        modified_lines = [line.replace("{age}", str(age)).replace("{U_init}", str(U_init)).replace("{U_end}", str(U_end)).replace("{U_step}", str(U_step)).replace("{nh_init}", str(nh_init)).replace("{nh_end}", str(nh_end)).replace("{nh_step}", str(nh_step)).replace("{Z_init}", str(Z_init)).replace("{covfac_init}", str(covfac_init)).replace("{Nh_init}", str(Nh_init)) for line in init_lines] #Simple version without grid 

        # Modify and write the lines in the CLOUDY input file with the dictionary variables
        with open(cloudy_input_file, "w") as file:
            file.writelines(modified_lines)

        # Run CLOUDY using the generated input file
        command = [cloudy_executable, "SBmodel.in"]
        subprocess.run(command)

        #os.chdir("..") 
        os.chdir(current_working_directory)
    
    end_time = time.time()
    elapsed_time = (end_time - start_time) / 60  # Convert to minutes

    # Final time report
    print(f'_______________________')
    print(f'Total time taken: {elapsed_time:.2f} minutes (~{elapsed_time // 60} hours)')
    print('_______________________')




#######################################################
#######################################################




if extractor_switch==1:
    
    
    for age_value in age: 
        # Create the output file in every folder for each combination of parameters
        folder_name = f"SBmodel_age{age_value:.2f}_U{U_init}to{U_end}inc{U_step}_nh{nh_init}to{nh_end}inc{nh_step}_Z{Z_init}_covfac{covfac_init}_NH{Nh_init}_ctedens"
        
        
        lineoutput_path = os.path.join(folder_name,  'SBmodellineoutput.txt')        
        try:
            with open(lineoutput_path, 'r') as lineoutput:
                lines = []    
    
                def parse_line_output(file):
                    for line in file:
                        if line.startswith('iteration'):
                            parts = line.split()
                            lines.append(parts[2:])
                            print(f"Intensities extracted: {parts[1:]}")  

                parse_line_output(lineoutput)
                
                
                
                # Add the print statement to ensure that all lines have been read
                print(f"Lines read from file: {len(lines)}")  # Ensure this prints all lines

        except FileNotFoundError:
            print('_______________________')
            print(f"File not found: {lineoutput_path}. Run the CLOUDY models for that parameter combination.")
            sys.exit(1)
            
        # Filter combinations for the current age
        age_combinations = [combo for combo in combinations if combo[0] == age_value]


        # Check how many combinations and lines were extracted for this age
        print(f"Lines extracted for age {age_value}: {len(lines)}")
        print(f"Combinations for age {age_value}: {len(age_combinations)}")

        
        

        # Loop for each dessired line      
        for i, item in zip(range(len(desired_lines)), desired_lines):
            item = item.replace(' ', '_').replace('/', '_')  # Replace spaces and slashes
            file_name = os.path.join(output_folder, item+'.json')

            # Initialize lists to store parameter values and intensity values for each line
            param_names = ['age', 'U', 'nh', 'Z', 'cov', 'NH', 'intensity']
            param_values = {param: [] for param in param_names}       
    
            
            
            # Assign values to keys
            for j, (a, b, c, d, e, f) in enumerate(age_combinations): # Assign intensity values by matching each combination
                
                
                param_values['age'].append(a)
                param_values['U'].append(b)
                param_values['nh'].append(c)
                param_values['Z'].append(d)
                param_values['cov'].append(e)
                param_values['NH'].append(f)
                
                # Now assign the corresponding intensity for each line from `lines`
                if j < len(lines) and i < len(lines[j]):  # Ensure index j is within bounds
                    param_values['intensity'].append(float(lines[j][i]))  # Assign intensity for the current line (i)
                else:
                    print(f"Warning: No intensity found for combination {j}, age {age_value}")

            # Debugging: Check what has been computed so far
            print(f"Computed models for {item} (age {age_value}) are:")
            print(param_values)
            
            
            print('_______________________')
            print(f'Computed models for {item} are:')
            print(param_values)
            

            # Check if the file already exists, and if so load the existing dictionary
            if os.path.exists(file_name):
                with open(file_name, 'r') as json_file:
                    existing_data = json.load(json_file)
            else:
                existing_data = {'age':[],'U':[],'nh':[],'Z':[],'cov':[],'NH':[],'intensity':[]}

            print('_______________________')
        
            # Convert any numpy types to regular Python types
            param_values = convert_numpy_types(param_values)

            # Define the function to merge with existing data
            def merge(existing_data, param_values):
                
                # Check if the combination already exists in existing data
                existing_combinations = list(zip(
                    existing_data['age'],
                    existing_data['U'],
                    existing_data['nh'],
                    existing_data['Z'],
                    existing_data['cov'],
                    existing_data['NH']
                ))
            
                # Append new combinations and intensities
                for j, combination in enumerate(zip(
                    param_values['age'],
                    param_values['U'],
                    param_values['nh'],
                    param_values['Z'],
                    param_values['cov'],
                    param_values['NH']
                )):
                    if combination not in existing_combinations:
                        for key in param_values:
                            if len(param_values[key]) > j:
                                existing_data[key].append(param_values[key][j])
                            else:
                                print(f"Skipping {key} for combination {combination} due to missing values")


           # Update existing data with new
            merge(existing_data, param_values)   

            # Convert existing_data to JSON-serializable types
            existing_data = convert_numpy_types(existing_data)

            with open(file_name, 'w') as json_file:
                json.dump(existing_data, json_file, indent=4)

            print('_______________________')
            print(f'Data saved to {file_name}')
            print('_______________________')
 



#######################################################
#######################################################




if bestmodel_switch == 1:
    # Iterate over each observational file to compute chi-squared independently
    for obs_file in obs_files:
        print(f"Processing observational file: {obs_file}")
        
        # Read the observational file
        obs_linenames = []
        obs_linevalues = []
        obs_lineerrors = []

        with open(obs_file, 'r') as obslines_file:
            for line in obslines_file:
                if line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) == 3:  # Ensure the line has the required 3 columns
                    obs_linenames.append(parts[0])
                    obs_linevalues.append(float(parts[1]))
                    obs_lineerrors.append(float(parts[2]))

        # Skip files with no valid observational data
        if not obs_linenames:
            print(f"Warning: Observational file {obs_file} contains no valid data. Skipping.")
            continue

        # Convert to numpy arrays for calculations
        obs_linevalues = np.array(obs_linevalues)
        obs_lineerrors = np.array(obs_lineerrors)

        # Load the model data for the lines in this file
        model_data = {}
        for item in obs_linenames:
            file_name = os.path.join(output_folder, f"{item}.json")
            if os.path.exists(file_name):
                with open(file_name, 'r') as json_file:
                    model_data[item] = json.load(json_file)
            else:
                print(f"The line {item} does not exist in the database")

        # Skip files with no matching model data
        if not model_data:
            print(f"Warning: No model data found for observational file {obs_file}. Skipping.")
            continue

        # Prepare a flux matrix for this file
        n_lines = len(obs_linenames)
        n_combinations = len(combinations)
        model_fluxes = np.zeros((n_combinations, n_lines))

        # Populate model_fluxes
        for line_idx, line_name in enumerate(obs_linenames):
            if line_name in model_data:
                line_data = model_data[line_name]
                for comb_idx in range(len(combinations)):
                    if comb_idx < len(line_data['intensity']):
                        model_fluxes[comb_idx, line_idx] = line_data['intensity'][comb_idx]
                    else:
                        print(f"Warning: Missing intensity for line {line_name}, combination {comb_idx}")

        # Perform scaling and chi-squared calculations for this file
        numerator = np.nansum(
            obs_linevalues * model_fluxes / (obs_lineerrors**2), axis=1
        )
        denominator = np.nansum(
            model_fluxes**2 / (obs_lineerrors**2), axis=1
        )
        # Handle invalid scaling gracefully
        with np.errstate(invalid='ignore', divide='ignore'):
            scal = numerator / denominator

        # Avoid division errors by ensuring valid scal values
        if np.all(np.isnan(scal)):
            print(f"Warning: Unable to compute scaling factors for file {obs_file}. Skipping.")
            continue

        residuals = obs_linevalues - scal[:, np.newaxis] * model_fluxes
        chi_sq = np.nansum((residuals**2) / (obs_lineerrors**2), axis=1)

        # Find the model with the smallest chi-squared value for this file
        min_chi2_idx = np.argmin(chi_sq)
        min_chi2_value = chi_sq[min_chi2_idx]

        # Extract corresponding parameters for the best model for this file
        b_age = model_data[obs_linenames[0]]["age"][min_chi2_idx]
        b_U = model_data[obs_linenames[0]]["U"][min_chi2_idx]
        b_nh = model_data[obs_linenames[0]]["nh"][min_chi2_idx]
        b_Z = model_data[obs_linenames[0]]["Z"][min_chi2_idx]
        b_cov = model_data[obs_linenames[0]]["cov"][min_chi2_idx]
        b_NH = model_data[obs_linenames[0]]["NH"][min_chi2_idx]

        # Print the results for this file
        print(f"File: {obs_file}")
        print(f"Smallest chi-squared value: {min_chi2_value}")
        print(f"Best parameters (age, U, nh, Z, cov, NH): ({b_age}, {b_U}, {b_nh}, {b_Z}, {b_cov}, {b_NH})")
        print('--------------------------')

        # Save the results to a text file
        output_file = obs_file.replace("_combined.txt", "_best_model.txt")
        with open(output_file, 'w') as f:
            f.write(f"chi-squared   {min_chi2_value}\n")
            f.write(f"age           {b_age}\n")
            f.write(f"U             {b_U}\n")
            f.write(f"nh            {b_nh}\n")
            f.write(f"Z             {b_Z}\n")
            f.write(f"cov           {b_cov}\n")
            f.write(f"NH            {b_NH}\n")

        print(f"Results saved to: {output_file}")



#######################################################
#######################################################
   
    
