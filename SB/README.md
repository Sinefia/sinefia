# **SINEFIA_SB**

SINEFIA_SB is a Python-based tool for automating the modeling of observational data using CLOUDY simulations and SED templates. It facilitates the generation of parameter grids, execution of CLOUDY models and statistical comparison of results with observational data. The tool is specifically designed to analyze starburst galaxy spectra and emission line intensities, identifying the best-fit parameters for given observational datasets.

---

## **Features**

- **SED Template Integration**: Uses stellar templates (e.g., BPASS) for astrophysical modeling.
- **Grid Generation**: Creates parameter combinations for age, ionization parameter (U), hydrogen density (nh), metallicity (Z), covering factor (cov), and column density (NH).
- **CLOUDY Automation**: Prepares input files, runs simulations, and manages outputs.
- **Best Model Identification**: Computes chi-squared values to identify the best-fit parameters for each observational dataset.

---

## **File Organization**

The repository is organized as follows:

- `sinefia_SB.py`: Main script for running the SINEFIA pipeline.
- `galaxy_data_loader.py`: Script to load and parse observational data.
- `linelist.dat`: List of emission lines used in CLOUDY simulations.
- `BPASS/`: Directory containing SED templates for CLOUDY.
  - `BPASSv2_imf135_300_burst_single.ascii`: Example BPASS stellar template file.
- `Spitzer_galaxies_for_CLOUDY_comparison/`: Folder for observational data.
  - Each galaxy has its own folder containing:
    - `<GalaxyName_combined.txt>`: Input file with emission line data.
    - `<GalaxyName_best_model.txt>`: Best-fit result file.
- `models/`: Directory where CLOUDY outputs and JSON results are stored.
  - This folder will contain generated JSON files for model outputs.
- `README.md`: Documentation file explaining the project.

---


### **Comments on File Organization**

#### `sinefia_SB.py`
The main script that drives the entire SINEFIA workflow.  
Handles parameter grid generation, CLOUDY execution, data extraction, and chi-squared analysis.

#### `galaxy_data_loader.py`
Responsible for reading and parsing observational data files.  
Designed to be easily customizable for different data formats.

#### `linelist.dat`
Specifies the emission lines of interest for CLOUDY simulations.

#### `BPASS/`
Contains SED templates required by CLOUDY.  
Example file: `BPASSv2_imf135_300_burst_single.ascii`.

#### `Spitzer_galaxies_for_CLOUDY_comparison/`
Observational data organized by galaxy.  
Each galaxy folder contains:
- A `*_combined.txt` file with emission line data.
- A `*_best_model.txt` file with best-fit results.

#### `models/`
Directory for storing CLOUDY results and generated JSON files.

