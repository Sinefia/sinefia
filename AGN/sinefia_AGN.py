###########################################################
#######################             #######################
####################### SINEFIA AGN #######################
#######################             #######################
###########################################################

import numpy as np
import matplotlib.pyplot as plt
#from scipy.stats import chi2
from multiprocessing import Pool
from functools import partial
import os, sys, shutil, subprocess, json, itertools, traceback, csv, re, glob, math

#######################################################
# User parameters (edit as needed)
#######################################################

Ncpus = 8  # Number of CPUs used for the parallelization
cloudy_executable = "/Users/roman/Documents/PhD/cloudy/c17.03/source/cloudy.exe"  # Path to cloudy executable

# For non-grid parameters set step to 0
aox_init = -1.4 # X-ray to UV-optical ratio
aox_end = -1.4
aox_step = 0

U_init = -3.5 # Ionization parameter (log)
U_end = -3.0
U_step = 0.5

density_law = 'constant density'

nH_init = 1.0 # Hydrogen density at the iluminated face of the cloud
nH_end = 1.0
nH_step = 0.0

Z_init = 1.0 # Solar metallicity of both metals and grain
Z_end = 1.0
Z_step = 0

covfac_init = 0.3 # Covering factor
covfac_end = 0.6
covfac_step = 0.3

NH_init = 23.0 # Column density (log), the stopping criteria
NH_end = 23.0
NH_step = 0

obs_lines_file = "ObsLineList.txt"
list_lines_file = "LineList.dat"
cloudy_input_file = "agnmodel.in"

# Switch to run CLOUDY
cloudy_switch = 0
# Switch to extract the model results and sotre them
extractor_switch = 0
# Find best model from database comparing with observations
best_switch = 1

# If model-folder exists: overwrite or skip?
overwrite_existing = True

# Output folder for all model folders
output_folder = 'models'

#######################################################
# Build parameter arrays
#######################################################
if aox_step == 0:
    aox = np.arange(aox_init, aox_init + 1, 1)
else:
    aox = np.arange(aox_init, aox_end + aox_step, aox_step)

if U_step == 0:
    U = np.arange(U_init, U_init + 1, 1)
else:
    U = np.arange(U_init, U_end + U_step, U_step)

if nH_step == 0:
    nH = np.arange(nH_init, nH_init + 1, 1)
else:
    nH = np.arange(nH_init, nH_end + nH_step, nH_step)

if Z_step == 0:
    Z = np.arange(Z_init, Z_init + 1, 1)
else:
    Z = np.arange(Z_init, Z_end + Z_step, Z_step)

if covfac_step == 0:
    cov = np.arange(covfac_init, covfac_init + 1, 1)
else:
    cov = np.arange(covfac_init, covfac_end + covfac_step, covfac_step)

if NH_step == 0:
    NH = np.arange(NH_init, NH_init + 1, 1)
else:
    NH = np.arange(NH_init, NH_end + NH_step, NH_step)

# All combinations
combinations = list(itertools.product(aox, U, nH, Z, cov, NH))

#######################################################
# Read desired lines list (keeps same behaviour)
#######################################################
with open(list_lines_file, 'r') as linelist_file:
    desired_lines = [line.strip().replace(' ', '_') for line in linelist_file if not line.startswith('#')]

#######################################################
# Prepare output folder
#######################################################
def sanitize_val_for_name(val): # Convert negative signs -> 'm', decimal points -> 'p'
    s = f"{val:.6g}"  # compact representation
    s = s.replace('-', 'm').replace('.', 'p').replace('+', '')
    return s

def make_run_folder_name(aox, U, nH, Z, cov, NH):
    return f"model_aox{sanitize_val_for_name(aox)}_U{sanitize_val_for_name(U)}_nH{sanitize_val_for_name(nH)}_Z{sanitize_val_for_name(Z)}_cov{sanitize_val_for_name(cov)}_NH{sanitize_val_for_name(NH)}"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

#######################################################
# Worker function for a single combination
#######################################################
# combo: tuple (aox, U, nH, Z, cov, NH)
# This function creates a run folder, copies the input file(s), modifies placeholders, writes params.json and executes CLOUDY with subprocess.run(..., cwd=run_folder)
def run_single_combination(combo, cloudy_exe, cloudy_input_template, list_lines_file_local, output_root, overwrite=overwrite_existing):
    aox_c, U_c, nH_c, Z_c, cov_c, NH_c = combo
    run_folder = os.path.join(output_root, make_run_folder_name(aox_c, U_c, nH_c, Z_c, cov_c, NH_c))

    # Skip if folder exists and overwrite is False
    if os.path.exists(run_folder) and not overwrite:
        return {"combo": combo, "status": "skipped", "folder": run_folder}


    # Create run folder (if exists and overwrite=True, we clear it)
    if os.path.exists(run_folder) and overwrite:
        shutil.rmtree(run_folder)
    os.makedirs(run_folder, exist_ok=True)

    # Copy input files into run folder
    shutil.copy(cloudy_input_template, run_folder)
    shutil.copy(list_lines_file_local, run_folder)

    # Compute abundances for this Z (C, N, and He)
    def compute_abundances(Z_current):
        oh = 3.19e-4
        c_abundance = np.log10(10**(-0.8) + 10**(np.log10(oh * Z_current) + 2.72)) + np.log10(oh * Z_current)
        n_abundance = np.log10(10**(-1.732) + 10**(np.log10(oh * Z_current) + 2.19)) + np.log10(oh * Z_current)
        he_abundance = -1.0783 + np.log10(1 + 0.1703 * Z_current)
        return c_abundance, n_abundance, he_abundance
    c_abundance, n_abundance, he_abundance = compute_abundances(Z_c)

    # Replace variables in input file agnmodel.in
    template_path = os.path.join(run_folder, os.path.basename(cloudy_input_template))
    with open(template_path, 'r') as f:
        init_lines = f.readlines()

    modified_lines = [
        line.replace("{aox_init}", str(aox_c))
            .replace("{U_init}", str(U_c))
            .replace("{nH_init}", str(nH_c))
            .replace("{Z_init}", str(Z_c))
            .replace("{covfac_init}", str(cov_c))
            .replace("{NH_init}", str(NH_c))
            .replace("{c_abundance}", str(c_abundance))
            .replace("{n_abundance}", str(n_abundance))
            .replace("{he_abundance}", str(he_abundance))
            .replace("{density_law}", str(density_law))
        for line in init_lines
    ]

    with open(template_path, 'w') as f:
        f.writelines(modified_lines)

    # Save params for each run
    params = {
        "aox": float(aox_c),
        "U": float(U_c),
        "nH": float(nH_c),
        "Z": float(Z_c),
        "cov": float(cov_c),
        "NH": float(NH_c),
        "c_abundance": float(c_abundance),
        "n_abundance": float(n_abundance),
        "he_abundance": float(he_abundance),
        "cloudy_input": os.path.basename(cloudy_input_template)
    }
    #with open(os.path.join(run_folder, "params.json"), "w") as jf:
    #    json.dump(params, jf, indent=2)
    txt_path = os.path.join(run_folder, "params.txt")
    with open(txt_path, "w") as f:
        f.write("# model parameters\n")
        for k, v in params.items():
            f.write(f"{k} = {v}\n")

    # Run CLOUDY inside the run folder
    subprocess.run([cloudy_exe, os.path.basename(cloudy_input_template)], cwd=run_folder)
    

#######################################################
# Extracting and storing models
#######################################################

def sanitize_line_name(name): # To store emission line names
    s = name.strip()
    s = re.sub(r'[^\w\s\-\.]', '', s)       # remove some characters
    s = re.sub(r'\s+', '_', s)              # replace space with underscore
    return s

def params_txt(path): # Read params.txt in models folders
    params = {}
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '=' in line:
                k, v = [x.strip() for x in line.split('=', 1)]
                try:
                    params[k] = float(v)
                except Exception:
                    params[k] = v
    return params
 
def agnmodellineoutput_txt(path): # Read model lines fluxes from agnmodellineoutput.txt
    with open(path, 'r') as f:
        lines = [ln.rstrip('\n') for ln in f if ln.strip()]

    # model lines header
    header_line = next((ln for ln in lines if ln.lstrip().startswith('#')), None)
    modellines = header_line.lstrip('#').split('\t')
    #if '\t' in header_line:
    #    modellines = header_line.lstrip('#').split('\t')
    #else:
    #    modellines = re.split(r'\s{2,}', header_line.lstrip('#'))

    header = [t.strip().replace(' ', '_') for t in modellines[1:] if t.strip()]
    # pick last iteration and then the flux values
    iter_lines = [ln for ln in lines if ln.lower().lstrip().startswith('iteration')]
    last_iter = iter_lines[-1].strip()

    parts = re.split(r'\t|\s+', last_iter)
    flux_values = parts[2:] 

    values = []
    for v in flux_values:
        try:
            values.append(float(v))
        except Exception:
            values.append(v)
    return header, values

def update_lines_database(output_root='models'):
    lines_db_dir = os.path.join(output_root, 'lines_db') # Lines database
    os.makedirs(lines_db_dir, exist_ok=True)

    # gather models
    model_entries = []
    for entry in sorted(os.listdir(output_root)):
        folder = os.path.join(output_root, entry)
        if not os.path.isdir(folder):
            continue
        params_path = os.path.join(folder, 'params.txt')
        lineoutput_path = os.path.join(folder, 'agnmodellineoutput.txt')
        if not (os.path.exists(params_path) and os.path.exists(lineoutput_path)):
            continue

        params = params_txt(params_path)
        header, values = agnmodellineoutput_txt(lineoutput_path)
        if not header or not values:
            continue

        n = min(len(header), len(values))
        line_map = { header[i]: values[i] for i in range(n) }

        # Parameters
        def get_param(pdict, keys, default=math.nan):
            for kk in keys:
                if kk in pdict:
                    return pdict[kk]
            return default

        aox_v = get_param(params, ['aox'])
        U_v   = get_param(params, ['U'])
        nH_v  = get_param(params, ['nH'])
        Z_v   = get_param(params, ['Z'])
        cov_v = get_param(params, ['cov'])
        NH_v  = get_param(params, ['NH'])  
        
        keytuple = (float(aox_v), float(U_v), float(nH_v), float(Z_v), float(cov_v), float(NH_v))
        
        model_entries.append({'folder': folder, 'params': params, 'keytuple': keytuple, 'line_map': line_map})

    # collect line names
    all_lines = set()
    for me in model_entries:
        all_lines.update(me['line_map'].keys())

    # update CSV per line
    for line_name in sorted(all_lines):
        safe_name = sanitize_line_name(line_name)
        csv_path = os.path.join(lines_db_dir, f"{safe_name}.csv")

        # load existing keys
        existing_keys = set()
        if os.path.exists(csv_path):
            with open(csv_path, 'r', newline='') as cf:
                reader = csv.DictReader(cf)
                for row in reader:
                    try:
                        k = (
                            float(row['aox']),
                            float(row['U']),
                            float(row['nH']),
                            float(row['Z']),
                            float(row['cov']),
                            float(row['NH']),
                        )
                        existing_keys.add(k)
                    except Exception:
                        pass

        write_header = not os.path.exists(csv_path)
        appended = 0
        with open(csv_path, 'a', newline='') as cf:
            fieldnames = ['aox', 'U', 'nH', 'Z', 'cov', 'NH', 'intensity']
            writer = csv.DictWriter(cf, fieldnames=fieldnames)
            if write_header:
                writer.writeheader()

            for me in model_entries:
                if line_name not in me['line_map']:
                    continue
                key = me['keytuple']
                if key is None:
                    continue
                if key in existing_keys:
                    continue
                aox_v, U_v, nH_v, Z_v, cov_v, NH_v = key
                intensity = me['line_map'][line_name]
                try:
                    intensity = float(intensity)
                except Exception:
                    intensity = math.nan
                writer.writerow({
                    'aox': f"{aox_v:.6g}",
                    'U': f"{U_v:.6g}",
                    'nH': f"{nH_v:.6g}",
                    'Z': f"{Z_v:.6g}",
                    'cov': f"{cov_v:.6g}",
                    'NH': f"{NH_v:.6g}",
                    'intensity': f"{intensity:.6g}"
                })
                appended += 1
        if appended:
            print(f"Appended {appended} rows to {csv_path}")

if extractor_switch == 1:
    update_lines_database(output_folder)
    
#######################################################
# Best model comparing with observations
#######################################################

def round_key_tuple(ktuple, ndigits=8): # Round to same number of digits
    try:
        return tuple([round(float(x), ndigits) for x in ktuple])
    except Exception:
        return None

def find_best_model(obs_file,
                    lines_db_dir='models/lines_db',
                    models_root='models',
                    param_order=('aox','U','nH','Z','cov','NH'),
                    rounding=8,
                    topn=5):
    # read obs
    obs_names = []; obs_fluxes = []; obs_errors = []
    with open(obs_file, 'r') as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln or ln.startswith('#'):
                continue
            parts = ln.split()
            if len(parts) < 3:
                raise ValueError(f"Bad line in obs file: {ln!r} (expected: line flux error)")
            obs_names.append(parts[0])
            obs_fluxes.append(float(parts[1]))
            obs_errors.append(float(parts[2]))
    if len(obs_names) == 0:
        print("No observed lines found in", obs_file)
        return None

    # observed lines in the database
    available = []
    missing = []
    line_csv_paths = {}
    for name in obs_names:
        fname = sanitize_line_name(name) + '.csv'
        p = os.path.join(lines_db_dir, fname)
        if os.path.exists(p):
            available.append(name)
            line_csv_paths[name] = p
        else:
            missing.append(name)

    if missing:
        print("Warning: the following observed lines are NOT in the database and will be ignored for this fit:")
        for m in missing:
            print("  -", m)
        print("Caution: line names need one _ for two letter names (e.g.: Ar), and two __ for single letter names (e.g.: S)")

    if not available:
        print("No observed lines are present in the database. Aborting.")
        return None

    # load available model lines
    line_maps = {}
    for name, path in line_csv_paths.items():
        mp = {}
        with open(path, 'r', newline='') as cf:
            reader = csv.DictReader(cf)
            for row in reader:
                try:
                    key = (
                        float(row['aox']),
                        float(row['U']),
                        float(row['nH']),
                        float(row['Z']),
                        float(row['cov']),
                        float(row['NH'])
                    )
                except Exception:
                    continue
                rkey = round_key_tuple(key, ndigits=rounding)
                try:
                    intensity = float(row.get('intensity', math.nan))
                except Exception:
                    intensity = math.nan
                mp[rkey] = intensity
        line_maps[name] = mp
        

    # models that have ALL available lines
    key_sets = [set(mp.keys()) for mp in line_maps.values()]
    common_keys = set.intersection(*key_sets) if key_sets else set()

    if not common_keys:
        print("No single model (parameter combination) contains all the available observed lines simultaneously.")
        return None
    

    # model_fluxes (n_models x n_used_lines)
    sorted_keys = sorted(common_keys)
    n_models = len(sorted_keys)
    n_lines = len(available)
    model_fluxes = np.zeros((n_models, n_lines), dtype=float)
    for j, name in enumerate(available):
        mp = line_maps[name]
        for i, k in enumerate(sorted_keys):
            model_fluxes[i, j] = mp.get(k, math.nan)

    obs_fluxes_arr = np.array([v for name,v in zip(obs_names, obs_fluxes) if name in available], dtype=float)
    obs_errors_arr = np.array([v for name,v in zip(obs_names, obs_errors) if name in available], dtype=float)
      
    # compute best-fit single scaling per model and chi2
    #    scal = sum(obs * model / err^2) / sum(model^2 / err^2)
    if np.any(obs_errors_arr == 0):
        raise ValueError("Some observation errors are zero — cannot weight by 1/sigma^2.")
    w = 1.0 / (obs_errors_arr**2)
    numerator = np.nansum((obs_fluxes_arr * model_fluxes) * w, axis=1)
    denominator = np.nansum((model_fluxes**2) * w, axis=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        scal = numerator / denominator
    residuals = obs_fluxes_arr[np.newaxis, :] - (scal[:, np.newaxis] * model_fluxes)
    chi2 = np.nansum((residuals**2) * w, axis=1)

    # choose best
    best_idx = int(np.nanargmin(chi2))
    best_key = sorted_keys[best_idx]
    best_chi2 = float(chi2[best_idx])
    
    try:
        folder_name = make_run_folder_name(*best_key)
        folder_path = os.path.join(models_root, folder_name)
        if not os.path.exists(folder_path):
            folder_path = None
    except Exception:
        folder_path = None

    # report
    print("Best-fit model using available lines:")
    print(f"  used lines (N={len(available)}): {', '.join(available)}")
    print(f"  best chi2 = {best_chi2:.6g}")
    print("  best parameters (aox, U, nh, Z, cov, NH):")
    print("   ", best_key)
    if folder_path:
        print("  best model folder:", folder_path)
    else:
        print("  best model folder not found")

    # top-N list of best models
    order = np.argsort(chi2)
    print(f"\nTop {min(topn, n_models)} models by chi2:")
    for r in range(min(topn, n_models)):
        idx = int(order[r])
        print(f" {r+1:2d}) chi2={chi2[idx]:.6g}  params={sorted_keys[idx]}")

    return {
        'available_lines': available,
        'missing_lines': missing,
        'best_key': best_key,
        'best_chi2': best_chi2,
        'best_folder': folder_path,
        'sorted_keys': sorted_keys,
        'chi2_array': chi2,
        'scale_array': scal,
        'model_fluxes': model_fluxes,
        'obs_fluxes_used': obs_fluxes_arr,
        'obs_errors_used': obs_errors_arr
    }


#######################################################
# CLOUDY parallel running
#######################################################
if __name__ == "__main__":
    if best_switch:
    # obs_lines_file is the variable you used earlier (e.g. "ObsLineList.txt")
        obs_file = obs_lines_file  # or put a path string: "my_obs.txt"
        result = find_best_model(
            obs_file,
            lines_db_dir=os.path.join(output_folder, 'lines_db'),
            models_root=output_folder
            )

        if result is None:
            print("No best model found.")
        else:
            print("Best model params:", result['best_key'])
            print("Best model folder:", result.get('best_folder'))
            print("Best chi2:", result['best_chi2'])
            
    if cloudy_switch != 1:
        print("cloudy_switch != 1 -> not running CLOUDY. Exiting.")
        sys.exit(0)

    # Check for existing combinations in output folder
    existing_folders = set(os.listdir(output_folder))
    duplicates = []
    for combo in combinations:
        name = make_run_folder_name(*combo)
        if name in existing_folders:
            duplicates.append(name)
    if duplicates:
        print("The following models already exist:")
        for d in duplicates[:100]:
            print("  ", d)
        if not overwrite_existing:
            print("They will be skipped unless you set overwrite_existing = True at top of script.")

    # Copy template input and linelist into output root
    try:
        shutil.copy(cloudy_input_file, output_folder)
    except Exception:
        pass
    try:
        shutil.copy(list_lines_file, output_folder)
    except Exception:
        pass

    # Run parallel jobs and show progress
    try:
        from tqdm import tqdm
        HAS_TQDM = True
    except Exception:
        HAS_TQDM = False

    worker = partial(run_single_combination,
                     cloudy_exe=cloudy_executable,
                     cloudy_input_template=cloudy_input_file,
                     list_lines_file_local=list_lines_file,
                     output_root=output_folder,
                     overwrite=overwrite_existing)

    n_jobs = min(max(1, Ncpus), len(combinations))
    
    # If there are fewer models than Ncpus, reduce processes
    print(f"Starting parallel run with {n_jobs} processes...")

    with Pool(processes=n_jobs) as pool:
        iterator = pool.imap_unordered(worker, combinations)

        # If tqdm available use it
        if HAS_TQDM:
            for _ in tqdm(iterator, total=len(combinations), desc="Running models", unit="model"):
                pass
        else:
            for _ in iterator:
                pass

    print("Parallel run finished.")
    




