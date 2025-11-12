import os
import glob

def parse_observation_file(file_path):
    """
    Parse a single observational file to extract data.
    :param file_path: Path to the observational file
    :return: List of tuples (line_name, flux, error)
    """
    galaxy_data = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        # Check if the file has no useful data (only comments or empty lines)
        if all(line.strip().startswith("#") or not line.strip() for line in lines):
            print(f"Warning: File {file_path} contains no useful data. Skipping.")
            return galaxy_data
        
        for line in lines:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) == 3:
                line_name = parts[0]
                flux = float(parts[1])
                error = float(parts[2])
                galaxy_data.append((line_name, flux, error))
    return galaxy_data

def load_galaxy_data(galaxy_folder, galaxy_names, return_obs_files=False):
    """
    Load observational data for a list of galaxies.
    :param galaxy_folder: Base folder containing galaxy subfolders
    :param galaxy_names: List of galaxy names
    :param return_obs_files: Whether to return the list of observational files
    :return: Dictionary of galaxy data and optionally a list of all observational files
    """
    galaxies_data = {}
    all_obs_files = []

    for galaxy_name in galaxy_names:
        galaxy_subfolder = os.path.join(galaxy_folder, galaxy_name)
        obs_files = glob.glob(os.path.join(galaxy_subfolder, "*.txt"))
        all_obs_files.extend(obs_files)  # Collect all files in one list

        for obs_file in obs_files:
            file_name = os.path.basename(obs_file)
            galaxy_data = parse_observation_file(obs_file)
            if galaxy_data:  # Only add non-empty, useful data
                galaxies_data[f"{galaxy_name}/{file_name}"] = galaxy_data
            else:
                print(f"Skipping file: {obs_file}, no useful data found.")

    if return_obs_files:
        return galaxies_data, all_obs_files
    return galaxies_data


