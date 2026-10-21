# SINEFIA

SINEFIA runs grids of [CLOUDY](https://www.nublado.org) photoionization models and compares them with observed emission-line fluxes. The ionizing source can be a stellar population (starburst / HII regions) or an AGN. For every galaxy it finds the best models by χ² and gives the Bayesian posterior of every grid parameter.

## Files

| File | What it is |
|---|---|
| `sinefia_parameters.py` | **The only file you edit.** Every setting, with an explanation: the grid, the ionizing source, the gas, the stops, the paths and the fit options. There is no code in it. |
| `sinefia.py` | All the code. You don't need to edit it. |
| `linelist.dat` | An example CLOUDY line list (mid- and far-infrared lines). |
| `requirements.txt` | The Python packages needed. |

## Requirements

- CLOUDY (tested with c22.02).
- Python 3 with `numpy`, `h5py` and `matplotlib`: `pip install -r requirements.txt`
- For a starburst grid, a stellar population template compiled for CLOUDY (for example PopStar, in a `PopStar/` folder next to the two `.py` files). See SECTION 2a of `sinefia_parameters.py`.

## Quick start

1. Put `sinefia.py`, `sinefia_parameters.py`, your line list and your stellar template in one folder.
2. Edit `sinefia_parameters.py`. At least set `cloudy_executable` and `Ncpus` (SECTION 1), and choose the stages to run in SECTION 0.
3. Check the settings. This prints the grid, the number of models and an example CLOUDY input, and runs nothing:

   ```
   python sinefia.py --check
   ```

4. Run:

   ```
   python sinefia.py
   ```

Another settings file can be given as an argument: `python sinefia.py my_agn_grid.py`.

## The three stages

Each stage is switched on or off in SECTION 0 of `sinefia_parameters.py`.

1. **`RUN_CLOUDY`** writes one CLOUDY input per model, each in its own folder, and runs `Ncpus` models at a time. Finished models are skipped, so an interrupted or extended grid continues where it stopped.
2. **`EXTRACT_LINES`** stores the line fluxes of every finished model in one HDF5 file (`DATABASE_FILE`).
3. **`FIT_OBSERVATIONS`** compares the database with your observed fluxes and writes the results.

Stages 2 and 3 can be run again at any time without running CLOUDY.

## Folder layout

```
models/<GRID_NAME>/<one folder per model>/
    sfmodel.in               CLOUDY input  (agnmodel.* for an AGN grid)
    sfmodel.out              CLOUDY output
    sfmodellineoutput.txt    line fluxes
    params.txt               the model's parameters and settings
models/lines_db.h5           the line database (DATABASE_FILE)
results/                     the fit results (RESULTS_FOLDER)
```

## Observation files

One text file per galaxy in `OBS_FOLDER`, with three columns: line name, flux and error. The line name is the CLOUDY label with the spaces replaced by `_`. Lines starting with `#` are ignored.

```
# line           flux       error
Ne_2_12.8101m    4.4e-16    8.7e-17
O_3_88.3323m     2.1e-16    4.0e-17
```

Lines are matched by species and wavelength, so `Ne_2_12.8101m` and `Ne_2_12.81010m` are the same line. Lines with a nan or negative error are skipped and reported.

## Results

In `RESULTS_FOLDER`:

- `fit_summary.txt`: the best model of every galaxy, one line each.
- `<galaxy>_best_model.txt`: the median and 68% range of every parameter, the best model line by line, the top models and any notes.
- `<galaxy>_posterior_<parameter>.png`: the posterior of each parameter.
- `<galaxy>_all_models.txt`: χ² and probability of every model.
