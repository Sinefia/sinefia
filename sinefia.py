###########################################################
#######################             #######################
#######################   SINEFIA   #######################
#######################             #######################
###########################################################
#
#  Runs the CLOUDY grid defined in sinefia_parameters.py (starburst or
#  AGN ionizing source), stores the line fluxes of every model in an
#  HDF5 database (one entry per model, keyed by the md5 hash of its
#  parameters), and compares the models with observed line fluxes:
#  best models by chi^2 and Bayesian posteriors of every parameter.
#
#  You do not need to edit this file: every setting is in
#  sinefia_parameters.py. Keep both files in the same folder and run
#
#      python sinefia.py --check     check the settings, print the grid
#                                    and an example CLOUDY input
#      python sinefia.py             run the stages switched on in
#                                    SECTION 0 of the parameters
#
#  Another settings file can be given as the last argument, e.g.
#      python sinefia.py --check my_agn_grid.py
#
#  The three stages:
#    1) RUN_CLOUDY        write one CLOUDY input per model (one folder
#                         each) and run them, Ncpus at a time
#    2) EXTRACT_LINES     line fluxes of every finished model -> HDF5
#    3) FIT_OBSERVATIONS  chi^2 ranking + 1-D posteriors per galaxy
#
#  Database layout:  <DATABASE_FILE> / <SB or AGN> / <md5 key>
#      attrs   : the model's parameters, its folder, the line-output
#                modification time, its grid name and settings (JSON)
#      datasets: 'line' (line names) and 'flux' (their fluxes); each
#                model stores only the lines it has
#
#  Needs: numpy, h5py, matplotlib (only for the plots).
#
###########################################################

import os
import re
import sys
import glob
import json
import math
import time
import shutil
import hashlib
import itertools
import subprocess
import importlib.util
from functools import partial
from multiprocessing import Pool, TimeoutError as PoolTimeout

import numpy as np
import h5py

P = None                # the settings module, loaded in main()
BASE_PATH = None        # folder of the settings file

CLOUDY_OK_TEXT = "Cloudy exited OK"

# CLOUDY input / output names in every model folder, per source type
MODEL_PREFIX = {"SB": "sfmodel", "AGN": "agnmodel"}

# The AGN continuum shapes of SECTION 2b, and the CLOUDY line of each
AGN_SED_TYPES = {
    "AGN":             "AGN T={agnT:.6g} K, a(ox)={aox:g}, a(uv)={auv:g}, a(x)={ax:g}",
    "table AGN":       "table AGN break {agn_break:.6g} microns",
    "table power law": "table power law slope {pl_slope:g}, {ir_break:.6g} micron break",
    "table SED":       'table SED "{sed_file}"',
}

# The brightness choices of SECTION 3:
#   mode -> (axis name, CLOUDY command, needs the source distance)
BRIGHTNESS_MODES = {
    "U":                ("U",         "ionization parameter {val}", False),
    "phi(H)":           ("phiH",      "phi(H) {val}",               False),
    "intensity":        ("intensity", "intensity {val}",            False),
    "xi":               ("xi",        "xi {val}",                   False),
    "Q(H)":             ("QH",        "Q(H) {val}",                 True),
    "luminosity":       ("lum",       "luminosity {val}",           True),
    "luminosity total": ("lumtot",    "luminosity total {val}",     True),
}

# Names of the built-in grid axes (CUSTOM_AXES may not reuse them)
BUILTIN_AXES = ({m[0] for m in BRIGHTNESS_MODES.values()} |
                {"nh", "Zgas", "Zstar", "age", "agnT", "aox", "auv", "ax", "agn_break",
                 "pl_slope", "cov", "radius", "cr", "dust", "turb", "bfield",
                 "ff", "stop"})


#######################################################
#  SETTINGS
#######################################################

def load_parameters(path):
    """Import the settings file (sinefia_parameters.py by default)."""
    global P, BASE_PATH
    path = os.path.abspath(path)
    if not os.path.isfile(path):
        sys.exit(f"Settings file not found: {path}")
    sys.dont_write_bytecode = True      # no __pycache__ next to the settings
    spec = importlib.util.spec_from_file_location("sinefia_parameters", path)
    P = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(P)
    BASE_PATH = os.path.dirname(path)


def resolve_path(path):
    """A relative path is taken relative to the settings file."""
    if os.path.isabs(path):
        return path
    return os.path.join(BASE_PATH, path)


def log_scale_axes():
    """Axes given in linear units but spanning decades: their posteriors
    are computed and plotted in log10."""
    return ("age", "agnT", "agn_break") + tuple(
        name for name, spec in P.CUSTOM_AXES.items()
        if isinstance(spec, dict) and spec.get("log_plot"))


def model_files():
    """(input, output, line output) file names inside a model folder."""
    prefix = MODEL_PREFIX[P.SOURCE_TYPE]
    return prefix + ".in", prefix + ".out", prefix + "lineoutput.txt"


#######################################################
#  GRID AND CLOUDY INPUT FILES
#######################################################

def _axis(init, end, step):
    """Return the grid values for one axis.
    step == 0  -> single fixed value [init].
    step != 0  -> init, init+step, ... up to end (inclusive), without
                  the floating-point overshoot of np.arange.
    """
    if step == 0:
        return np.array([float(init)])
    n = (end - init) / step
    if n < -1e-9:
        raise ValueError(
            f"Axis from {init} to {end} with step {step} is empty: "
            "the step has the wrong sign.")
    n = int(math.floor(n + 1e-9))
    return np.round(init + step * np.arange(n + 1), 10)


def _values(spec):
    """Grid values from a list [v1, v2, ...] or a tuple (init, end, step)."""
    if isinstance(spec, tuple):
        if len(spec) != 3:
            raise ValueError(f"{spec}: a tuple must be (init, end, step).")
        return _axis(*spec)
    return np.array(list(spec), dtype=float)


def sed_uses_age():
    """True if the chosen stellar template takes an age."""
    return P.SOURCE_TYPE == "SB" and "{age" in P.SED_LIBRARY[P.SED_TYPE]["command"]


def brightness_axis():
    """Name of the brightness axis (e.g. "U") for BRIGHTNESS_MODE."""
    return BRIGHTNESS_MODES[P.BRIGHTNESS_MODE][0]


def build_grid_axes():
    """The active grid axes, as (dict name -> values, list of names in order)."""
    sw = P.PARAM_SWITCHES
    axes = {}
    axes[brightness_axis()] = _axis(P.brightness_init, P.brightness_end, P.brightness_step)
    axes["nh"] = _axis(P.nh_init, P.nh_end, P.nh_step)
    if sw["metallicity"]:
        axes["Zgas"] = _axis(P.Zgas_init, P.Zgas_end, P.Zgas_step)
    if P.SOURCE_TYPE == "SB":
        sed = P.SED_LIBRARY[P.SED_TYPE]
        # a stellar Z that does not follow the gas is stored as "Zstar"
        if sed.get("z_mode") and not (P.LINK_STELLAR_TO_GAS and "Zgas" in axes):
            axes["Zstar"] = np.array([float(P.Zstar_fixed)])
        if sed_uses_age():
            axes["age"] = np.array(P.age, dtype=float)
        for name, spec in sed.get("axes", {}).items():
            axes[name] = _values(spec)
    elif P.AGN_SED == "AGN":
        axes["agnT"] = np.array(P.agnT, dtype=float)
        axes["aox"] = _axis(P.aox_init, P.aox_end, P.aox_step)
        axes["auv"] = _axis(P.auv_init, P.auv_end, P.auv_step)
        axes["ax"] = _axis(P.ax_init, P.ax_end, P.ax_step)
    elif P.AGN_SED == "table AGN":
        axes["agn_break"] = np.array(P.agn_break, dtype=float)
    elif P.AGN_SED == "table power law":
        axes["pl_slope"] = _axis(P.pl_slope_init, P.pl_slope_end, P.pl_slope_step)
    if sw["covering_factor"]:
        axes["cov"] = _axis(P.covfac_init, P.covfac_end, P.covfac_step)
    if sw["radius"]:
        axes["radius"] = _axis(P.radius_init, P.radius_end, P.radius_step)
    if sw["cosmic_rays"] and P.COSMIC_RAY_MODE == "rate":
        axes["cr"] = _axis(P.cr_init, P.cr_end, P.cr_step)
    if sw["grains"] and sw["dust_to_metals"]:
        axes["dust"] = _axis(P.dust_init, P.dust_end, P.dust_step)
    if sw["turbulence"]:
        axes["turb"] = _axis(P.turb_init, P.turb_end, P.turb_step)
    if sw["magnetic_field"]:
        axes["bfield"] = _axis(P.bfield_init, P.bfield_end, P.bfield_step)
    if sw["filling_factor"]:
        axes["ff"] = _axis(P.ff_init, P.ff_end, P.ff_step)
    for name, spec in P.CUSTOM_AXES.items():
        axes[name] = _values(spec["values"])
    axes["stop"] = _axis(P.stop_init, P.stop_end, P.stop_step)
    return axes, list(axes.keys())


def all_combinations():
    """List of dicts, one per model: {axis name: value}."""
    axes, order = build_grid_axes()
    return [dict(zip(order, (float(v) for v in combo)))
            for combo in itertools.product(*[axes[k] for k in order])]


def validate_parameters():
    """Stop early, with a readable message, on inconsistent settings."""
    problems = []
    sw = P.PARAM_SWITCHES
    if P.SOURCE_TYPE not in ("SB", "AGN"):
        problems.append(f'SOURCE_TYPE must be "SB" or "AGN", not "{P.SOURCE_TYPE}".')
    elif P.SOURCE_TYPE == "AGN":
        if P.AGN_SED not in AGN_SED_TYPES:
            problems.append(f'AGN_SED must be one of {list(AGN_SED_TYPES)}, not "{P.AGN_SED}".')
        elif P.AGN_SED == "AGN":
            if len(P.agnT) == 0:
                problems.append("The agnT list is empty.")
            elif min(P.agnT) <= 10:
                problems.append("agnT must be given in K (linear), e.g. [1.5e5]: "
                                "CLOUDY would read a value <= 10 as log T.")
            if max(P.aox_init, P.aox_end) > 0:
                print("Warning: aox > 0. CLOUDY has no implicit minus sign; "
                      "a typical AGN has aox around -1.4.")
        elif P.AGN_SED == "table AGN":
            if len(P.agn_break) == 0:
                problems.append("The agn_break list is empty.")
            elif min(P.agn_break) < 0.455 or max(P.agn_break) > 3.0e7:
                problems.append("agn_break must be between 0.46 and 3e7 microns "
                                "(0.2 to 3.04e-9 Ryd, Hazy 1 6.12.2).")
        elif P.AGN_SED == "table power law":
            if P.PL_IR_BREAK_MICRONS <= 0:
                problems.append("PL_IR_BREAK_MICRONS must be a positive wavelength in microns.")
        elif P.AGN_SED == "table SED":
            if "/" in P.AGN_SED_FILE and not os.path.exists(resolve_path(P.AGN_SED_FILE)):
                problems.append(f"AGN SED file not found: {resolve_path(P.AGN_SED_FILE)}")
    elif P.SED_TYPE not in P.SED_LIBRARY:
        problems.append(f"SED_TYPE '{P.SED_TYPE}' is not a key of SED_LIBRARY.")
    else:
        sed = P.SED_LIBRARY[P.SED_TYPE]
        if sed.get("z_mode") not in ("parameter", "file", None):
            problems.append("SED z_mode must be 'parameter', 'file' or None.")
        if sed.get("z_mode") == "parameter" and sed.get("z_units") not in ("log_abs", "abs"):
            problems.append("SED z_units must be 'log_abs' or 'abs'.")
        if sed_uses_age() and len(P.age) == 0:
            problems.append("The age list is empty.")
        if (sed.get("z_mode") != "file" and "/" in sed["sed_file"]
                and not os.path.exists(resolve_path(sed["sed_file"]))):
            problems.append(f"Stellar template not found: {resolve_path(sed['sed_file'])}")
        for name in sed.get("axes", {}):
            if not name.isidentifier() or name in BUILTIN_AXES or name in ("sed_file", "Zstar"):
                problems.append(f"SED axis name '{name}' is not allowed (use a new, simple name).")
            elif "{" + name not in sed["command"]:
                problems.append(f"SED axis '{name}' does not appear as {{{name}}} in the command.")
    if P.BRIGHTNESS_MODE not in BRIGHTNESS_MODES:
        problems.append(f'BRIGHTNESS_MODE must be one of {list(BRIGHTNESS_MODES)}, '
                        f'not "{P.BRIGHTNESS_MODE}".')
    elif BRIGHTNESS_MODES[P.BRIGHTNESS_MODE][2] and not sw["radius"]:
        problems.append(f'BRIGHTNESS_MODE "{P.BRIGHTNESS_MODE}" is the source luminosity: CLOUDY '
                        "also needs the distance, so set PARAM_SWITCHES['radius'] = True.")
    taken = set(BUILTIN_AXES)
    if P.SOURCE_TYPE == "SB" and P.SED_TYPE in P.SED_LIBRARY:
        taken |= set(P.SED_LIBRARY[P.SED_TYPE].get("axes", {}))
    for name, spec in P.CUSTOM_AXES.items():
        if not name.isidentifier() or name in taken:
            problems.append(f"CUSTOM_AXES name '{name}' is not allowed (use a new, simple name).")
        elif not isinstance(spec, dict) or "values" not in spec or "command" not in spec:
            problems.append(f"CUSTOM_AXES '{name}' needs a \"values\" and a \"command\".")
        elif "{value" not in spec["command"]:
            problems.append(f"CUSTOM_AXES '{name}': the command has no {{value}}.")
    if P.DENSITY_LAW not in ("hden_constant", "constant_pressure",
                             "constant_gaspressure", "powerlaw_density"):
        problems.append(f"DENSITY_LAW '{P.DENSITY_LAW}' is not one of the listed options.")
    if P.DENSITY_LAW == "powerlaw_density" and not sw["radius"]:
        problems.append("DENSITY_LAW 'powerlaw_density' needs PARAM_SWITCHES['radius'] = True.")
    if P.SPHERE and not sw["radius"]:
        problems.append("SPHERE = True needs PARAM_SWITCHES['radius'] = True.")
    if P.STOP_MODE not in P.STOP_COMMANDS:
        problems.append(f"STOP_MODE '{P.STOP_MODE}' is not a key of STOP_COMMANDS.")
    if P.COSMIC_RAY_MODE not in ("background", "rate"):
        problems.append("COSMIC_RAY_MODE must be 'background' or 'rate'.")
    if not (P.LOWEST_TEMPERATURE in (None, "off") or isinstance(P.LOWEST_TEMPERATURE, (int, float))):
        problems.append('LOWEST_TEMPERATURE must be a number (K), "off" or None.')
    if P.ITERATE not in (None, "iterate", "iterate to convergence"):
        problems.append("ITERATE must be 'iterate', 'iterate to convergence' or None.")
    if sw["filling_factor"] and max(P.ff_init, P.ff_end) > 0:
        problems.append("The log filling factor must be <= 0 (filling factor <= 1).")
    if sw["covering_factor"] and not (0 < min(P.covfac_init, P.covfac_end)
                                      and max(P.covfac_init, P.covfac_end) <= 1):
        problems.append("The covering factor must be between 0 and 1.")
    if not os.path.exists(resolve_path(P.list_lines_filename)):
        problems.append(f"Line list not found: {resolve_path(P.list_lines_filename)}")
    if (P.FIT_OBSERVATIONS and P.FIT_LINELIST_FILENAME
            and not os.path.exists(resolve_path(P.FIT_LINELIST_FILENAME))):
        problems.append(f"Fit line list not found: {resolve_path(P.FIT_LINELIST_FILENAME)} "
                        "(FIT_LINELIST_FILENAME, SECTION 12; None uses every observed line)")
    if int(P.MIN_LINES) < 1:
        problems.append("MIN_LINES must be at least 1.")
    if not problems:
        try:
            build_grid_axes()
            if P.SOURCE_TYPE == "SB":
                for combo in all_combinations():
                    sed_command(combo)
        except ValueError as err:
            problems.append(str(err))
    if problems:
        sys.exit("SINEFIA parameter problems:\n  - " + "\n  - ".join(problems))


def stellar_metallicity(pdict):
    """Relative stellar metallicity log10(Z*/Zsun) for one model."""
    if "Zstar" in pdict:
        return pdict["Zstar"]
    if P.LINK_STELLAR_TO_GAS and "Zgas" in pdict:
        return pdict["Zgas"]
    return P.Zstar_fixed


def sed_command(pdict):
    """The 'table star' line for one model."""
    sed = P.SED_LIBRARY[P.SED_TYPE]
    fmt = {name: pdict[name] for name in sed.get("axes", {})}
    if "age" in pdict:
        fmt["age"] = pdict["age"]
    zrel = stellar_metallicity(pdict)
    z_abs = sed["Zsun_template"] * 10.0 ** zrel
    sed_file = sed["sed_file"]

    if sed.get("z_mode") == "parameter":
        if sed["z_units"] == "log_abs":
            fmt["Zstar"] = f"{math.log10(z_abs):.5f}"
        else:
            fmt["Zstar"] = f"{z_abs:.6g}"
    elif sed.get("z_mode") == "file":
        available = sed.get("available_Z", [])
        matches = [z for z in available if abs(z - z_abs) / z < 0.01]
        if available and not matches:
            raise ValueError(
                f"Stellar Z = {z_abs:.4g} (log Z/Zsun = {zrel:.3f}) has no "
                f"template file; available Z: {available}")
        z_use = matches[0] if matches else z_abs
        sed_file = sed_file.replace("{Z4}", f"{int(round(z_use * 1e4)):04d}")

    if "/" in sed_file:
        sed_file = resolve_path(sed_file)
    fmt["sed_file"] = sed_file
    return sed["command"].format(**fmt)


def source_command(pdict):
    """The line that sets the shape of the ionizing continuum."""
    if P.SOURCE_TYPE == "AGN":
        sed_file = P.AGN_SED_FILE
        if "/" in sed_file:
            sed_file = resolve_path(sed_file)
        return AGN_SED_TYPES[P.AGN_SED].format(
            ir_break=P.PL_IR_BREAK_MICRONS, sed_file=sed_file, **pdict)
    return sed_command(pdict)


def model_folder_name(pdict):
    """Readable, unique folder name built from every parameter value."""
    return "_".join(f"{k}{v:.6g}" for k, v in pdict.items())


def build_input_deck(pdict):
    """The CLOUDY input of ONE model, as a list of lines."""
    sw = P.PARAM_SWITCHES
    lines = []
    lines.append("title SINEFIA " + " ".join(f"{k}={v:.6g}" for k, v in pdict.items()))

    # ---- ionizing source (stars or AGN) and its intensity ----
    # The CMB sets its own shape AND intensity, so it must come after the
    # source's shape + intensity pair, never between them (Hazy 1, 4.4.2).
    lines.append(source_command(pdict))
    lines.append(BRIGHTNESS_MODES[P.BRIGHTNESS_MODE][1].format(val=pdict[brightness_axis()]))
    if P.ADD_CMB:
        lines.append(f"CMB redshift {P.CMB_REDSHIFT}")

    # ---- cosmic rays ----
    if sw["cosmic_rays"]:
        if P.COSMIC_RAY_MODE == "background":
            if P.COSMIC_RAY_BACKGROUND_FACTOR == 1.0:
                lines.append("cosmic rays background")
            else:
                lines.append(f"cosmic rays background {P.COSMIC_RAY_BACKGROUND_FACTOR:.6g} linear")
        elif P.COSMIC_RAY_MODE == "rate":
            lines.append(f"cosmic ray rate {pdict['cr']}")

    # ---- density / pressure law ----
    if P.DENSITY_LAW == "powerlaw_density":
        lines.append(f"hden {pdict['nh']}, power = {P.HDEN_POWERLAW_INDEX}")
    else:
        lines.append(f"hden {pdict['nh']}")
    if P.DENSITY_LAW == "hden_constant":
        lines.append("constant density")
    elif P.DENSITY_LAW == "constant_pressure":
        lines.append("constant pressure")
    elif P.DENSITY_LAW == "constant_gaspressure":
        lines.append("constant gas pressure")

    # ---- chemical composition ----
    if sw["metallicity"]:
        lines.append(f"abundances {P.ABUNDANCE_SET} no grains")
        # CLOUDY reads a number <= 0 as log10 and a positive number as a
        # LINEAR factor (Hazy 1, 7.10). So a positive log value is written
        # as its factor.
        z = pdict["Zgas"]
        if z > 0:
            lines.append(f"metals and grains {10.0 ** z:.6g}")
        else:
            lines.append(f"metals and grains {z}")
    for element, value in P.ELEMENT_ABUNDANCES.items():
        lines.append(f"element abundance {element} {value}")
    if sw["grains"]:
        for g in P.GRAIN_TYPES:
            if "dust" in pdict:
                lines.append(f"grains {g} {10.0 ** pdict['dust']:.6g} linear")
            else:
                lines.append(f"grains {g}")

    # ---- geometry ----
    if "radius" in pdict:
        lines.append(f"radius {pdict['radius']}")
    if P.SPHERE:
        lines.append("sphere")
    if "cov" in pdict:
        lines.append(f"covering factor {pdict['cov']}")
    if "ff" in pdict:
        lines.append(f"filling factor {pdict['ff']} {P.FILLING_INDEX}")

    # ---- extra microphysics ----
    if "turb" in pdict:
        lines.append(f"turbulence {pdict['turb']} km/s")
    if "bfield" in pdict:
        lines.append(f"magnetic field {pdict['bfield']}")

    # ---- your own axes (SECTION 9b) ----
    for name, spec in P.CUSTOM_AXES.items():
        lines.append(spec["command"].format(value=pdict[name]))

    # ---- stopping criteria and iterations ----
    lines.append(P.STOP_COMMANDS[P.STOP_MODE].format(val=pdict["stop"]))
    if P.LOWEST_TEMPERATURE == "off":
        lines.append("stop temperature off")
    elif P.LOWEST_TEMPERATURE is not None:
        lines.append(f"stop temperature {P.LOWEST_TEMPERATURE} K linear")
    lines.extend(P.EXTRA_STOPS)
    if P.ITERATE:
        lines.append(P.ITERATE)
    lines.extend(P.EXTRA_COMMANDS)

    # ---- output ----
    lineoutput = model_files()[2]
    lines.append(f'save last line list {P.SAVE_LINE_OPTIONS} '
                 f'"{lineoutput}" "{os.path.basename(P.list_lines_filename)}"')
    if P.SAVE_OVERVIEW:
        lines.append(f'save last overview "{MODEL_PREFIX[P.SOURCE_TYPE]}overview.txt"')
    return lines


def grid_summary():
    """Print a human-readable summary of the grid in the settings."""
    axes, order = build_grid_axes()
    total = 1
    print("=" * 55)
    print(" SINEFIA grid configuration")
    print("=" * 55)
    print(f" Grid name           : {P.GRID_NAME}")
    print(f" Ionizing source     : {P.SOURCE_TYPE}")
    print(f" Brightness          : {P.BRIGHTNESS_MODE} (axis '{brightness_axis()}')")
    if P.SOURCE_TYPE == "AGN":
        extra = {"table power law": f" (infrared break {P.PL_IR_BREAK_MICRONS:g} microns)",
                 "table SED": f' "{P.AGN_SED_FILE}"'}.get(P.AGN_SED, "")
        print(f" AGN continuum       : {P.AGN_SED}{extra}")
    else:
        print(f" SED type            : {P.SED_TYPE}")
    if P.ADD_CMB:
        print(f" CMB added           : redshift {P.CMB_REDSHIFT}")
    print(f" Density law         : {P.DENSITY_LAW}")
    print(f" Abundance set       : {P.ABUNDANCE_SET}")
    if P.SOURCE_TYPE == "SB" and P.SED_LIBRARY[P.SED_TYPE].get("z_mode"):
        print(f" Stellar Z           : "
              f"{'linked to gas Z' if P.LINK_STELLAR_TO_GAS else f'fixed at {P.Zstar_fixed}'}")
    print(f" Stop mode           : {P.STOP_MODE}")
    print(f" Lowest temperature  : "
          f"{P.LOWEST_TEMPERATURE if P.LOWEST_TEMPERATURE is not None else 'CLOUDY default (4000 K)'}")
    print(f" Iterations          : {P.ITERATE}")
    print("-" * 55)
    print(" Grid axes (name : n points : values)")
    for name in order:
        vals = axes[name]
        total *= len(vals)
        shown = ", ".join(f"{v:g}" for v in vals[:6])
        if len(vals) > 6:
            shown += ", ..."
        print(f"   {name:9s} : {len(vals):3d} : [{shown}]")
    print("-" * 55)
    print(f" TOTAL MODELS        : {total}")
    print("=" * 55)
    return total


def grid_settings():
    """Every fixed (non-gridded) choice that changes the CLOUDY physics.
    It is stored with each model (params.txt and the database), so models
    made with different settings are never confused."""
    sw = P.PARAM_SWITCHES
    s = {"SOURCE_TYPE": P.SOURCE_TYPE, "BRIGHTNESS_MODE": P.BRIGHTNESS_MODE}
    if P.SOURCE_TYPE == "AGN":
        s["AGN_SED"] = P.AGN_SED
        if P.AGN_SED == "table power law":
            s["PL_IR_BREAK_MICRONS"] = P.PL_IR_BREAK_MICRONS
        elif P.AGN_SED == "table SED":
            s["AGN_SED_FILE"] = P.AGN_SED_FILE
    else:
        sed = P.SED_LIBRARY[P.SED_TYPE]
        s["SED_TYPE"] = P.SED_TYPE
        s["SED"] = {k: sed.get(k) for k in ("command", "sed_file", "z_mode",
                                            "z_units", "Zsun_template")}
        s["LINK_STELLAR_TO_GAS"] = P.LINK_STELLAR_TO_GAS
    s["CMB_REDSHIFT"] = P.CMB_REDSHIFT if P.ADD_CMB else None
    s["COSMIC_RAYS"] = (None if not sw["cosmic_rays"] else
                        P.COSMIC_RAY_BACKGROUND_FACTOR if P.COSMIC_RAY_MODE == "background"
                        else "rate")
    s["DENSITY_LAW"] = P.DENSITY_LAW
    if P.DENSITY_LAW == "powerlaw_density":
        s["HDEN_POWERLAW_INDEX"] = P.HDEN_POWERLAW_INDEX
    s["ABUNDANCE_SET"] = P.ABUNDANCE_SET if sw["metallicity"] else None
    s["ELEMENT_ABUNDANCES"] = P.ELEMENT_ABUNDANCES
    s["GRAIN_TYPES"] = P.GRAIN_TYPES if sw["grains"] else None
    s["SPHERE"] = P.SPHERE
    if sw["filling_factor"]:
        s["FILLING_INDEX"] = P.FILLING_INDEX
    s["CUSTOM_AXES"] = {k: v["command"] for k, v in P.CUSTOM_AXES.items()}
    s["STOP"] = P.STOP_COMMANDS[P.STOP_MODE]
    s["LOWEST_TEMPERATURE"] = P.LOWEST_TEMPERATURE
    s["EXTRA_STOPS"] = P.EXTRA_STOPS
    s["ITERATE"] = P.ITERATE
    s["EXTRA_COMMANDS"] = P.EXTRA_COMMANDS
    s["SAVE_LINE_OPTIONS"] = P.SAVE_LINE_OPTIONS
    # JSON round trip, so the dict compares equal to one read back from disk
    return json.loads(json.dumps(s))


def settings_differences(old, new):
    """Readable list of the settings that differ between two dicts."""
    return [f"{k}: {old.get(k)!r} -> {new.get(k)!r}"
            for k in sorted(set(old) | set(new)) if old.get(k) != new.get(k)]


#######################################################
#  SMALL HELPERS
#######################################################

def sanitize_line_name(name):
    """CLOUDY line label -> name stored in the database.
    'O  3 5006.84A' -> 'O_3_5006.84A'"""
    s = name.strip()
    s = re.sub(r'[^\w\s\-\.]', '', s)       # remove some characters
    s = re.sub(r'\s+', '_', s)              # replace spaces with one underscore
    return s


def line_key(name):
    """Name-independent identity of a line, so that 'Ne 2 12.8101m',
    'Ne_2_12.8101m' and 'Ne_2_12.81010m' all match:
    (species, unit, wavelength to 6 significant digits)."""
    tokens = re.split(r"[\s_]+", name.strip())
    m = re.match(r"^([0-9.]+)([A-Za-z]*)$", tokens[-1]) if len(tokens) > 1 else None
    if not m:
        return (name.strip().lower(), "", "")
    return (" ".join(tokens[:-1]).lower(), m.group(2).lower(), f"{float(m.group(1)):.6g}")


def read_line_list(path):
    """CLOUDY line labels from a line-list file (skips blanks, # lines
    and # comments at the end of a line)."""
    labels = []
    with open(path, "r") as f:
        for line in f:
            line = line.split("#")[0].strip()
            if line:
                labels.append(line)
    return labels


class Progress:
    """A one-line progress bar, redrawn in place:
      [##########----------]  45/189  24%  0:01:02 elapsed, ~0:03:10 left  <item>"""

    def __init__(self, total, title, guess_seconds=None):
        """guess_seconds: expected total time, used for the time left
        until the first item is finished."""
        self.total = max(1, total)
        self.done = 0
        self.start = time.time()
        self.guess = guess_seconds
        self.item = ""
        print(title)
        if guess_seconds is not None:
            self.draw()

    @staticmethod
    def _hms(seconds):
        seconds = int(round(seconds))
        return f"{seconds // 3600}:{seconds % 3600 // 60:02d}:{seconds % 60:02d}"

    def update(self, item=""):
        """Call once per finished item; 'item' is shown at the end of the bar."""
        self.done += 1
        self.item = item
        self.draw()

    def draw(self):
        """Redraw the bar (also between items, so the clock keeps running)."""
        elapsed = time.time() - self.start
        if self.done:
            left = f"~{self._hms(elapsed / self.done * (self.total - self.done))} left"
        elif self.guess is not None:
            left = f"~{self._hms(max(0.0, self.guess - elapsed))} left (rough estimate)"
        else:
            left = ""
        n = int(20 * self.done / self.total)
        line = (f"  [{'#' * n}{'-' * (20 - n)}] {self.done:>{len(str(self.total))}}/{self.total} "
                f"{100 * self.done / self.total:3.0f}%  {self._hms(elapsed)} elapsed, "
                f"{left}  {self.item}")
        width = shutil.get_terminal_size((100, 20)).columns - 1
        sys.stdout.write("\r" + line[:width].ljust(width))
        sys.stdout.flush()

    def close(self):
        sys.stdout.write("\n")
        sys.stdout.flush()


def cloudy_command():
    """The CLOUDY executable: a path with a "/" is taken relative to the
    settings file, a bare name is looked up on the system PATH."""
    exe = P.cloudy_executable
    if os.sep in exe or "/" in exe:
        exe = resolve_path(exe)
        return exe if os.path.isfile(exe) else None
    return shutil.which(exe)


def grid_folder():
    return os.path.join(resolve_path(P.MODELS_FOLDER), P.GRID_NAME)


def is_finished(folder, out_file, lineoutput_file):
    """A model is finished if CLOUDY ended normally and wrote its lines."""
    out = os.path.join(folder, out_file)
    if not (os.path.exists(out) and os.path.exists(os.path.join(folder, lineoutput_file))):
        return False
    with open(out, "r", errors="replace") as f:
        return CLOUDY_OK_TEXT in f.read()


def write_params_txt(folder, pdict, settings):
    """params.txt: one 'name = value' line per parameter, plus the
    settings (JSON) and the grid name."""
    with open(os.path.join(folder, "params.txt"), "w") as f:
        f.write("# SINEFIA model parameters\n")
        for k, v in pdict.items():
            f.write(f"{k} = {v!r}\n")
        f.write(f"grid_name = {P.GRID_NAME}\n")
        f.write(f"settings = {json.dumps(settings, sort_keys=True)}\n")


def params_txt(path):
    """Read params.txt of a model folder (numbers as floats)."""
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
                except ValueError:
                    params[k] = v
    return params


def split_params(info):
    """(grid parameters, settings dict, grid name) from a params.txt dict."""
    params = {k: v for k, v in info.items() if k not in ("settings", "grid_name")}
    try:
        settings = json.loads(info.get("settings", "null"))
    except ValueError:
        settings = None
    return params, settings, info.get("grid_name", "")


#######################################################
#  STAGE 1 : RUN CLOUDY
#######################################################

def prepare_model(pdict, settings):
    """Create the model folder, its CLOUDY input, params.txt and line list."""
    folder = os.path.join(grid_folder(), model_folder_name(pdict))
    os.makedirs(folder, exist_ok=True)
    # remove the outputs of an earlier run, so a failed rerun never
    # leaves old fluxes behind
    for old in glob.glob(os.path.join(folder, "*")):
        if os.path.isfile(old):
            os.remove(old)
    input_file = model_files()[0]
    with open(os.path.join(folder, input_file), "w") as f:
        f.write("\n".join(build_input_deck(pdict)) + "\n")
    write_params_txt(folder, pdict, settings)
    shutil.copy(resolve_path(P.list_lines_filename), folder)
    return folder


def run_single_model(folder, cloudy_exe, input_file, out_file, lineoutput_file):
    """Run CLOUDY inside the model folder:  cloudy.exe < model.in > model.out"""
    with open(os.path.join(folder, input_file)) as fin, \
            open(os.path.join(folder, out_file), 'w') as fout:
        r = subprocess.run([cloudy_exe], stdin=fin, stdout=fout,
                           stderr=subprocess.STDOUT, cwd=folder)
    return folder, is_finished(folder, out_file, lineoutput_file), r.returncode


def run_models():
    exe = cloudy_command()
    if exe is None:
        sys.exit(f"Error: the CLOUDY executable '{P.cloudy_executable}' was not found "
                 "(cloudy_executable, SECTION 1).")
    input_file, out_file, lineoutput_file = model_files()
    combos = all_combinations()
    settings = grid_settings()

    # A finished model made with other fixed settings has the same folder
    # name but different physics: never count it as done.
    todo, skipped = [], 0
    for pdict in combos:
        folder = os.path.join(grid_folder(), model_folder_name(pdict))
        if is_finished(folder, out_file, lineoutput_file) and not P.OVERWRITE_EXISTING:
            pfile = os.path.join(folder, "params.txt")
            old = split_params(params_txt(pfile))[1] if os.path.isfile(pfile) else None
            if old != settings:
                print(f"Error: {grid_folder()} already holds finished models made with "
                      "different settings:")
                for d in settings_differences(old or {}, settings):
                    print("    " + d)
                sys.exit("Use a new GRID_NAME for the new settings, or set "
                         "OVERWRITE_EXISTING = True to run these models again.")
            skipped += 1
            continue
        todo.append(pdict)

    n_jobs = min(max(1, int(P.Ncpus)), max(1, len(todo)))
    hours = len(todo) * P.AVERAGE_MINUTES_PER_MODEL / n_jobs / 60
    if not todo:
        print(f"\n1) CLOUDY: all {len(combos)} models are already finished.")
        return

    folders = [prepare_model(pdict, settings) for pdict in todo]
    worker = partial(run_single_model, cloudy_exe=exe, input_file=input_file,
                     out_file=out_file, lineoutput_file=lineoutput_file)
    bar = Progress(len(folders), f"\n1) CLOUDY: running {len(todo)} models with {n_jobs} "
                                 f"processes ({skipped} already finished; rough estimate "
                                 f"{Progress._hms(hours * 3600)} h:mm:ss)")
    bar.guess = hours * 3600
    bar.draw()
    failed = []
    with Pool(processes=n_jobs) as pool:
        results = pool.imap_unordered(worker, folders)
        while bar.done < len(folders):
            try:
                folder, ok, code = results.next(timeout=5)
            except PoolTimeout:          # no model finished yet: just move the clock
                bar.draw()
                continue
            if not ok:
                failed.append(folder)
            bar.update("last finished: " + os.path.basename(folder) + ("" if ok else "  FAILED"))
    bar.close()
    if failed:
        print(f"   {len(failed)} model(s) failed or gave no line output (check their {out_file}):")
        for folder in failed[:100]:
            print("     ", os.path.relpath(folder, BASE_PATH))


#######################################################
#  STAGE 2 : EXTRACT LINES -> HDF5 DATABASE
#######################################################

def lineoutput_txt(path):
    """(line names, fluxes) of the last iteration in a CLOUDY line output.
    The names come from the file's own '#' header line."""
    with open(path, 'r') as f:
        lines = [ln.rstrip('\n') for ln in f if ln.strip()]
    header_line = next((ln for ln in lines if ln.lstrip().startswith('#')), None)
    iter_lines = [ln for ln in lines if ln.lower().lstrip().startswith('iteration')]
    if header_line is None or not iter_lines:
        return [], []
    header = [sanitize_line_name(t) for t in header_line.lstrip('#').split('\t')[1:]]
    values = []
    for v in iter_lines[-1].split()[2:]:            # skip 'iteration N'
        try:
            values.append(float(v))
        except ValueError:
            values.append(math.nan)
    return header, values


def model_key(params, settings):
    """Unique model id: md5 hash of its parameters and fixed settings."""
    vals = {k: (f"{v:.6g}" if isinstance(v, float) else str(v)) for k, v in params.items()}
    vals["_settings"] = settings
    return hashlib.md5(json.dumps(vals, sort_keys=True).encode()).hexdigest()


def update_lines_database():
    """Store every finished model of this grid in the HDF5 database.
    Models already stored and unchanged are skipped; rerun models (newer
    line output) are replaced."""
    input_file, out_file, lineoutput_file = model_files()
    root = grid_folder()
    if not os.path.isdir(root):
        print(f"\n2) No models found in {root}. Run CLOUDY first (RUN_CLOUDY = True).")
        return
    db_file = resolve_path(P.DATABASE_FILE)
    os.makedirs(os.path.dirname(db_file) or ".", exist_ok=True)

    added = updated = unchanged = not_finished = removed = 0
    warnings, duplicates = [], []
    entries = [e for e in sorted(os.listdir(root))
               if os.path.isfile(os.path.join(root, e, 'params.txt'))]
    bar = Progress(len(entries), f"\n2) Storing the line fluxes of {len(entries)} model "
                                 f"folders in {os.path.relpath(db_file, BASE_PATH)}")
    with h5py.File(db_file, 'a') as db:
        # which stored models came from which folder
        stored = {}
        for gname in db:
            for k in db[gname]:
                stored.setdefault(str(db[gname][k].attrs['folder']), []).append((gname, k))

        for entry in entries:
            bar.update(entry)
            folder = os.path.join(root, entry)
            params_path = os.path.join(folder, 'params.txt')
            rel_folder = os.path.relpath(folder, BASE_PATH)
            finished = is_finished(folder, out_file, lineoutput_file)
            params, settings, grid_name = split_params(params_txt(params_path))
            key = model_key(params, settings) if settings else None
            # A folder that was rerun (failed, still running, or with other
            # settings) no longer holds the fluxes stored from it: drop them.
            for gname, k in stored.get(rel_folder, []):
                if not (finished and k == key):
                    del db[gname][k]
                    removed += 1
            if not finished:
                not_finished += 1
                continue
            if not settings:
                warnings.append(f"{rel_folder}: params.txt has no settings, skipped")
                continue
            grp = db.require_group(settings["SOURCE_TYPE"])
            lineoutput_path = os.path.join(folder, lineoutput_file)
            mtime = os.path.getmtime(lineoutput_path)
            if key in grp:
                if grp[key].attrs['folder'] != rel_folder:
                    duplicates.append((rel_folder, str(grp[key].attrs['folder'])))
                    continue
                if grp[key].attrs['mtime'] == mtime:      # already stored and unchanged
                    unchanged += 1
                    continue

            header, values = lineoutput_txt(lineoutput_path)
            line_map = {h: v for h, v in zip(header, values) if h and np.isfinite(v)}
            if not line_map:
                warnings.append(f"{rel_folder}: no line fluxes read, skipped")
                continue

            if key in grp:
                del grp[key]
                updated += 1
            else:
                added += 1
            g = grp.create_group(key)
            g.attrs['folder'] = rel_folder
            g.attrs['mtime'] = mtime
            g.attrs['grid_name'] = grid_name
            g.attrs['settings'] = json.dumps(settings, sort_keys=True)
            g.attrs['param_names'] = json.dumps(list(params))
            for k, v in params.items():
                g.attrs[k] = v
            g.create_dataset('line', data=list(line_map.keys()), dtype=h5py.string_dtype())
            g.create_dataset('flux', data=np.array(list(line_map.values()), dtype=float))
        ntot = len(db[P.SOURCE_TYPE]) if P.SOURCE_TYPE in db else 0

    bar.close()
    print(f"   {added} added, {updated} updated, {unchanged} unchanged, "
          f"{not_finished} not finished (skipped); {ntot} {P.SOURCE_TYPE} models in the database")
    if removed:
        print(f"   Removed {removed} stored model(s) whose folder was run again "
              "and holds no finished model with the same settings any more")
    if duplicates:
        print(f"   Warning: {len(duplicates)} model(s) skipped: a model with the same parameters "
              f"and settings is already stored from another folder, e.g.\n"
              f"      {duplicates[0][0]}\n      is the same as {duplicates[0][1]}")
    for w in warnings[:5]:
        print("   Warning: " + w)
    if len(warnings) > 5:
        print(f"   ... and {len(warnings) - 5} more warnings like these")


def load_models():
    """Models of the database made with the CURRENT settings and grid
    parameters (from any GRID_NAME). Returns (param names, list of param
    dicts, list of {line_key: (name, flux)}, folders, n in group)."""
    db_file = resolve_path(P.DATABASE_FILE)
    if not os.path.exists(db_file):
        print("Database not found:", db_file)
        return None
    settings = json.dumps(grid_settings(), sort_keys=True)
    _, param_names = build_grid_axes()
    params, lines, folders = [], [], []
    with h5py.File(db_file, 'r') as db:
        if P.SOURCE_TYPE not in db or len(db[P.SOURCE_TYPE]) == 0:
            print(f"No {P.SOURCE_TYPE} models in", db_file)
            return None
        grp = db[P.SOURCE_TYPE]
        n_group = len(grp)
        for k in sorted(grp.keys()):
            g = grp[k]
            if g.attrs.get('settings') != settings:
                continue
            if sorted(json.loads(g.attrs['param_names'])) != sorted(param_names):
                continue
            params.append({p: float(g.attrs[p]) for p in param_names})
            lines.append({line_key(n): (n, f) for n, f in
                          zip(g['line'].asstr()[()], g['flux'][()])})
            folders.append(str(g.attrs['folder']))
    # sort by parameter values, so the output tables are easy to read
    order = sorted(range(len(params)), key=lambda i: [params[i][p] for p in param_names])
    params = [params[i] for i in order]
    lines = [lines[i] for i in order]
    folders = [folders[i] for i in order]
    return param_names, params, lines, folders, n_group


#######################################################
#  STAGE 3 : COMPARISON WITH OBSERVATIONS
#######################################################

def read_observations(obs_file, allowed):
    """(names, fluxes, errors, notes) of the allowed lines in one galaxy
    file; notes lists the rows that could not be used."""
    names, values, errors, notes = [], [], [], []
    with open(obs_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 3:
                continue
            name = parts[0]
            if allowed is not None and line_key(name) not in allowed:
                continue
            try:
                flux, err = float(parts[1]), float(parts[2])
            except ValueError:
                notes.append(f"'{line.strip()}' is not a row of: name flux error (skipped)")
                continue
            if not (np.isfinite(flux) and np.isfinite(err)) or err <= 0:
                notes.append(f"{name}: flux {flux}, error {err} not usable (skipped)")
                continue
            names.append(name)
            values.append(flux)
            errors.append(err)
    return names, np.array(values), np.array(errors), notes


def weighted_percentiles(values, probs, qs=(0.16, 0.50, 0.84)):
    """Percentiles of a discrete distribution (grid values + probabilities)."""
    cdf = np.cumsum(probs) - 0.5 * probs      # mid-point CDF
    return [float(np.interp(q, cdf, values)) for q in qs]


def chi2_of(obs, sigma, model, has):
    """Best scale (if FREE_SCALING) and chi^2 of every model, using only
    the lines each model has (has = boolean mask, n_models x n_lines)."""
    w = np.where(has, 1.0 / sigma ** 2, 0.0)
    F = np.where(has, model, 0.0)
    with np.errstate(divide='ignore', invalid='ignore'):
        if P.FREE_SCALING:
            scale = np.sum(obs * F * w, axis=1) / np.sum(F ** 2 * w, axis=1)
        else:
            scale = np.ones(len(F))
        chi2 = np.sum(w * (obs - scale[:, None] * F) ** 2, axis=1)
    chi2[~np.isfinite(chi2)] = np.inf
    return scale, chi2


def galaxy_name(obs_file):
    """The galaxy name: the file name without the OBS_FILE_PATTERN suffix."""
    galaxy = os.path.basename(obs_file)
    suffix = P.OBS_FILE_PATTERN.replace("*", "")
    return galaxy[:-len(suffix)] if suffix and galaxy.endswith(suffix) else os.path.splitext(galaxy)[0]


def fit_galaxy(obs_file, allowed, grid):
    """Fit one galaxy and write its result files. Returns (best-model
    parameters or None if skipped, reduced chi^2, number of lines used,
    notes)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    galaxy = galaxy_name(obs_file)
    names, obs, err, notes = read_observations(obs_file, allowed)
    if not names:
        return None, np.nan, 0, notes + ["no usable observed lines: skipped"]
    param_names, params, model_lines, folders, _ = grid

    # model fluxes (n_models x n_obs_lines), NaN where a model lacks a line
    keys = [line_key(n) for n in names]
    model = np.array([[ml[k][1] if k in ml else np.nan for k in keys] for ml in model_lines])
    has = np.isfinite(model)
    missing = [n for j, n in enumerate(names) if not has[:, j].any()]
    if missing:
        notes.append("observed lines in no model (ignored): " + ", ".join(missing))

    sigma = np.sqrt(err ** 2 + (P.SYSTEMATIC_ERROR * obs) ** 2)
    n_scale = 1 if P.FREE_SCALING else 0

    # ---- 1) rank the models by reduced chi^2 (models may have different lines)
    scale, chi2 = chi2_of(obs, sigma, model, has)
    nlines = has.sum(axis=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        chi2_red = np.where((nlines >= int(P.MIN_LINES)) & (nlines - n_scale > 0),
                            chi2 / (nlines - n_scale), np.nan)
    chi2_red[~np.isfinite(chi2_red)] = np.nan
    n_ranked = int(np.isfinite(chi2_red).sum())
    if n_ranked == 0:
        return None, np.nan, 0, notes + [f"no model has at least {P.MIN_LINES} of the observed "
                                         "lines (MIN_LINES): skipped"]
    order = np.argsort(chi2_red)               # NaN go last
    best = int(order[0])
    used_mask = has[best]
    used = [n for n, u in zip(names, used_mask) if u]
    not_in_best = [n for j, n in enumerate(names) if not used_mask[j] and n not in missing]
    if not_in_best:
        notes.append("observed lines not in the best model (not used): " + ", ".join(not_in_best))

    # ---- 2) posteriors over the models that have every line of the best model,
    #         with chi^2 on exactly those lines (flat prior on the grid points)
    complete = np.all(has[:, used_mask], axis=1)
    sub_has = np.broadcast_to(used_mask, has.shape) & complete[:, None]
    sub_scale, sub_chi2 = chi2_of(obs, sigma, np.where(sub_has, model, np.nan), sub_has)
    sub_chi2[~complete] = np.inf
    like = np.exp(-0.5 * (sub_chi2 - sub_chi2[best]))
    like[~np.isfinite(like)] = 0.0
    like /= like.sum()

    keys_arr = np.array([[p[n] for n in param_names] for p in params])
    sel = complete
    # Free parameters: those that add new models. A parameter that only
    # moves together with another one (e.g. Zstar with Zgas when two
    # grids with different metallicities share the database) is not free.
    varied, n_combos = [], 1
    for j, p in enumerate(param_names):
        cols = [param_names.index(v) for v in varied] + [j]
        n = len({tuple(row) for row in keys_arr[sel][:, cols]})
        if n > n_combos:
            varied.append(p)
            n_combos = n
    n_free = len(varied) + n_scale
    dof = len(used) - n_free
    if dof <= 0:
        notes.append(f"{len(used)} line(s) for {n_free} free parameters "
                     f"({', '.join(varied)}{' + flux scale' if n_scale else ''}): the fit cannot "
                     "constrain them; the posteriors mostly reflect the grid")

    results_folder = resolve_path(P.RESULTS_FOLDER)
    os.makedirs(results_folder, exist_ok=True)
    summary = {}
    for j, p in enumerate(param_names):
        values = np.array(sorted(set(keys_arr[sel, j])))
        post = np.array([like[sel & (keys_arr[:, j] == v)].sum() for v in values])
        mle = float(keys_arr[best, j])
        if len(values) == 1:
            summary[p] = (mle, mle, mle, mle)
            continue
        in_log = p in log_scale_axes() and values.min() > 0   # e.g. age in log years
        x = np.log10(values) if in_log else values
        lo, med, hi = weighted_percentiles(x, post)
        if in_log:
            lo, med, hi = 10 ** lo, 10 ** med, 10 ** hi
        # round away floating-point noise (1e-16 -> 0) before writing
        lo, med, hi = (float(np.round(v, 10)) + 0.0 for v in (lo, med, hi))
        summary[p] = (med, lo, hi, mle)

        plt.figure()
        plt.plot(values, post, "ko-", markersize=6, label="Posterior")
        plt.axvline(med, color="red", label=f"Median: {med:.3g}")
        plt.axvspan(lo, hi, color="blue", alpha=0.15, label="68% interval")
        plt.axvline(mle, color="orange", label=f"Best model: {mle:.3g}")
        if in_log:
            plt.xscale("log")
        plt.xlabel(p)
        plt.ylabel(f"P({p})")
        plt.title(f"{galaxy}\nposterior of {p}", fontsize=9)
        plt.grid(True)
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(os.path.join(results_folder, f"{galaxy}_posterior_{p}.png"))
        plt.close()

    top = [int(i) for i in order[:min(int(P.TOP_N), n_ranked)]]

    # ---- report: files
    out_file = os.path.join(results_folder, f"{galaxy}_best_model.txt")
    with open(out_file, "w") as f:
        f.write(f"# SINEFIA fit of {obs_file}\n")
        f.write(f"# lines used (N={len(used)}): {', '.join(used)}\n")
        for note in notes:
            f.write(f"# note: {note}\n")
        f.write(f"# {len(params)} models with these settings, {n_ranked} with >= "
                f"{P.MIN_LINES} of the observed lines\n")
        f.write(f"# systematic error: {P.SYSTEMATIC_ERROR}   free scaling: {P.FREE_SCALING}   "
                f"min lines: {P.MIN_LINES}\n")
        f.write(f"# best model folder: {folders[best]}\n")
        f.write(f"# min chi2: {chi2[best]:.4g}   reduced chi2 (chi2 / (N_lines - "
                f"{n_scale})): {chi2_red[best]:.4g}\n")
        f.write(f"# free parameters: {', '.join(varied)}{' + flux scale' if n_scale else ''}"
                f"   dof: {dof}   chi2/dof: "
                f"{chi2[best] / dof if dof > 0 else float('nan'):.4g}\n")
        f.write(f"# best-model flux scale: {scale[best]:.6g}\n")
        f.write(f"# posteriors over {int(complete.sum())} models that have all the used lines\n")
        f.write("parameter median p16 p84 best_model\n")
        for p, (med, lo, hi, mle) in summary.items():
            f.write(f"{p} {med:.6g} {lo:.6g} {hi:.6g} {mle:.6g}\n")
        f.write("#\n# best model, line by line (model = scale x CLOUDY flux;"
                " residual in units of the total error)\n")
        f.write("# line observed error_total model residual\n")
        for j, n in enumerate(names):
            if used_mask[j]:
                m = scale[best] * model[best, j]
                f.write(f"# {n} {obs[j]:.4g} {sigma[j]:.4g} {m:.4g} "
                        f"{(obs[j] - m) / sigma[j]:+.2f}\n")
        f.write(f"#\n# top {len(top)} models by reduced chi2\n")
        f.write("# rank chi2_red N_lines " + " ".join(param_names) + "\n")
        for r, i in enumerate(top, start=1):
            f.write(f"# {r} {chi2_red[i]:.4g} {nlines[i]} "
                    + " ".join(f"{params[i][p]:.6g}" for p in param_names) + "\n")

    # every model with its chi^2 and probability, for later plots
    table_file = os.path.join(results_folder, f"{galaxy}_all_models.txt")
    header = " ".join(param_names) + " N_lines scale chi2 chi2_red probability"
    np.savetxt(table_file, np.column_stack([keys_arr, nlines, scale, chi2, chi2_red, like]),
               header=header, fmt="%.6g")
    return params[best], float(chi2_red[best]), len(used), notes


def fit_observations():
    if P.FIT_LINELIST_FILENAME:
        fit_lines = read_line_list(resolve_path(P.FIT_LINELIST_FILENAME))
        allowed = {line_key(l) for l in fit_lines}
        lines_text = ", ".join(fit_lines)
    else:
        allowed = None
        lines_text = "every usable line in each observation file"
    obs_files = sorted(glob.glob(os.path.join(resolve_path(P.OBS_FOLDER), P.OBS_FILE_PATTERN)))
    if P.GALAXIES_TO_FIT:
        obs_files = [f for f in obs_files
                     if any(g in os.path.basename(f) for g in P.GALAXIES_TO_FIT)]
    if not obs_files:
        print(f"\n3) No observation files matching {P.OBS_FILE_PATTERN} in "
              f"{resolve_path(P.OBS_FOLDER)}")
        return
    grid = load_models()
    if grid is None:
        return
    if not grid[1]:
        print(f"\n3) The database has {grid[4]} {P.SOURCE_TYPE} models, but none made with the "
              "settings and grid parameters of this settings file.")
        return

    param_names = grid[0]
    bar = Progress(len(obs_files), f"\n3) Fitting {len(obs_files)} "
                                   f"{'galaxy' if len(obs_files) == 1 else 'galaxies'} with "
                                   f"{len(grid[1])} models")
    rows = []
    for obs_file in obs_files:
        rows.append((galaxy_name(obs_file),) + fit_galaxy(obs_file, allowed, grid))
        bar.update(galaxy_name(obs_file))
    bar.close()

    # one line per galaxy: the best model (details in <galaxy>_best_model.txt)
    results_folder = resolve_path(P.RESULTS_FOLDER)
    os.makedirs(results_folder, exist_ok=True)
    summary_file = os.path.join(results_folder, "fit_summary.txt")
    skipped = [r for r in rows if r[1] is None]
    with open(summary_file, "w") as f:
        f.write(f"# SINEFIA fit summary: best model of every galaxy (grid {P.GRID_NAME})\n")
        f.write(f"# lines allowed in the fit: {lines_text}\n")
        f.write("# details, posteriors and plots: <galaxy>_best_model.txt and "
                "<galaxy>_posterior_<parameter>.png\n")
        f.write("galaxy N_lines chi2_red " + " ".join(param_names) + " N_notes\n")
        for galaxy, best, chi2_red, n_used, notes in rows:
            vals = (" ".join(f"{best[p]:.6g}" for p in param_names) if best is not None
                    else " ".join("nan" for _ in param_names))
            f.write(f"{galaxy} {n_used} {chi2_red:.4g} {vals} {len(notes)}\n")
        for galaxy, best, chi2_red, n_used, notes in rows:
            for note in notes:
                f.write(f"# {galaxy}: {note}\n")
    print(f"   {len(rows) - len(skipped)} of {len(rows)} "
          f"{'galaxy' if len(rows) == 1 else 'galaxies'} fitted.")
    for galaxy, *_, notes in skipped:
        print(f"   Not fitted: {galaxy} ({notes[-1].replace(': skipped', '')})")
    print(f"   The results are in the folder {os.path.relpath(results_folder, BASE_PATH)}/ :")
    print("     fit_summary.txt              the best model of every galaxy, one line each")
    print("     <galaxy>_best_model.txt      the parameters (median and 68% range), the best")
    print("                                  model line by line, the top models and the notes")
    print("     <galaxy>_posterior_<par>.png the posterior of each parameter")
    print("     <galaxy>_all_models.txt      chi^2 and probability of every model")


#######################################################
#  MAIN
#######################################################

def main(argv):
    check_only = "--check" in argv
    files = [a for a in argv if not a.startswith("--")]
    default = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sinefia_parameters.py")
    load_parameters(files[0] if files else default)

    if check_only:
        validate_parameters()
        grid_summary()
        first = all_combinations()[0]
        print("\nExample CLOUDY input for the first model")
        print(f"(folder: {os.path.join(P.MODELS_FOLDER, P.GRID_NAME, model_folder_name(first))})\n")
        print("\n".join(build_input_deck(first)))
        return

    validate_parameters()
    source = P.AGN_SED if P.SOURCE_TYPE == "AGN" else P.SED_TYPE
    print(f"SINEFIA: grid {P.GRID_NAME} ({P.SOURCE_TYPE}, {source}), "
          f"{len(all_combinations())} models   [python sinefia.py --check shows the details]")
    if P.RUN_CLOUDY:
        run_models()
    if P.EXTRACT_LINES:
        update_lines_database()
    if P.FIT_OBSERVATIONS:
        fit_observations()


if __name__ == "__main__":
    main(sys.argv[1:])
