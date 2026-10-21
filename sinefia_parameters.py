#######################################################
#######################################################
###                                                 ###
###          SINEFIA  -  USER PARAMETER BLOCK       ###
###                                                 ###
###   Define here EVERY parameter of the CLOUDY     ###
###   grid you want to run, and how the models are  ###
###   compared with your observations. Read the     ###
###   comments: each one tells you what the         ###
###   parameter does, its CLOUDY command, typical   ###
###   values, and units. You only edit this file -- ###
###   sinefia.py reads it automatically. There is   ###
###   no code in here, only settings.               ###
###                                                 ###
###   Two kinds of ionizing source can be modelled, ###
###   chosen with SOURCE_TYPE (SECTION 2):          ###
###     "SB"  -> stellar populations (HII regions / ###
###              starbursts), SED from an SSP grid  ###
###     "AGN" -> an AGN continuum (CLOUDY 'AGN')    ###
###                                                 ###
#######################################################
#######################################################
#
# -------------------- HOW THIS WORKS --------------------
#
#  * Every grid parameter is defined by THREE numbers:
#        <name>_init  -> first value of the grid (log unless noted)
#        <name>_end   -> last value of the grid
#        <name>_step  -> increment (in the same units)
#    and, for the optional ones, an ON/OFF flag in PARAM_SWITCHES.
#
#  * To keep a parameter FIXED (not gridded) set:  _step = 0
#    -> the code then uses only <name>_init.
#
#  * A few parameters whose values usually span decades (the stellar
#    ages, the AGN temperature, the sub-mm break of 'table AGN') are
#    given as a LIST of the values you want instead, e.g. [1e6, 3e6].
#    A one-element list keeps them fixed.
#
#  * To IGNORE an optional parameter entirely (use CLOUDY's default
#    and not even write the command) set its switch to False in
#    PARAM_SWITCHES. This keeps input files clean and avoids
#    forcing CLOUDY defaults you didn't ask for.
#
#  * Total number of models = product of the lengths of every
#    grid axis. Watch the combinatorics! Run
#        python sinefia.py --check
#    to check these settings and print the grid, the number of models
#    and an example CLOUDY input file BEFORE launching anything.
#    (python sinefia.py on its own runs the stages of SECTION 0.)
#
#  * Each model runs as its own CLOUDY job, in its own folder:
#        <MODELS_FOLDER>/<GRID_NAME>/<one folder per model>/
#    The folder name lists the parameter values of that model.
#
#  * Paths written as relative paths (e.g. "PopStar/...") are taken
#    relative to the folder that contains this file.
#
#  * The file is split into numbered sections (0, 0b, 1, 2a ... 12).
#    "SECTION n" in a comment points to that section of this file.
#
# --------------------------------------------------------


# =======================================================
# 0)  WHAT TO RUN
# =======================================================
# The pipeline (sinefia.py) has three stages. Switch each one on (True)
# or off (False):
#   RUN_CLOUDY      -> write the CLOUDY input files and run the models
#   EXTRACT_LINES   -> read the line fluxes of every finished model and
#                      store them in the HDF5 line database
#                      (DATABASE_FILE, SECTION 12)
#   FIT_OBSERVATIONS-> compare the database with your observed line
#                      fluxes: best models by chi^2 and Bayesian
#                      posteriors of every parameter (SECTION 12)

RUN_CLOUDY       = False
EXTRACT_LINES    = True
FIT_OBSERVATIONS = False

# GRID_NAME: a name you choose for this set of models (any short text,
# no spaces). It is simply the name of the folder where the CLOUDY runs
# of this grid are written:
#     <MODELS_FOLDER>/<GRID_NAME>/<one folder per model>/
# Why it is useful:
#   * Each model folder is named only by its grid values (U, n_H, age
#     ...). If you change a FIXED setting (say the density law or the
#     stellar template) the folder names would stay the same and new
#     runs would land on top of the old ones. Giving the new setup its
#     own GRID_NAME keeps both sets of runs, side by side:
#         models/SB_popstar_2.5Z/...
#         models/SB_popstar_2.5Z_constant_density/...
#   * It lets you keep several grids (SB, AGN, other metallicities) in
#     one MODELS_FOLDER and re-extract or extend any of them later.
#   * If you forget, the code notices: it stops when this folder already
#     holds finished models made with different settings.
# Re-running with the SAME name and the same settings is safe: finished
# models are skipped, so you can extend a grid or resume a crashed run.
GRID_NAME = "SB_popstar_2.5Z"



# =======================================================
# 0b) MASTER SWITCHES for the OPTIONAL parameters
# =======================================================
# ("SECTION n" everywhere in this file means the numbered section of
# THIS file where that parameter is set, e.g. SECTION 5b below.)
#
# The ionizing source (SECTION 2), the brightness
# (SECTION 3) and the gas density (SECTION 4) are ALWAYS used: CLOUDY
# needs them, so they have no switch. So are the stellar age (SECTION 6)
# for SOURCE_TYPE = "SB" and the AGN continuum parameters (SECTION 2b)
# for "AGN".
#
# True  -> the command is written to the CLOUDY input
#          (gridded if its _step != 0, fixed if _step == 0)
# False -> the command is omitted and CLOUDY uses its own default.

PARAM_SWITCHES = {
    "metallicity":    True,   # gas abundances + metallicity scaling (SECTION 5)
    "grains":         True,   # ISM / PAH grains                      (SECTION 5b)
    "dust_to_metals": False,  # vary grains independently of metals  (SECTION 5c)
    "cosmic_rays":    True,   # cosmic-ray background / rate         (SECTION 7)
    "covering_factor":True,   # geometric covering factor            (SECTION 8a)
    "radius":         False,  # inner radius of the cloud            (SECTION 8b)
    "filling_factor": False,  # clumping / filling factor            (SECTION 8d)
    "turbulence":     False,  # microturbulence velocity             (SECTION 9)
    "magnetic_field": False,  # magnetic field strength              (SECTION 9)
}


# =======================================================
# 1)  RUN / ENVIRONMENT SETTINGS
# =======================================================

# Number of CLOUDY models run at the same time (one CPU each).
Ncpus = 4

# Path to the cloudy executable. A path with a "/" is taken relative to
# this file (or give the full path); a bare name such as "cloudy.exe"
# is looked up on your system PATH.
cloudy_executable = "cloudy.exe"

# Where the model folders are created (relative to this file).
MODELS_FOLDER = "models"

# If a model folder already contains a finished run ("Cloudy exited
# OK" and a line output), it is skipped. Set True to run it again.
OVERWRITE_EXISTING = False

# Only used to print a rough time estimate before launching.
AVERAGE_MINUTES_PER_MODEL = 5


# =======================================================
# 2)  IONIZING SOURCE
# =======================================================
# Which kind of source lights the cloud:
#   "SB"  -> a stellar population: SED from a population-synthesis grid
#            (SECTION 2a), with a stellar AGE axis (SECTION 6) and a
#            stellar metallicity (SECTION 5)
#   "AGN" -> an AGN continuum (SECTION 2b): CLOUDY's 'AGN' command, with
#            its temperature and three slopes as grid axes, or one of
#            CLOUDY's tabulated AGN SEDs. No age and no stellar
#            metallicity.
# Everything else (log U, density, gas composition, geometry, stops,
# output and the fit) is shared by both.

SOURCE_TYPE = "SB"

# -------------------------------------------------------
# 2a) STELLAR POPULATION SED  (SOURCE_TYPE = "SB";  Hazy 1, 'table star')
# -------------------------------------------------------
# Choose ONE stellar template for the incident continuum. Set SED_TYPE
# (just after SED_LIBRARY below) to one of the keys of SED_LIBRARY. To use another stellar
# population (or a grid of single-star atmospheres), copy an entry,
# give it a new name and edit it: nothing else needs to change.
#
# These are grids compiled into CLOUDY with the 'compile stars' command
# (Hazy 1, 6.13). 'table star available' in CLOUDY lists the installed
# grids and their parameters. A template whose command contains {age}
# uses the stellar AGE axis (SECTION 6); the ages must lie inside the
# template grid.
#
# Each entry has:
#   "command"  : the CLOUDY line. {sed_file}, {age} and {Zstar} are
#                filled in automatically for each model, and so is
#                every name listed in "axes".
#   "sed_file" : the template file. A path with a "/" is taken relative
#                to this file; a bare file name is looked up by CLOUDY
#                in its own data folder.
#   "z_mode"   : how the STELLAR metallicity is chosen:
#                  "parameter" -> the template has a metallicity axis
#                                 and the command passes Z = {Zstar}
#                  "file"      -> one template file per metallicity;
#                                 the file name contains {Z4}, which is
#                                 the absolute Z x 10^4 on 4 digits
#                                 (Z = 0.02 -> "0200", Z = 0.05 -> "0500")
#                  None        -> the template has no metallicity choice
#   "z_units"  : (only for "parameter") the units of the template's
#                metallicity axis: "log_abs" = log10(Z), "abs" = Z
#   "Zsun_template": absolute Z the template calls solar. Used to turn
#                the RELATIVE metallicity of SECTION 5 into the absolute
#                Z that the template expects.
#   "available_Z": (only for "file") the absolute Z values you have files
#                for. The code stops with a clear message if you ask
#                for one you don't have.
#   "axes"     : (optional) any OTHER parameters of the template that you
#                want to grid, as {name: values}. Values are a list, or
#                (init, end, step) like every other axis. Each name is
#                filled into the command and becomes a grid axis, e.g.
#                    "axes": {"teff": [3.5e4, 4.0e4, 4.5e4]}
#                Danger (Hazy 1, 6.13): never write "log g" in a command,
#                CLOUDY then reads the temperature as a log. Write log(g).

SED_LIBRARY = {

    # PopStar, a single grid with an age AND a metallicity axis.
    # e.g. sp-kro-z0200.ascii, whose header lists
    # two parameters, "Age" and "log(Z)", so Z is given as log10(Z).
    "popstar": {
        "command": 'table star "{sed_file}" age = {age:.5f} years Z = {Zstar}',
        "sed_file": "PopStar/sp-kro-z0200.ascii",
        "z_mode": "parameter",
        "z_units": "log_abs",
        "Zsun_template": 0.02,
    },

    # PopStar, one file per metallicity (sp-kro-z0040, -z0200, -z0500 ...).
    # Use this one if each of your compiled files has an age axis only.
    "popstar_one_file_per_Z": {
        "command": 'table star "{sed_file}" age = {age:.5f} years',
        "sed_file": "PopStar/sp-kro-z{Z4}.ascii",
        "z_mode": "file",
        "Zsun_template": 0.02,
        "available_Z": [0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05],
    },

    # BPASS binary stellar populations.
    "bpass": {
        "command": 'table star "{sed_file}" age = {age:.5f} years',
        "sed_file": "BPASS/bpass_v2.ascii",
        "z_mode": None,
        "Zsun_template": 0.02,
    },

    # Starburst99 single-population output.
    "starburst99": {
        "command": 'table star "{sed_file}" age = {age:.5f} years',
        "sed_file": "SB99/sb99.ascii",
        "z_mode": None,
        "Zsun_template": 0.02,
    },

    # EXAMPLE of a template without ages: a grid of single-star
    # atmospheres with two parameters, Teff and log(g), both gridded
    # through "axes". (The file name is a placeholder: use your own.)
    "star_atmospheres": {
        "command": 'table star "{sed_file}" Teff = {teff:.6g} K, log(g) = {logg}',
        "sed_file": "Stars/my_atmospheres.ascii",
        "z_mode": None,
        "Zsun_template": 0.02,
        "axes": {"teff": [3.5e4, 4.0e4, 4.5e4], "logg": (4.0, 4.0, 0.0)},
    },
}

# Which entry of SED_LIBRARY above to use for this grid:
SED_TYPE = "popstar"   # <-- pick a key from SED_LIBRARY

# -------------------------------------------------------
# 2b) AGN CONTINUUM  (SOURCE_TYPE = "AGN";  Hazy 1, chapter 6)
# -------------------------------------------------------
# Choose the SHAPE of the AGN continuum with AGN_SED. Its brightness is
# set in SECTION 3 (e.g. log U), as for the stars. Every number of the chosen
# shape is a grid axis; only the parameters of the chosen shape are
# used, the others are ignored.
#
#   "AGN"             -> CLOUDY's parametric AGN continuum (Hazy 1, 6.2).
#                        The most flexible choice: a temperature and
#                        three slopes, all griddable.
#   "table AGN"       -> the fixed Mathews & Ferland (1987) radio-quiet
#                        AGN SED (Hazy 1, 6.12.2). One parameter: the
#                        wavelength of its sub-millimetre break.
#   "table power law" -> a power law that is well behaved at low and high
#                        energies (Hazy 1, 6.12.9). One parameter: the
#                        slope between the infrared break and 50 keV.
#   "table SED"       -> a tabulated SED file (Hazy 1, 6.11), e.g. one of
#                        the observed AGN SEDs that come with CLOUDY in
#                        data/SED/: "NGC5548.sed", "Mrk509.sed",
#                        "akn120.sed", or your own file. No parameters.
#
# CLOUDY never uses an implicit minus sign: every slope below is
# negative for a typical AGN.

AGN_SED = "AGN"

# --- "AGN":   AGN T=<agnT> K, a(ox)=<aox>, a(uv)=<auv>, a(x)=<ax> ---
# The continuum is a "Big Bump" (a power law nu^auv with an exponential
# cut-off at the temperature agnT; it peaks near 1 Ryd) plus an X-ray
# power law nu^ax, scaled so that the optical (2500 A) to X-ray (2 keV)
# slope is aox (Hazy 1, eq. 6.2):
#     f_nu = nu^auv exp(-h nu/k agnT) exp(-k T_IR/h nu) + a nu^ax
#
#   agnT  temperature of the Big Bump [K, linear], as a LIST of the
#         values you want. Typical 1e5 - 2e6 K: a hotter bump gives a
#         harder ionizing continuum. Give T in K: CLOUDY would read a
#         number <= 10 as log T.
#   aox   optical-to-X-ray slope alpha_ox. Typical -1.8 to -1.0; -1.4 is
#         the classic value (Zamorani et al. 1981). A very hot bump also
#         emits X-rays, and the real alpha_ox then comes out more
#         negative than the number given here.
#   auv   low-energy slope of the Big Bump. CLOUDY default -0.5
#         (Elvis et al. 1994; Francis 1993).
#   ax    slope of the X-ray power law. CLOUDY default -1.
#
# The defaults below are  AGN T=1.5e5 K, a(ox)=-1.4, a(uv)=-0.5, a(x)=-1.
# The 'AGN kirk' continuum of Korista et al. (1997) is agnT = [1e6]
# with the same three slopes.

agnT = [1.5e5]                                   # K, linear (a list)
aox_init, aox_end, aox_step = -1.4, -1.4, 0.0
auv_init, auv_end, auv_step = -0.5, -0.5, 0.0
ax_init,  ax_end,  ax_step  = -1.0, -1.0, 0.0

# --- "table AGN":   table AGN break <agn_break> microns ---
#   agn_break  wavelength of the sub-millimetre break [microns, linear],
#              as a LIST. Longward of it the SED falls as nu^2.5
#              (self-absorbed synchrotron). CLOUDY's default is 10.
#              Allowed: 0.46 to 3e7 microns; 3e7 is the same as
#              'table AGN no break'. The break matters most for dense,
#              high-U clouds (free-free heating).
agn_break = [10.0]

# --- "table power law":   table power law slope <pl_slope>, <break> micron break ---
#   pl_slope  f_nu ~ nu^pl_slope between the infrared break and 50 keV
#             (CLOUDY default -1; AGN are often modelled with -1.2 to
#             -2). Below the break the SED rises as nu^2.5, above 50 keV
#             it falls as nu^-2.
#   PL_IR_BREAK_MICRONS  wavelength of the infrared break (fixed;
#             CLOUDY's default is 10 microns).
pl_slope_init, pl_slope_end, pl_slope_step = -1.4, -1.4, 0.0
PL_IR_BREAK_MICRONS = 10.0

# --- "table SED":   table SED "<AGN_SED_FILE>" ---
# A bare file name is looked up by CLOUDY in its data folder; a path
# with a "/" is taken relative to this file.
AGN_SED_FILE = "NGC5548.sed"


# -------------------------------------------------------
# 2c) COSMIC MICROWAVE BACKGROUND  (both source types;  Hazy 1, 6.6)
# -------------------------------------------------------
# Adds the CMB blackbody at redshift CMB_REDSHIFT on top of the source.
# Hazy 1 (6.2) recommends it for the "AGN" continuum: the Big Bump
# falls to zero at long wavelengths, so CLOUDY warns "continuum has
# zero intensity ... Problems???" but still
# computes the model. The CMB fills that gap. It slightly changes the
# physics, so keep the same choice for grids you want to compare.
ADD_CMB = False
CMB_REDSHIFT = 0.0   # used only if ADD_CMB is True


# =======================================================
# 3)  BRIGHTNESS OF THE SOURCE   (Hazy 1, chapter 5)
# =======================================================
# SECTION 2 sets the SHAPE of the ionizing SED; this sets how BRIGHT it
# is. CLOUDY needs exactly one such command. Always used.
#
# Pick the quantity with BRIGHTNESS_MODE. The three numbers below are
# its value (all are log10), gridded like any other axis. The grid
# axis, the folder names and the database use the "axis" name.
#
# BRIGHTNESS_MODE      CLOUDY command          meaning                         axis
# ---------------      --------------          -------                         ----
# "U"              ->  ionization parameter    log U = log[phi(H) / (n_H c)]   U
# "phi(H)"         ->  phi(H)                  log flux of H-ionizing photons  phiH
#                                              at the cloud [cm^-2 s^-1]
# "intensity"      ->  intensity               log 4 pi J in H-ionizing        intensity
#                                              radiation [erg cm^-2 s^-1]
# "xi"             ->  xi                      log xi = L(1-1000 Ryd)/(n_H r^2) xi
#                                              [erg cm s^-1] (X-ray / AGN work)
# "Q(H)"           ->  Q(H)                    log number of H-ionizing        QH
#                                              photons of the source [s^-1]
# "luminosity"     ->  luminosity              log luminosity in H-ionizing    lum
#                                              radiation [erg s^-1]
# "luminosity total"-> luminosity total        log bolometric luminosity       lumtot
#                                              [erg s^-1]
#
# The first four give the brightness AT THE CLOUD ("intensity case").
# "Q(H)", "luminosity" and "luminosity total" give the brightness OF THE
# SOURCE ("luminosity case"): CLOUDY then also needs the distance, so
# switch "radius" on in PARAM_SWITCHES (SECTION 8b). log U is the usual
# choice for nebular work: models with the same U and different
# densities are similar (Hazy 1, 5.8).

BRIGHTNESS_MODE = "U"

brightness_init = -4.5
brightness_end  = -3.5
brightness_step = 0.5


# =======================================================
# 4)  GAS DENSITY  /  PRESSURE  LAW  (Hazy 1)
# =======================================================
# Choose how the density structure of the cloud behaves with depth.
#   "hden_constant"        -> constant hydrogen density  (simplest)
#   "constant_pressure"    -> constant *total* pressure
#   "constant_gaspressure" -> constant gas (nkT) pressure
#   "powerlaw_density"     -> n_H(r) = n_0 (r/r0)^index
#                             (needs an inner radius: switch "radius"
#                             on in PARAM_SWITCHES, SECTION 8b)
#
# In every case the grid axis below sets the *starting* hydrogen
# density log n_H at the illuminated face. Always used.

DENSITY_LAW = "constant_pressure"

# --- log n_H : hydrogen density at the illuminated face [cm^-3] ---
nh_init = 2.0
nh_end  = 4.0
nh_step = 1.0

# extra parameter for the power-law density law (ignored otherwise)
HDEN_POWERLAW_INDEX = -1.0


# =======================================================
# 5)  CHEMICAL COMPOSITION  (Hazy 1)
# =======================================================
#
# There are TWO metallicities and they should not be confused:
#
#   * GAS metallicity (Zgas): the metallicity of the nebular gas being
#     illuminated. It scales every metal (and the grains) of the
#     abundance set through the CLOUDY command 'metals and grains'.
#
#   * STELLAR metallicity (Zstar): the metallicity of the stars that
#     produce the ionizing SED. It selects the stellar template.
#     (SOURCE_TYPE = "SB" only; an AGN has no stellar metallicity.)
#
# BOTH are given here in the SAME, RELATIVE units:
#
#       Z = log10( metallicity / reference )
#
#   0       -> the reference itself (solar-like)
#   0.39794 -> log10(2.5)  = 2.5 x solar
#  -0.30103 -> log10(0.5)  = half solar
#  -1       -> 1/10 solar
#
# For the gas the reference is ABUNDANCE_SET below; for the stars it
# is "Zsun_template" of the chosen SED (SECTION 2). The code converts
# the stellar value to whatever the template needs (for PopStar:
# 2.5 x 0.02 = Z 0.05 -> log10(0.05) = -1.301), so you never have to.
#
# LINK_STELLAR_TO_GAS controls whether the stars follow the gas:
#   True   -> Zstar = Zgas for every model (one axis)
#   False  -> Zstar fixed at Zstar_fixed below

LINK_STELLAR_TO_GAS = True

# --- 5a) GAS metallicity, log10(Z / reference)  (grid axis) ---
Zgas_init = 0.39794   # log10(2.5) -> 2.5 x the abundance set
Zgas_end  = 0.39794
Zgas_step = 0.0

# --- 5a') STELLAR metallicity, log10(Z / Zsun_template) ---
# Used ONLY when LINK_STELLAR_TO_GAS is False. It is stored with every
# model as the parameter "Zstar", so grids run with different stellar
# metallicities can share one database without being confused.
Zstar_fixed = 0.39794

# Base abundance set written before scaling.
#   "ISM", "HII region", "PN", "solar", "primordial", "starburst",
#   "Orion", "Crab", or your own file name in quotes, e.g. '"myabund.abn"'
ABUNDANCE_SET = "ISM"

# Optional: set individual elements by hand, e.g. to scale He, C and N
# with metallicity at low Z (Nicholls et al. 2017). Each entry writes
#     element abundance <element> <value>
# with value = log10(n_element / n_H), after the metallicity scaling.
# Leave empty {} to keep the abundance set as it is. Example:
#     ELEMENT_ABUNDANCES = {"helium": -1.01, "carbon": -3.9, "nitrogen": -4.6}
ELEMENT_ABUNDANCES = {}

# --- 5b) Grains ---
# Activated by PARAM_SWITCHES["grains"]. One 'grains <type>' command is
# written per entry. Common combo for ISM-like dust + PAHs:
GRAIN_TYPES = ["ISM", "PAH"]

# --- 5c) Independent dust-to-metals (grain) scaling ---
# Activated by PARAM_SWITCHES["dust_to_metals"]. log10 of an EXTRA factor
# applied to every grain type on top of the gas metallicity:
#   0 -> grains follow the metals exactly,  -0.3 -> half as much dust.
dust_init, dust_end, dust_step = 0.0, 0.0, 0.0


# =======================================================
# 6)  STELLAR POPULATION AGE  (SOURCE_TYPE = "SB" only)
# =======================================================
# List the ages (in YEARS, linear) you want. They must exist in (or lie
# inside) your stellar template grid. Ignored for an AGN, and for a
# template whose command has no {age}.
age = [1.0e6, 1.0e7, 19952623.1496888]   # 1, 10, ~20 Myr


# =======================================================
# 7)  COSMIC RAYS  (Hazy 1)
# =======================================================
# Activated by PARAM_SWITCHES["cosmic_rays"]. COSMIC_RAY_MODE:
#   "background" -> the mean Galactic cosmic-ray background, times
#                   COSMIC_RAY_BACKGROUND_FACTOR (linear; 1 = Galactic)
#   "rate"       -> a fixed H0 ionization rate, gridded below
#                   (log s^-1, typical -16 to -15)
COSMIC_RAY_MODE = "background"
COSMIC_RAY_BACKGROUND_FACTOR = 1.0   # LINEAR multiple of the Galactic background

# grid axis for the "rate" mode (log primary ionization rate per H)
cr_init, cr_end, cr_step = -16.0, -15.0, 1.0


# =======================================================
# 8)  GEOMETRY  (Hazy 1)
# =======================================================

# --- 8a) Covering factor (fraction of 4pi covered by gas, linear, 0-1) ---
# Note: when fitting observations with a free flux scale (SECTION 12)
# the covering factor cannot be constrained; it only rescales fluxes.
covfac_init = 0.8
covfac_end  = 0.8
covfac_step = 0.0

# --- 8b) Inner radius of the cloud, log10(r / cm) ---
# Activated by PARAM_SWITCHES["radius"]. Not needed for a plane-parallel
# cloud defined by log U, but REQUIRED for SPHERE = True, for the
# "powerlaw_density" law, and when the brightness is a luminosity
# ("Q(H)", "luminosity", "luminosity total" in SECTION 3).
radius_init, radius_end, radius_step = 19.0, 19.0, 0.0

# --- 8c) Sphere vs open geometry ---
# False -> plane-parallel / open (default). True -> closed sphere
# (needs the radius above).
SPHERE = False

# --- 8d) Filling / clumping factor, log10(ff), must be <= 0 ---
# Activated by PARAM_SWITCHES["filling_factor"]. 0 -> ff = 1 (no clumping).
ff_init, ff_end, ff_step = 0.0, 0.0, 0.0
FILLING_INDEX = 0.0   # radial power-law index of the filling factor


# =======================================================
# 9)  EXTRA MICROPHYSICS (optional axes)
# =======================================================

# --- Turbulence, km/s (linear). Typical HII regions: 2-10 km/s ---
turb_init, turb_end, turb_step = 5.0, 5.0, 0.0

# --- Magnetic field, log10(B / Gauss). Diffuse ISM ~ -5 (10 microG) ---
bfield_init, bfield_end, bfield_step = -5.0, -5.0, 0.0


# =======================================================
# 9b)  YOUR OWN GRID AXES  (any CLOUDY command)
# =======================================================
# Grid ANY other CLOUDY command without touching the code. Each entry
#     "name": {"values": ..., "command": "... {value} ..."}
# becomes a grid axis called "name": for every model, {value} in the
# command is replaced by that model's value, and the line is written to
# the CLOUDY input. The axis is stored in the database and fitted like
# the built-in ones.
#
#   "values"   : a list [v1, v2, ...] or (init, end, step); a single
#                value, [v1], writes the command with that fixed value.
#   "command"  : the CLOUDY line. Check its syntax and units in Hazy,
#                and which numbers CLOUDY reads as log or linear.
#   "log_plot" : (optional) True for a linear quantity spanning decades:
#                its posterior is then computed and plotted in log10.
#
# Use a name that is not already a grid axis (U, nh, Zgas, age, cov ...).
# Examples:
#   "He":    {"values": (-1.1, -0.9, 0.1),
#             "command": "element abundance helium {value}"},
#   "angle": {"values": [0, 30, 60],
#             "command": "illumination angle {value} deg"},
#   "Tdust": {"values": [20, 40], "log_plot": False,
#             "command": "constant grain temperature {value} K linear"},
CUSTOM_AXES = {
}


# =======================================================
# 10)  STOPPING CRITERIA AND ITERATIONS  (Hazy 1)
# =======================================================
# How CLOUDY decides the outer edge of the cloud. Pick ONE primary
# stop with STOP_MODE; you can list additional safety stops in
# EXTRA_STOPS.
#
# The stop VALUE is generic -- it is "the number that goes with
# whatever STOP_MODE you chose". Its meaning and units depend on the
# mode (see the table below). Set stop_step != 0 to grid it.
#
# STOP_MODE          meaning of the stop value                   units
# ----------         -------------------------                   -----
# "H2_column"    ->  H2 column density                           log cm^-2
# "column"       ->  total hydrogen column density               log cm^-2
# "neutral_col"  ->  neutral hydrogen column density             log cm^-2
# "Av"           ->  visual extinction A_V                       mag
# "temperature"  ->  stop when T drops below this value          K (linear)
# "thickness"    ->  physical thickness of the cloud             log cm
# "zone"         ->  number of zones computed                    integer
STOP_MODE = "H2_column"

stop_init = 20.0
stop_end  = 20.0
stop_step = 0.0

# Mapping of STOP_MODE -> CLOUDY command ({val} is the stop value).
# To add a new stopping criterion, add an entry here and document it in
# the table above -- nothing else needs to change.
# Tip: after a first test run, open the .out file and check, in the
# echo of the input at the top and in the "stop" reason at the end,
# that CLOUDY understood the stop the way you meant it.
STOP_COMMANDS = {
    "H2_column":   "stop H2 column density {val}",
    "column":      "stop column density {val}",
    "neutral_col": "stop neutral column density {val}",
    "Av":          "stop Av {val}",
    "temperature": "stop temperature {val} K linear",
    "thickness":   "stop thickness {val}",
    "zone":        "stop zone {val}",
}

# LOWEST TEMPERATURE. By default CLOUDY stops a model when the gas
# cools below 4000 K, i.e. right at the ionization front. A stop deeper
# in the cloud (an H2 column, a large A_V ...) is then never reached.
# Set a lower value (K, linear) to let the model continue into the
# neutral / molecular gas, "off" to remove the temperature limit
# completely ('stop temperature off'), or None to keep CLOUDY's
# 4000 K default.
LOWEST_TEMPERATURE = 10.0

# Iterations. Emergent line fluxes need at least a second iteration so
# CLOUDY knows the optical depths of the whole cloud.
#   "iterate"                -> two iterations (usual choice)
#   "iterate to convergence" -> iterate until the optical depths converge
#   None                     -> a single iteration (not recommended)
ITERATE = "iterate"

# Additional always-on stops (written verbatim), e.g.
#   "stop Av 100"
EXTRA_STOPS = []

# Any other CLOUDY commands, written verbatim at the end of the input,
# e.g. ["failures 3", "print last iteration"].
EXTRA_COMMANDS = []


# =======================================================
# 11)  OUTPUT (line list to save)
# =======================================================
# File listing the emission lines CLOUDY should report (one per line,
# CLOUDY c22 label + wavelength, e.g. "O  3 5006.84A"). Relative to this
# file; it is copied into every model folder.
list_lines_filename = "linelist.dat"
SAVE_LINE_OPTIONS = "emergent absolute"   # how fluxes are reported

# Also save CLOUDY's overview (T, n_e, ionization vs depth) for every
# model. Useful to check models later; costs some disk space.
SAVE_OVERVIEW = True


# =======================================================
# 12)  COMPARISON WITH OBSERVATIONS
# =======================================================
# The line database: ONE HDF5 file holding every extracted model.
# Each model is stored once, under a unique key (the md5 hash of its
# parameters and fixed settings), with its parameter values, its model
# folder and the fluxes of the lines it has. Grids run at different
# times (other GRID_NAMEs, other settings, SB and AGN) can all live in
# the same file: a rerun model is updated, nothing is duplicated, and
# the fit only uses the models made with the settings in THIS file.
# EXTRACT_LINES writes it and FIT_OBSERVATIONS reads it. Relative to
# this file.
DATABASE_FILE = "models/lines_db.h5"

# Your observed line fluxes: one text file per galaxy, three columns
#     line_name   flux   error
# where line_name is the CLOUDY label with the spaces replaced by "_"
# (e.g. "C_2_157.636m"). Lines starting with "#" are ignored.
OBS_FOLDER = "observations"
OBS_FILE_PATTERN = "*_combined_final.txt"

# Only fit these galaxies (names that appear in the file names).
# Leave empty [] to fit every file in OBS_FOLDER.
GALAXIES_TO_FIT = []

# Which lines are used in the fit (same format as the CLOUDY line list;
# commented lines are skipped). Relative to this file. Set None to use
# every usable line of each observation file.
# Line names are matched by species and wavelength, so "Ne 2 12.8101m",
# "Ne_2_12.8101m" and "Ne_2_12.81010m" are recognised as the same line.
# Observed lines with a nan or negative error are skipped (and reported).
FIT_LINELIST_FILENAME = None

# Extra fractional error added in quadrature to every observed flux
# (calibration + model uncertainties). 0.15 = 15 %. Use 0 to weight
# by the observational errors only.
SYSTEMATIC_ERROR = 0.15

# True  -> each model gets its own best flux scale factor
#          (observed = scale x model), so only line RATIOS matter.
# False -> the absolute model fluxes are compared directly.
FREE_SCALING = True

# A model is only ranked if it predicts at least this many of the
# observed lines (models may lack a line, e.g. if it was not in their
# line list).
MIN_LINES = 3

# How many of the best models (lowest reduced chi^2) are printed and
# written to the results file.
TOP_N = 5

# Where the results are written: for every galaxy the best models, the
# posterior (median and 68% interval) of every parameter with one plot
# each, and a table of every model with its chi^2 and probability.
RESULTS_FOLDER = "results"


#######################################################
###                                                 ###
###   END OF THE SETTINGS.  Run:                    ###
###       python sinefia.py --check   (check only)  ###
###       python sinefia.py           (run)         ###
###                                                 ###
#######################################################
