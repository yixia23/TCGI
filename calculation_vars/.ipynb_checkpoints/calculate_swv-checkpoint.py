"""
calculate_swv.py
================
Calculates saturated water vapor (swv).
 
Usage
-----
Edit the USER SETTINGS section below, then run:
    python calculate_swv.py
 
Requirements
------------
- PSL_{year}.nc  : monthly surface pressure (Pa/hPa)
- TS_{year}.nc   : monthly surface temperature (K/degC)
- T_{year}.nc    : monthly multi-level temperature profile (K/degC)
- Units will be adjusted accordingly in the script.
- functions.py must be in the same directory
"""
 

import os
import sys
import time
import xarray as xr
import functions


# ====================================================================
# USER SETTINGS — modify these for your dataset
# ====================================================================
 
# First and last year to process
cal_start, cal_end = 1951, 2100
 
# Directory containing the input files
input_dir = "INPUT_DIRECTORY"
 
# Directory where output files will be saved.
# This corresponds to WORK_DIR in the TCGI calculation notebooks.
output_dir = "OUTPUT_DIRECTORY"
 
# File name templates for input files. Use 'YYYY' as a placeholder for the
# year, e.g. 'PSL_YYYY.nc' will be resolved to 'PSL_1951.nc', etc.
psl_filename = "PSL_YYYY.nc"
ts_filename  = "TS_YYYY.nc"
t_filename   = "T_YYYY.nc"
 
# Variable names inside each NetCDF file
psl_varname = "PSL"   # surface pressure
ts_varname  = "TS"    # surface temperature
t_varname   = "T"     # multi-level temperature profile
 
# Name of the vertical coordinate in the temperature file (must be in Pa)
plev_name = "plev"


# ====================================================================
# MAIN LOOP
# ====================================================================
 
os.makedirs(output_dir, exist_ok=True)
 
for iyear in range(cal_start, cal_end + 1):
 
    # --- Skip year if output file already exists ---
    swv_out = f"{output_dir}/swv_{iyear}.nc"
    if os.path.exists(swv_out):
        print(f"[{iyear}] Output file already exists, skipping.")
        continue
 
    # --- Check that all input files are present before loading ---
    psl_path = f"{input_dir}/{psl_filename.replace('YYYY', str(iyear))}"
    ts_path  = f"{input_dir}/{ts_filename.replace('YYYY', str(iyear))}"
    t_path   = f"{input_dir}/{t_filename.replace('YYYY', str(iyear))}"
 
    missing = [p for p in [psl_path, ts_path, t_path] if not os.path.exists(p)]
    if missing:
        print(f"[{iyear}] ERROR: Missing input file(s):")
        for p in missing:
            print(f"  {p}")
        sys.exit(1)
 
    # --- Load input fields ---
    PSL_ds = xr.open_dataset(psl_path)[psl_varname]
    TS_ds  = xr.open_dataset(ts_path)[ts_varname]
    T_ds   = xr.open_dataset(t_path)[t_varname]
 
    # --- Compute saturated water vapor ---
    t0 = time.time()
    result = functions.saturated_water_vapor(T_ds, PSL_ds, TS_ds, plev_name)
    result.rename('swv').to_netcdf(
        swv_out,
        encoding={'swv': {'zlib': True, 'complevel': 5}}, mode='w'
    )
 
    print(f"[{iyear}] swv saved ({time.time() - t0:.1f}s)")
