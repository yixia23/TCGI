"""
calculate_avort_ws.py
=====================
Calculates absolute vorticity (avort) and vertical wind shear (ws).
 
Usage
-----
Edit the USER SETTINGS section below, then run:
    python calculate_avort_ws.py
 
Requirements
------------
- ua_{year}.nc : monthly multi-level zonal wind (must include 200 and 850 hPa)
- va_{year}.nc : monthly multi-level meridional wind (must include 200 and 850 hPa)
- Pressure levels must be in Pa (e.g. 85000 for 850 hPa)
- functions.py must be in the same directory
"""

import os
import sys
import time
import numpy as np
import xarray as xr
import functions


# ====================================================================
# USER SETTINGS — modify these for your dataset
# ====================================================================
 
# First and last year to process
cal_start, cal_end = 1951, 2100
 
# Directory containing the input wind files (ua_{year}.nc, va_{year}.nc)
input_dir = "INPUT_DIRECTORY"
 
# Directory where output files will be saved.
# This corresponds to WORK_DIR in the TCGI calculation notebooks.
output_dir = "OUTPUT_DIRECTORY"

# File name templates for the wind input files. Use 'YYYY' as a placeholder
# for the year, e.g. 'Ua_YYYY.nc' will be resolved to 'Va_1951.nc', etc.
u_filename = "Ua_YYYY.nc"
v_filename = "Va_YYYY.nc"
 
# Variable names for zonal and meridional wind inside the NetCDF files
u_varname = "U"
v_varname = "V"
 
# Name of the vertical coordinate in the input files (must be in Pa)
plev_name = "plev"


# ====================================================================
# MAIN LOOP
# ====================================================================
 
os.makedirs(output_dir, exist_ok=True)
 
for iyear in range(cal_start, cal_end + 1):
 
    # --- Skip year if both output files already exist ---
    avort_out = f"{output_dir}/avort_{iyear}.nc"
    ws_out    = f"{output_dir}/ws_{iyear}.nc"
    if os.path.exists(avort_out) and os.path.exists(ws_out):
        print(f"[{iyear}] Output files already exist, skipping.")
        continue
 
    # --- Load input wind fields ---
    u_path = f"{input_dir}/{u_filename.replace('YYYY', str(iyear))}"
    v_path = f"{input_dir}/{v_filename.replace('YYYY', str(iyear))}"
 
    if not os.path.exists(u_path) or not os.path.exists(v_path):
        print(f"[{iyear}] ERROR: Missing input file(s). Expected:\n  {u_path}\n  {v_path}")
        sys.exit(1)
 
    U_ds = xr.open_dataset(u_path)[u_varname]
    V_ds = xr.open_dataset(v_path)[v_varname]
 
    # --- Absolute Vorticity at 850 hPa ---
    t0 = time.time() 
    result_avort = functions.absolute_vort(
        U_ds.sel({plev_name: 85000}),
        V_ds.sel({plev_name: 85000})
    )
    result_avort.rename('avort').to_netcdf(
        avort_out,
        encoding={'avort': {'zlib': True, 'complevel': 5}}, mode='w'
    ) 
    print(f"[{iyear}] avort saved ({time.time() - t0:.1f}s)")
 
    # --- Vertical Wind Shear (between 850 and 200 hPa) ---
    t0 = time.time() 
    result_ws = functions.wind_shear(U_ds, V_ds, plev_name)
    result_ws.to_netcdf(
        ws_out,
        encoding={'ws': {'zlib': True, 'complevel': 5}}, mode='w'
    )
    print(f"[{iyear}] ws saved ({time.time() - t0:.1f}s)\n")