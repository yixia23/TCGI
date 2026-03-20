"""
run_pi.py
=========
Calculates tropical cyclone potential intensity (PI).
Modified from the pyPI script by Daniel Gilford, PhD (daniel.gilford@rutgers.edu).
 
Usage
-----
Edit the USER SETTINGS section below, then run:
    python run_pi.py
 
Requirements
------------
- ts_{year}.nc   : monthly sea surface temperature (K or degC)
- psl_{year}.nc  : monthly sea level pressure (Pa or hPa)
- ta_{year}.nc   : monthly multi-level air temperature (K or °C, pressure levels in Pa or hPa)
- hus_{year}.nc  : monthly multi-level specific humidity (kg/kg or g/kg, pressure levels in Pa or hPa)
- Pressure levels must be named consistently across T and Q files
"""

import os
import sys
import time
import numpy as np
import xarray as xr
 
from tcpyPI.pi import pi
from tcpyPI.utilities import *


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
 
# File name templates. Use 'YYYY' as a placeholder for the year.
ts_filename  = "ts_YYYY.nc"
psl_filename = "psl_YYYY.nc"
ta_filename  = "ta_YYYY.nc"
hus_filename = "hus_YYYY.nc"
 
# Variable names inside each NetCDF file
ts_varname  = "ts"    # sea surface temperature
psl_varname = "psl"   # sea level pressure
ta_varname  = "ta"    # multi-level air temperature
hus_varname = "hus"   # multi-level specific humidity
 
# Name of the pressure level coordinate in the temperature and humidity files
plev_name = "plev"


# ====================================================================
# HELPER FUNCTION
# ====================================================================

def run_PI(ds, dim='p', CKCD=0.9):
    """ This function calculates PI over the given dataset using xarray """
    
    # calculate PI over the whole data set using the xarray universal function
    result = xr.apply_ufunc(
        pi,
        ds['sst'], ds['psl'], ds[dim], ds['t'], ds['q'],
        kwargs=dict(CKCD=CKCD, ascent_flag=0, diss_flag=1.0, V_reduc=0.8, ptop=30, miss_handle=0),
        input_core_dims=[[], [], [dim, ], [dim, ], [dim, ],],
        output_core_dims=[[], [], [], [], []],
        vectorize=True
    )

    vmax, pmin, ifl, t0, otl = result
    vmax = vmax.fillna(0)
    
    out_ds = xr.Dataset({'vmax': vmax, 'pmin': pmin, 'ifl': ifl, 't0': t0, 'otl': otl})
    
    out_ds.vmax.attrs['standard_name'],out_ds.vmax.attrs['units']='Maximum Potential Intensity','m/s'
    out_ds.pmin.attrs['standard_name'],out_ds.pmin.attrs['units']='Minimum Central Pressure','hPa'
    out_ds.ifl.attrs['standard_name']='pyPI Flag'
    out_ds.t0.attrs['standard_name'],out_ds.t0.attrs['units']='Outflow Temperature','K'
    out_ds.otl.attrs['standard_name'],out_ds.otl.attrs['units']='Outflow Temperature Level','hPa'

    return out_ds  


# ====================================================================
# MAIN LOOP
# ====================================================================

if __name__ == "__main__":

    os.makedirs(output_dir, exist_ok=True)
    
    for iyear in range(cal_start, cal_end + 1):
        
        # --- Skip year if output file already exists ---
        pi_out = f"{output_dir}/PI_{iyear}.nc"
        if os.path.exists(pi_out):
            print(f"[{iyear}] Output file already exists, skipping.")
            continue

        # --- Check that all input files are present before loading ---
        ts_path  = f"{input_dir}/{ts_filename.replace('YYYY', str(iyear))}"
        psl_path = f"{input_dir}/{psl_filename.replace('YYYY', str(iyear))}"
        ta_path  = f"{input_dir}/{ta_filename.replace('YYYY', str(iyear))}"
        hus_path = f"{input_dir}/{hus_filename.replace('YYYY', str(iyear))}"
        
        missing = [p for p in [ts_path, psl_path, ta_path, hus_path] if not os.path.exists(p)]
        if missing:
            print(f"[{iyear}] ERROR: Missing input file(s):")
            for p in missing:
                print(f"  {p}")
            sys.exit(1)

        # --- Load input fields ---
        TS_ds  = xr.open_dataset(ts_path)[ts_varname]
        PSL_ds = xr.open_dataset(psl_path)[psl_varname]
        T_ds   = xr.open_dataset(ta_path)[ta_varname]
        Q_ds   = xr.open_dataset(hus_path)[hus_varname]
        
        # --- Rename latitude/longitude to lat/lon if necessary ---
        rename_map = {}
        for ds_ref in [TS_ds, PSL_ds, T_ds, Q_ds]:
            if 'latitude' in ds_ref.dims and 'lat' not in ds_ref.dims:
                rename_map['latitude'] = 'lat'
                rename_map['longitude'] = 'lon'
        if rename_map:
            TS_ds  = TS_ds.rename(rename_map)
            PSL_ds = PSL_ds.rename(rename_map)
            T_ds   = T_ds.rename(rename_map)
            Q_ds   = Q_ds.rename(rename_map)
 
        # --- Ensure pressure levels are sorted from surface to top (descending) ---
        T_ds = T_ds.transpose("time", plev_name, "lat", "lon").sortby(plev_name, ascending=False)
        Q_ds = Q_ds.transpose("time", plev_name, "lat", "lon").sortby(plev_name, ascending=False)

        # --- Unit conversions (auto-detected from data values) ---
        # SST: K → °C
        ts  = TS_ds - 273.15 if float(TS_ds.isel(time=0).mean(skipna=True)) > 100 else TS_ds
        # Sea level pressure: Pa → hPa
        psl = PSL_ds / 100 if float(PSL_ds.isel(time=0).mean(skipna=True)) > 9000 else PSL_ds
        # Air temperature: K → °C; pressure coordinate: Pa → hPa
        t = T_ds - 273.15 if float(T_ds.isel(time=0, **{plev_name: 0}).mean(skipna=True)) > 100 else T_ds
        if float(t[plev_name][0]) > 9000:
            t = t.assign_coords({plev_name: t[plev_name].data / 100})
        # Specific humidity: kg/kg → mixing ratio in g/kg
        q  = Q_ds / (1 - Q_ds) * 1000 if float(Q_ds.isel(time=0, **{plev_name: 0}).mean(skipna=True)) < 0.1 \
             else Q_ds / (1 - Q_ds / 1000)
        if float(q[plev_name][0]) > 9000:
            q = q.assign_coords({plev_name: q[plev_name].data / 100})

        # --- Assemble input dataset for pyPI ---
        ds = xr.Dataset(
            data_vars={
                'sst': (['time', 'latitude', 'longitude'],          ts.data),
                'psl': (['time', 'latitude', 'longitude'],          psl.data),
                't':   (['time', 'p', 'latitude', 'longitude'],     t.data),
                'q':   (['time', 'p', 'latitude', 'longitude'],     q.data),
            },
            coords={
                'p':         ('p',         t[plev_name].data),
                'latitude':  ('latitude',  ts.lat.data),
                'longitude': ('longitude', ts.lon.data),
                'time':      ('time',      ts.time.data),
            }
        )
        

        # --- Run PI computation ---
        print(f"[{iyear}] Computing potential intensity...")
        t0 = time.time()
        ds = ds.load()
        result = run_PI(ds)
 
        # --- Save output ---
        encoding_settings = {v: {'zlib': True, 'complevel': 5} for v in result.data_vars}
        result = result.rename({'latitude': 'lat', 'longitude': 'lon'})
        result.to_netcdf(pi_out, encoding=encoding_settings, mode='w')
 
        print(f"[{iyear}] PI saved ({time.time() - t0:.1f}s)")
