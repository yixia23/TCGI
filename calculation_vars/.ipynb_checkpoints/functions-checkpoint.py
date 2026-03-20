import metpy.calc as mpcalc
import xarray as xr
import numpy as np
import xesmf as xe
from metpy.units import units


def wrap_ds(input_ds, var_name, var_attrs):
    """
    Wrap xarray datasets with attributes.

    Parameters:
    - input_ds (xr.Dataset): The input xarray dataset.
    - var_name (str): The name of the new variable to be added.
    - var_attrs (dict): Attributes to be added to the new variable.
      e,g. new_variable_attrs = {'description': 'New variable data', 'units': 'meters'}
    """
    # Create a new dataset and add the variable
    new_ds = xr.Dataset({f'{var_name}': input_ds.copy()})
    # Add attributes to the variable in the dataset 
    new_ds[f'{var_name}'].attrs.update(var_attrs)

    return new_ds


def get_dim_index(ds):
    # Find the index of the "lat" dimension
    if 'lat' in ds.dims:
        lat_index = ds.dims.index('lat')
    elif 'latitude' in ds.dims:
        lat_index = ds.dims.index('latitude')
    else:
        print('Neither "lat" nor "latitude" dimension found in the dataset.')
        return

    if 'lon' in ds.dims:
        lon_index = ds.dims.index('lon')
    elif 'longitude' in ds.dims:
        lon_index = ds.dims.index('longitude')
    else:
        print('Neither "lon" nor "longitude" dimension found in the dataset.')
        return

    return [lon_index, lat_index]


def relative_vort(u_wind, v_wind):
    # calculate relative vorticity [s-1]
    x_dim_index, y_dim_index = get_dim_index(u_wind)  # (x_dim:lon, y_dim:lat)
    u_wind = u_wind * units('m/s')

    v_wind = v_wind * units('m/s')
    relative_vort = mpcalc.vorticity(u_wind, v_wind, x_dim=x_dim_index, y_dim=y_dim_index)
    relative_vort_ds = wrap_ds(relative_vort.metpy.dequantify(), 'rvort', {'long_name': 'relative vorticity'})
    return relative_vort_ds


def absolute_vort(u_wind, v_wind, level_dim='level'):
    # calculate absolute vorticity [s-1]
    if u_wind.sortby(level_dim, ascending=True)[level_dim][-1].item() < 9000:
        u_wind = u_wind.assign_coords(**{level_dim: u_wind[level_dim].data * 100})  # change unit to Pa if nessesary
        v_wind = v_wind.assign_coords(**{level_dim: v_wind[level_dim].data * 100})
    x_dim_index, y_dim_index = get_dim_index(u_wind)  # (x_dim:lon, y_dim:lat)
    u_wind = u_wind * units('m/s')
    v_wind = v_wind * units('m/s')
    absolute_vort = mpcalc.absolute_vorticity(u_wind.sel({level_dim: 85000}), v_wind.sel({level_dim: 85000}), x_dim=x_dim_index,
                                              y_dim=y_dim_index)  # (x_dim:lon, y_dim:lat)
    absolute_vort_ds = wrap_ds(absolute_vort.metpy.dequantify(), 'avort', {'long_name': 'absolute vorticity'})
    return absolute_vort_ds


def wind_shear(u_wind, v_wind, level_dim='level', lower_level=20000, upper_level=85000):
    # calculate wind shear [m s-1]
    # The magnitude of the wind shear between two layers
    if u_wind.sortby(level_dim, ascending=True)[level_dim][-1].item() < 9000:
        u_wind = u_wind.assign_coords(**{level_dim: u_wind[level_dim].data * 100})  # change unit to Pa if nessesary
        v_wind = v_wind.assign_coords(**{level_dim: v_wind[level_dim].data * 100})
    upper_wind = u_wind.sel({f'{level_dim}': upper_level}) + 1j * v_wind.sel({f'{level_dim}': upper_level})
    lower_wind = u_wind.sel({f'{level_dim}': lower_level}) + 1j * v_wind.sel({f'{level_dim}': lower_level})
    wind_shear = np.abs(upper_wind - lower_wind)
    wind_shear_ds = wrap_ds(wind_shear, 'ws', {'unit': 'm s-1',
                                               'long_name': f'magnitude of the wind shear between {upper_level}- and {lower_level}-Pa'})
    return wind_shear_ds


def saturation_vapor_pressure(ta):
    '''
    Calculate saturation water vapor pressure (es) based on formula given in Bolton (1980). If T is replaced by the dew point Td, this function will return the water vapor pressure (e).

    Parameters:
    - T: temperature at multi-layers [C]
    - es: saturation water vapor pressure [hPa]
    '''
    if ta.sortby(level_dim, ascending=False)[0, -1, 0, 0] > 100:
        ta = ta - 273.15  # change unit to Calcius if nessesary
    es = 6.112 * np.exp(17.67 * ta / (ta + 243.5))
    es_ds = wrap_ds(es, 'es', {'unit': 'hPa', 'long_name': 'saturation vapor pressure'})
    return es_ds


def saturation_specific_humidity(ta, level_dim='level'):
    '''
    Calculate saturation specific humidity (qs) based on formula given in Bolton (1980). If T is replaced by the dew point Td, this function will return the specific humidity (q).

    Parameters:
    - T: temperature at multi-layers [C]
    - p: pressure at multi-layers [Pa]
    - es: saturation water vapor pressure [Pa]
    - qs: saturation specific humidity [g kg-1]
    '''
    if ta.sortby(level_dim, ascending=False)[0, -1, 0, 0] > 100:
        ta = ta - 273.15  # change unit to Calcius if nessesary

    es, qs = [np.zeros_like(ta) for _ in range(2)]
    eps = 0.62197
    for i in range(len(ta[level_dim])):
        p = ta[level_dim][i].data
        T = ta[:, i]
        es[:, i] = 611.2 * np.exp(17.67 * T / (T + 243.5))
        qs[:, i] = (eps * es[:, i]) / (p - (1 - eps) * es[:, i])
    qs_ds = wrap_ds(xr.DataArray(qs, dims=ta.dims, coords=ta.coords), 'qs',
                    {'unit': 'kg kg-1', 'long_name': 'saturation specific humidity'})
    return qs_ds


def saturated_water_vapor(ta, psl, ts, level_dim='level'):
    '''
    Calculate saturated water vapor (swv [kg m-2]) 

    Input datasets:
    - ta: temperature at multi-layers [C/K]
    - psl: mean sea level pressure [Pa/hPa]
    - ts: skin temperature (SST) [C/K]
    '''
    ta = ta.sortby(level_dim, ascending=True)
    psl = psl
    ts = ts

    # adjust units
    if ta.sortby(level_dim, ascending=True)[level_dim][-1].item() < 9000:
        ta = ta.assign_coords(**{level_dim: ta[level_dim].data * 100})  # change unit to Pa if nessesary
    if ta.sortby(level_dim, ascending=True)[0, -6, 0, 0] > 100:
        ta = ta - 273.15  # change unit to Calcius if nessesary
    if psl[0, 0, 0] < 9000:
        psl = psl * 100  # change unit to Pa if nessesary
    if ts[0, 0, 0] > 100:
        ts = ts - 273.15  # change unit to Calcius if nessesary

    # Setting constants
    pubc = 15000  # Upper BC pressure at which we take qs = 0.
    pmin = 10000    # Original script pmin=10000, modified for CMIP6 calculation, set pmin=20000
    pmax = 100000  # Min & Max pressure level at which qs is computed [Pa]
    g = 9.81
    eps = 0.62197

    # Calculate saturation specific humidity
    qs = saturation_specific_humidity(ta, level_dim)
    qs = qs.qs.sel({f'{level_dim}': slice(pmin, pmax)})

    qsold, swv = [np.zeros_like(psl) for _ in range(2)]
    for ilevel in range(len(qs[level_dim].sel({f'{level_dim}': slice(pubc, pmax)})) - 1):
        ilevel = ilevel + np.where(qs[level_dim].values == pubc)[0].item()
        if ilevel == np.where(qs[level_dim].values == pubc)[0].item():
            '''
            Add wvp contribution from layer lev-0.5, 
            where we take qs=0 (can't use log approx)
            '''
            dpres = qs[level_dim].isel({f'{level_dim}': ilevel}).values - qs[level_dim].isel(
                {f'{level_dim}': (ilevel - 1)}).values
            swv = 0.5 * qs.sel({f'{level_dim}': pubc}) * dpres / g
        else:
            '''
            If we assume log(qs) is a linear function of p
            wvp contribution is qs*(dp/g)*f(x), f(x) = -x/log(1-x) 
            Use a Taylor series approx to f(x) to avoid roundoff
            problems near x = 0. The standard trapezoidal rule is equivalent
            to only keeping the linear term in the Taylor series. This
            typically leads to overestimates of 0.5-1.5 mm in swv.
            '''
            x = 1 - (qsold / qs.isel({f'{level_dim}': ilevel}))
            tay = 1 - (x / 2) - (x ** 2 / 12) - (x ** 3 / 24)
            dpres = qs[level_dim].isel({f'{level_dim}': ilevel}).values - qs[level_dim].isel(
                {f'{level_dim}': (ilevel - 1)}).values
            vv = qs.isel({f'{level_dim}': ilevel}) * tay * dpres / g
            vvz = np.nan_to_num(vv)  # Replace NaN values with 0 in vv
            swv = swv + vvz
        qsold = qs.isel({f'{level_dim}': ilevel})

    '''
    Add surface layer contribution. It is assumed that psl > 1000 hPa; 
    there can be nontrivial errors above elevated land surfaces.
    '''
    psl = psl.where(~(psl.isnull()).data, 100000)

    es_surface = 611.2 * np.exp(17.67 * ts / (ts + 243.5))
    qs_surface = (eps * es_surface) / (psl - (1 - eps) * es_surface)
    x = 1 - (qsold / qs_surface)
    tay = 1 - (x / 2) - (x ** 2 / 12) - (x ** 3 / 24)
    dpres = psl - qs[level_dim].isel({f'{level_dim}': ilevel}).values
    vv = qs_surface * tay * dpres / g
    swv = swv + vv
    swv_ds = wrap_ds(swv, 'swv', {'unit': 'kg m-2', 'long_name': 'saturated water vapor'})
    return swv_ds
