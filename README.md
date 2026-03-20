# TCGI: Tropical Cyclone Genesis Index

A Python-based workflow for computing the **Tropical Cyclone Genesis Index (TCGI)**, which estimates the probability of tropical cyclone genesis using historical statistical relationships between large-scale environmental conditions and observed TC genesis, derived by Poisson regression (Tippett et al. 2011). Two variants are supported in this repository — **TCGI-CRH** and **TCGI-SD**, with different choices of the environmental humidity parameter, see details in Camargo et al. (2014) and Lee et al. (2020).

---

## Table of Contents

- [Overview](#overview)
- [Package Structure](#package-structure)
- [Python Packages Requirements](#python-packages-requirements)
- [Input Variables](#input-variables)
- [Workflow](#workflow)
- [Output Files](#output-files)
- [References](#references)
- [License](#license)

---

## Overview

The TCGI expressed as a Poisson equation:

$$\mu = \exp\left(b_0 + b_\eta \cdot \eta + b_H \cdot H + b_T \cdot T + b_V \cdot V + \log(\cos\phi \cdot \Delta x \cdot \Delta y)\right)$$

where:

| Symbol | Variable | Units |
|---|---|---|
| $\eta$ | Absolute vorticity at 850 hPa | 10⁻⁵ s⁻¹ |
| $H$ | Moisture term: CRH (TCGI-CRH) or SD (TCGI-SD) | % or kg m⁻² |
| $T$ | Tropical Cyclone Potential Intensity | m s⁻¹ |
| $V$ | Vertical wind shear (between 850 and 200 hPa) | m s⁻¹ |
| $\phi$ | Latitude | deg |

A **mean-state bias correction** is applied to each model variable before TCGI calculation, removing mean state climatological biases relative to ERA5 over the 1981–2010 baseline period.

---

## Package Structure

```
.
├── README.md
│
├── calculation_vars/
│   ├── calculate_avort_ws.py     # Calculate absolute vorticity and wind shear
│   ├── calculate_swv.py          # Calculate saturated water vapor
│   ├── functions.py              # Calculation functions
│   └── tcpyPI/                   # Potential intensity package developed by Daniel Gilford (see ref below)
│       ├── __init__.py
│       ├── pi.py
│       ├── constants.py
│       ├── utilities.py
│       ├── run_PI.py             # Calculate tropical cyclone potential intensity
│       └── setup.py
│
│── ERA5_landmask_2deg.nc         # Land mask for ERA5 in 2° × 2° resolution
│
├── cal_TCGI_v2_part1.ipynb       # Compute model climatology (1981–2010)
├── cal_TCGI_v2_part2.ipynb       # Compute TCGI
│
├── clim_ERA5_2deg/               # ERA5 1981–2010 climatologies in 2° × 2° resolution
│   ├── avort_ERA5_clim_1981-2010_2deg.nc
│   ├── crh_ERA5_clim_1981-2010_2deg.nc
│   ├── PI_ERA5_clim_1981-2010_2deg.nc
│   ├── sd_ERA5_clim_1981-2010_2deg.nc
│   └── ws_ERA5_clim_1981-2010_2deg.nc
│
└── ERA5_minmax/                  # ERA5 variable bounds for physical constraints
    ├── crh_min_ERA5_1981-2010_2deg.nc
    ├── PI_min_ERA5_1981-2010_2deg.nc
    ├── sd_max_ERA5_1981-2010_2deg.nc
    └── ws_min_ERA5_1981-2010_2deg.nc
```

---

## Python Packages Requirements

`xarray`, `numpy`, `xesmf`, `metpy`


---

## Input Variables
 
All preprocessing scripts read **monthly** model output. The table below lists every required input file along with which scripts use it. File name templates and variable names are configurable in the **User Settings** section of each script. All scripts auto-detect units and convert as needed (e.g., K → °C, Pa → hPa).
 
| File Template | Variable | Description | Units |
|---|---|---|---|
| `Ua_YYYY.nc` | `U` | Zonal wind (200 & 850 hPa) | m s⁻¹ |
| `Va_YYYY.nc` | `V` | Meridional wind (200 & 850 hPa) | m s⁻¹ |
| `ta_YYYY.nc` | `ta` | Multi-level air temperature | K or °C |
| `hus_YYYY.nc`| `hus` | Multi-level specific humidity | kg kg⁻¹ or g kg⁻¹ |
| `psl_YYYY.nc`| `psl` | Sea level pressure | Pa or hPa |
| `ts_YYYY.nc` | `ts` | Sea surface temperature | K or °C |
| `prw_YYYY.nc`| `prw` | Total column water vapor | kg m⁻² |
| `Model_land _mask` | `lmask` | Native-resolution land mask (ocean = 0) |
 
> **File name templates:** Replace `YYYY` with the actual year. For example, `Ua_YYYY.nc` → `Ua_1981.nc`. It is recommended to organize input files such that each file contains a single variable for a single year. However, file naming templates are fully configurable in each script’s User Settings.


---

## Workflow

The full workflow consists of three stages. Run them in order.

### Stage 1 — Conpute Input Variables

Run the three scripts for all years of interest. Each script reads model output, performs unit conversion, and saves the result as a NetCDF file. Edit the **User Settings** block at the top of each script before running.

```bash
cd calculation_vars/
python calculate_avort_ws.py   # → avort_{year}.nc, ws_{year}.nc
python calculate_swv.py        # → swv_{year}.nc
python tcpyPI/run_PI.py        # → PI_{year}.nc
```

All output files are saved to the directory specified by `output_dir` in each script (this becomes `WORK_DIR` in later stages).

---

### Stage 2 — Compute Model Climatology (`cal_TCGI_v2_part1.ipynb`)

This notebook computes the **1981–2010 monthly climatology** of all five TCGI input variables, regridded to 2° × 2°. 

**Step 1** — Regrid the input variables computed in stage 1 and `prw` to 2° × 2° and apply the land mask.
- Land points in `PI`, `swv`, and `prw` are pre-filled with their zonal-mean climatology to avoid coastal artefacts during bilinear regridding.
- Derives `crh = prw / swv` and `sd = prw − swv`.
- Output: `{WORK_DIR}/2deg/{var}_{year}_2deg.nc`

**Step 2** — Compute the monthly climatology in the model.
- Groups all years by month and averages over 1981–2010.
- Output: `{WORK_DIR}/clim/{var}_1981-2010_clim_2deg.nc`

**Step 3** *(If have multiple model ensembles)* — Average across multiple ensemble members.
- Loads climatologies for every ensemble member and computes their ensemble mean.
- Output: `{PARENT_WORK_DIR}/{var}_1981-2010_clim_2deg_mem.nc`

---

### Stage 3 — Calculate TCGI (`cal_TCGI_v2_part2.ipynb`)

This notebook applies the bias correction and computes TCGI.

**Step 1** — Regrid the input variables computed in stage 1 and `prw` to 2° × 2° (same as Part 1; skips years already processed).

**Step 2** — Compute the mean-state bias correction term for the model:
$$\Delta X = X_{\mathrm{ERA5\ clim}} - X_{\mathrm{model\ clim}}$$

**Step 3** — Apply the bias correction month-by-month, then enforce physical bounds:

| Variable | Post-correction constraint |
|---|---|
| `avort` | Take absolute value; convert to 10⁻⁵ s⁻¹; clip at 3.7 × 10⁻⁵ s⁻¹ |
| `crh` | Replace negative values with ERA5 monthly minimum; convert to % |
| `sd` | Replace positive values with ERA5 monthly maximum |
| `PI` | Replace negative values with ERA5 monthly minimum |
| `ws` | Replace negative values with ERA5 monthly minimum |

The vorticity clipping at 3.7 × 10⁻⁵ s⁻¹ improves sensitivity to near-equatorial meridional gradients while preventing over-prediction at higher latitudes (Tippett et al. 2011).


## Output Files

| File | Description |
|---|---|
| `{WORK_DIR}/2deg/{var}_{year}_2deg.nc` | Regridded (2° × 2°) input variables needed for TCGI calculation, land-masked |
| `{WORK_DIR}/clim/{var}_1981-2010_clim_2deg.nc` | 12-month climatology for single ensemble |
| `{PARENT_WORK_DIR}/{var}_1981-2010_clim_2deg.nc` | 12-month multi-ensemble mean climatology (if applicable) |
| `{WORK_DIR}/unbiased/{var}_{year}_2deg_unbiased.nc` | Bias-corrected, physically constrained |
| `{WORK_DIR}/TCGI-CRH_{year}.nc` | **TCGI-CRH** (moisture = column relative humidity) |
| `{WORK_DIR}/TCGI-SD_{year}.nc` | **TCGI-SD** (moisture = saturation deficit) |

---

## References

Fosu, B. O., Lee, C.-Y., Camargo, S. J., & Sobel, A. H. (2024). Assessing future tropical cyclone risk using downscaled CMIP6 projections. *Journal of Catastrophe Risk and Resilience*, 2024(2), 1. https://doi.org/10.55424/jcrr24009

Lee, C.-Y., Camargo, S. J., Sobel, A. H., & Tippett, M. K. (2020). Statistical–dynamical downscaling projections of tropical cyclone activity in a warming climate: Two diverging genesis scenarios. *Journal of Climate*, 33(11), 4815–4834. https://doi.org/10.1175/JCLI-D-19-0452.1

Tippett, M. K., Camargo, S. J., & Sobel, A. H. (2011). A Poisson regression index for tropical cyclone genesis and the role of large-scale vorticity in genesis. *Journal of Climate*, 24(9), 2335–2357. https://doi.org/10.1175/2010JCLI3811.1

Gilford, D. M.: pyPI (v1.3): Tropical Cyclone Potential Intensity Calculations in Python, Geosci. Model Dev., 14, 2351–2369, https://doi.org/10.5194/gmd-14-2351-2021, 2021.

Gilford, D. M. 2020: pyPI: Potential Intensity Calculations in Python, pyPI v1.3. Zenodo. http://doi.org/10.5281/zenodo.3985975



---

## License

If you use this code in your research, please cite this repository and the associated paper:

Xia, Yi, et al. "Dependence of Global Tropical Cyclones on the Tropical Pacific Mean State in the HighResMIP Models." Journal of Climate 39.4 (2026): 1099-1118.