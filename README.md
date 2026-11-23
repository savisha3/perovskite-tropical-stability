# perovskite-tropical-stability

Code, data and a browser calculator for the review
**“Are lab stability tests tropical-ready? Mapping heat, humidity and UV stress in the tropics to perovskite solar cell degradation”**
(S. Mahalingam, submitted).

The **stress-budget model** converts an accelerated laboratory test of a perovskite solar cell
(damp heat, light and heat, dark heat or thermal cycling) into the number of **field years** it is
equivalent to at a given site and mounting, using ten years of hourly NASA POWER weather data.

**Calculator (no installation):** https://savisha3.github.io/perovskite-tropical-stability/

## What is in this repository

| Path | Content |
|---|---|
| `index.html`, `calculator/model.js` | Browser calculator. Runs entirely in your browser; works for the 20 built-in sites, for any NASA POWER hourly file you upload, and (where the browser allows it) for any latitude and longitude. |
| `stressbudget.py` | The stress-budget model in Python (used for all results in the review). |
| `data/weather/` | NASA POWER hourly data, 2015–2024, for the 17 Southeast Asian cities and three reference sites (Miami, Riyadh, Freiburg). |
| `data/cities.csv` | Site list: coordinates, module tilt, Köppen–Geiger class and data file. |
| `data/cities_table1.csv` | Table 1 of the review. |
| `data/literature/benchmark_results.csv` | The 388 published accelerated-test results (Dataset S1, sheet “Benchmark”), with field-year conversions. |
| `data/literature/outdoor_tests.csv` | The 179 outdoor tests of perovskite devices (Dataset S1, sheet “Outdoor”). |
| `tests/reproduce_table1.py` | Re-runs the Python model for all 20 sites and checks the results against Table 1. |
| `tests/check_js_model.js` | Checks the browser calculator against Table 1. |

## Model in brief

For every hour of the weather record:

1. Plane-of-array irradiance *G*: Erbs decomposition of global horizontal irradiance and an isotropic sky (pvlib), module facing the equator.
2. Module temperature *T*<sub>m</sub>: Faiman model; open rack *U*<sub>0</sub> = 25 W m<sup>−2</sup> K<sup>−1</sup>, *U*<sub>1</sub> = 6.84 W s m<sup>−3</sup> K<sup>−1</sup>; rooftop *U*<sub>0</sub> = 20, *U*<sub>1</sub> = 0; 10 m wind speed scaled to 2 m by 0.72.
3. Module-surface relative humidity: ambient water-vapour pressure held constant and evaluated at *T*<sub>m</sub> (Magnus form, WMO).

The annual stress is then expressed as equivalent laboratory hours (or cycles):

| Test | Model | Default parameters |
|---|---|---|
| Damp heat | Peck: (RH/85)<sup>*n*</sup> · Arrhenius | 85 °C, 85% RH, *n* = 2.66, *E*<sub>a</sub> = 0.6 eV |
| Light and heat | (*G*/1000)<sup>0.72</sup> · Arrhenius | 85 °C, 1 sun, *E*<sub>a</sub> = 0.61 eV (0.43 eV for tandems) |
| Dark heat | Arrhenius | 85 °C, *E*<sub>a</sub> = 0.6 eV |
| Thermal cycling | Coffin–Manson: (Δ*T*<sub>day</sub>/125 K)<sup>2</sup> | −40 to 85 °C |

Field years = test duration ÷ equivalent hours per year. A test that ended before the device failed
gives a lower bound, not a predicted lifetime. Full equations, checks and sensitivity analyses are in
Section 3 of the review and ESI Section S2.

## Running the Python model

```bash
pip install -r requirements.txt
python stressbudget.py data/weather/kuala_lumpur.csv --tilt 10
python tests/reproduce_table1.py      # checks all 20 sites against Table 1
```

Any NASA POWER hourly CSV with the parameters `T2M`, `RH2M`, `ALLSKY_SFC_SW_DWN` and `WS10M`
(local solar time or UTC) can be used. Files can be downloaded from the
[NASA POWER Data Access Viewer](https://power.larc.nasa.gov/data-access-viewer/).

To check the browser calculator against the same values: `node tests/check_js_model.js`.

## How to cite

Please cite the review and this repository (see `CITATION.cff`; an archived version with a DOI is
on Zenodo). Please also cite the NASA POWER Project when using the weather data.

## Licence

Code: MIT (see `LICENSE`). Literature datasets and tables: CC BY 4.0 (see `data/LICENSE.md`).
Weather data: NASA POWER, no restrictions on use; please cite the NASA POWER Project.

## Contact

Savisha Mahalingam, Institute of Sustainable Energy, Universiti Tenaga Nasional (UNITEN), Malaysia —
savisha@uniten.edu.my
