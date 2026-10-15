# Testing Kepler's Third Law with Known Exoplanets

**MAT2007 – Introduction to Programming**  
**Thijmen Graper (6320296)**

## Project overview

This project investigates how closely the orbital parameters of known exoplanets agree with Kepler's third law. Using data from the NASA Exoplanet Archive, I compare orbital period with semi-major axis, first on its own and then while also accounting for the mass of the host star.

For orbital period in years, semi-major axis in astronomical units (AU), and stellar mass in solar masses, the approximate relation is:

`P² = a³ / M★`

The analysis uses logarithmic linear regressions to compare the fitted exponents with the theoretical values of **+1.5** for semi-major axis and **−0.5** for stellar mass.

## Files

- `code.txt` — Python script containing the data cleaning, regressions, calculations, and plots. Can be found under the folder "programming project".
- `data.txt` — CSV dataset downloaded from the NASA Exoplanet Archive. Can be found under the folder "programming project".
- `README.md` — Instructions and information about the project.
- The report PDF — Scientific background, methods, results, and discussion (included separately in the repository).

## Requirements

- Python 3
- NumPy
- pandas
- Matplotlib

The project was developed to run locally in Visual Studio Code, but it can also be run from a terminal. No other Python packages are required.

Install the dependencies with:

```bash
python -m pip install numpy pandas matplotlib
```

On Windows, `py` can be used instead of `python` if needed.

## How to run

1. Download the repository and keep `code.txt` and `data.txt` in the **same folder**.
2. Open that folder in Visual Studio Code (or a terminal).
3. Run:

   ```bash
   python code.txt
   ```

Although the script is supplied as `code.txt`, it contains Python code and runs with the command above. It can also be renamed to `code.py` if preferred. The script reads `data.txt` from the current working folder.

The terminal prints the fitted coefficients, R², RMSE, and the ten planets with the largest absolute residuals from the simple fit. Three plots are displayed in Matplotlib windows. **The script displays the plots but does not automatically save image files.** Depending on the plotting backend, close a plot window to continue to the next one.

## Dataset and cleaning

**Source:** NASA Exoplanet Archive, *Planetary Systems Composite Parameters* (`pscomppars`) table.

**Exact CSV download:**  
https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select%20pl_name%2Cpl_orbper%2Cpl_orbsmax%2Cst_mass%20from%20pscomppars&format=csv

The downloaded dataset contains these columns:

- `pl_name`: planet name
- `pl_orbper`: orbital period in days
- `pl_orbsmax`: semi-major axis in AU
- `st_mass`: host-star mass in solar masses

The included dataset has **6,375 entries**. The script removes rows with missing values, excludes non-positive values, and removes duplicate planet names. This leaves **5,601 planets** for analysis. All cleaning is performed in the code, so no manual cleaning of `data.txt` is needed.

The NASA archive is updated over time, so downloading the query again may produce a different number of entries. Use the included `data.txt` to reproduce the results reported here.

## Analysis and results

Two models are fitted using NumPy:

1. **Simple fit:** `log10(P) = b0 + ba × log10(a)`
2. **Mass-aware fit:** `log10(P) = c0 + ca × log10(a) + cM × log10(M★)`

The code calculates R² and RMSE using **log10 of the orbital period**, not period in years. It also plots the period–distance relationship, the theoretical curve for a one-solar-mass star, and residuals from the simple model.

| Result | Simple fit | Mass-aware fit |
|---|---:|---:|
| Semi-major-axis exponent | 1.466 | 1.497 |
| Stellar-mass exponent | — | −0.474 |
| R² | 0.990 | 0.999 |
| RMSE (log10 period) | 0.093 | 0.035 |

Accounting for host-star mass reduces the RMSE by approximately **62%**. Some catalogue semi-major axes may themselves have been calculated using Kepler's law, so the analysis measures agreement among reported parameters rather than providing a fully independent test.
