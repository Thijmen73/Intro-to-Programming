import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# 1. LOAD DATA
# ============================================================

file_name = "data.txt"

df = pd.read_csv(file_name)

print("Columns in dataset:")
print(df.columns.tolist())

print("\nNumber of rows:", len(df))


# ============================================================
# 2. CLEAN DATA
# ============================================================

# Keep only the variables needed for the analysis
df = df[["pl_name", "pl_orbper", "pl_orbsmax", "st_mass"]].copy()

# Remove missing values
df = df.dropna()

# Keep only physically meaningful positive values
df = df[
    (df["pl_orbper"] > 0) &
    (df["pl_orbsmax"] > 0) &
    (df["st_mass"] > 0)
]

# Remove duplicate planets
df = df.drop_duplicates(subset="pl_name")

print("Number of planets after cleaning:", len(df))


# ============================================================
# 3. PREPARE VARIABLES
# ============================================================

# Orbital period: days -> years
P = df["pl_orbper"].to_numpy() / 365.25

# Semi-major axis: AU
a = df["pl_orbsmax"].to_numpy()

# Host-star mass: solar masses
M = df["st_mass"].to_numpy()

# Logarithms for power-law fitting
log_a = np.log10(a)
log_P = np.log10(P)
log_M = np.log10(M)


# ============================================================
# 4. SIMPLE KEPLER FIT
# ============================================================

# Fit:
# log(P) = slope * log(a) + intercept

slope, intercept = np.polyfit(log_a, log_P, 1)

log_P_fit = slope * log_a + intercept
P_fit = 10 ** log_P_fit

# R^2
SS_res = np.sum((log_P - log_P_fit) ** 2)
SS_tot = np.sum((log_P - np.mean(log_P)) ** 2)

R2 = 1 - SS_res / SS_tot

# RMSE
RMSE = np.sqrt(np.mean((log_P - log_P_fit) ** 2))


# ============================================================
# 5. MASS-AWARE MODEL
# ============================================================

# Physical relation:
#
# P^2 = a^3 / M
#
# Therefore:
#
# log(P) = 1.5 log(a) - 0.5 log(M) + constant
#
# Fit log(P) using both log(a) and log(M).

X = np.column_stack([
    np.ones(len(log_a)),
    log_a,
    log_M
])

coefficients, _, _, _ = np.linalg.lstsq(X, log_P, rcond=None)

intercept_mass = coefficients[0]
slope_a_mass = coefficients[1]
slope_M_mass = coefficients[2]

log_P_mass_fit = X @ coefficients

SS_res_mass = np.sum((log_P - log_P_mass_fit) ** 2)
R2_mass = 1 - SS_res_mass / SS_tot

RMSE_mass = np.sqrt(
    np.mean((log_P - log_P_mass_fit) ** 2)
)


# ============================================================
# 6. RESIDUALS
# ============================================================

residuals = log_P - log_P_fit
absolute_residuals = np.abs(residuals)

df["absolute_residual"] = absolute_residuals

largest_outliers = df.sort_values(
    "absolute_residual",
    ascending=False
).head(10)


# ============================================================
# 7. PRINT RESULTS
# ============================================================

print("\n======================================")
print("KEPLER'S THIRD LAW")
print("======================================")

print(f"Number of planets: {len(df)}")
print(f"Fitted slope:      {slope:.3f}")
print(f"Theoretical slope: 1.500")
print(f"Difference:        {slope - 1.5:.3f}")
print(f"R²:                {R2:.3f}")
print(f"RMSE:              {RMSE:.3f}")


print("\n======================================")
print("MASS-AWARE MODEL")
print("======================================")

print(f"Slope for log(a):  {slope_a_mass:.3f}")
print(f"Theoretical:       1.500")

print(f"Slope for log(M):  {slope_M_mass:.3f}")
print(f"Theoretical:      -0.500")

print(f"R²:                {R2_mass:.3f}")
print(f"RMSE:              {RMSE_mass:.3f}")


print("\n======================================")
print("LARGEST OUTLIERS")
print("======================================")

print(
    largest_outliers[
        [
            "pl_name",
            "pl_orbper",
            "pl_orbsmax",
            "st_mass",
            "absolute_residual"
        ]
    ].to_string(index=False)
)


# ============================================================
# 8. PLOT 1 — OBSERVED DATA AND FIT
# ============================================================

plt.figure(figsize=(8, 6))

plt.scatter(
    a,
    P,
    s=8,
    alpha=0.35,
    label="Known exoplanets"
)

sort_order = np.argsort(a)

plt.plot(
    a[sort_order],
    P_fit[sort_order],
    linewidth=2,
    label=f"Best fit: slope = {slope:.3f}"
)
a_theory = np.logspace(
    np.log10(a.min()),
    np.log10(a.max()),
    300
)

P_theory = a_theory ** 1.5

plt.plot(
    a_theory,
    P_theory,
    linewidth=2,
    label="Kepler prediction: P = a³ᐟ²"
)

plt.xscale("log")
plt.yscale("log")

plt.xlabel("Semi-major axis [AU]")
plt.ylabel("Orbital period [years]")

plt.title("Exoplanet orbital period vs. orbital distance")

plt.legend()
plt.grid(alpha=0.3)

plt.tight_layout()
plt.show()


# ============================================================
# 9. PLOT 2 — THEORETICAL KEPLER RELATION
# ============================================================

plt.figure(figsize=(8, 6))

plt.scatter(
    a,
    P,
    s=8,
    alpha=0.35,
    label="Known exoplanets"
)

# Kepler prediction for a 1 solar-mass star:
# P = a^(3/2)

a_theory = np.logspace(
    np.log10(a.min()),
    np.log10(a.max()),
    300
)

P_theory = a_theory ** 1.5

plt.plot(
    a_theory,
    P_theory,
    linewidth=2,
    label="Kepler prediction: P = a³ᐟ²"
)

plt.xscale("log")
plt.yscale("log")

plt.xlabel("Semi-major axis [AU]")
plt.ylabel("Orbital period [years]")

plt.title("Comparison with Kepler's third law")

plt.legend()
plt.grid(alpha=0.3)

plt.tight_layout()
plt.show()


# ============================================================
# 10. PLOT 3 — RESIDUALS
# ============================================================

plt.figure(figsize=(8, 6))

plt.scatter(
    a,
    residuals,
    s=8,
    alpha=0.35
)

plt.axhline(
    0,
    linewidth=2
)

plt.xscale("log")

plt.xlabel("Semi-major axis [AU]")
plt.ylabel("Residual in log10(P)")

plt.title("Residuals from the simple Kepler fit")

plt.grid(alpha=0.3)

plt.tight_layout()
plt.show()


# ============================================================
# 11. FINAL SUMMARY
# ============================================================

print("\n======================================")
print("FINAL SUMMARY")
print("======================================")

print(f"Planets analysed:       {len(df)}")
print(f"Simple fitted slope:    {slope:.3f}")
print(f"Theoretical slope:      1.500")
print(f"Simple R²:              {R2:.3f}")
print(f"Simple RMSE:            {RMSE:.3f}")

print(f"\nMass-aware slope:        {slope_a_mass:.3f}")
print(f"Mass-aware R²:           {R2_mass:.3f}")
print(f"Mass-aware RMSE:         {RMSE_mass:.3f}")

print("\nAnalysis complete!")