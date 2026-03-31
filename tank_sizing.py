import math
import scipy.optimize as opt

# ----------------------------------------------------------------------
# 1. Constants & Material Properties
# ----------------------------------------------------------------------
rho_fuel = 1021.0           # kg/m^3
rho_ox   = 1440.0           # kg/m^3
R_gas_He = 2077.0           # J/(kg·K)

P_storage = 22.7e6          # Pa (227 bar)
P_op      = 1.45e6          # Pa (14.5 bar)
T_storage = 273.15 + 59.5   # K
T_op      = 273.15 + 59.5   # K

# Material: Aluminum 6061-T6
sigma_yield = 240e6          # Pa
SF          = 2.0
sigma_allow = sigma_yield / SF

f_fuel = 0.95
f_ox   = 0.95
limit_mm = 1400.0            # total length limit in mm

# ----------------------------------------------------------------------
# 2. Masses and Volumes
# ----------------------------------------------------------------------
m_total = 320.0 * 1.03  # kg (320kg propellant + 3% margin)
O_over_F = 0.8 

m_fuel_liq = m_total / (1.0 + O_over_F)
m_ox_liq   = m_total - m_fuel_liq

V_fuel_tank = (m_fuel_liq / rho_fuel) / f_fuel
V_ox_tank   = (m_ox_liq   / rho_ox)   / f_ox

V_ullage_total = (V_fuel_tank - (m_fuel_liq / rho_fuel)) + \
                 (V_ox_tank   - (m_ox_liq   / rho_ox))

m_He_needed = (P_op * V_ullage_total) / (R_gas_He * T_op)

V_He_tank = (m_He_needed * R_gas_He * T_storage) / P_storage

# ----------------------------------------------------------------------
# 3. Size Helium Tank (Sphere)
# ----------------------------------------------------------------------
R_He = (3.0 * V_He_tank / (4.0 * math.pi)) ** (1.0/3.0)
L_tot_He = 2.0 * R_He * 1000.0            # mm

# ----------------------------------------------------------------------
# 4. Solve for Propellant Tank Radius (R_prop)
# ----------------------------------------------------------------------
L_available_prop = limit_mm - L_tot_He

# Upper bound: radius if oxidizer tank were a sphere
Rmax_ox = (3.0 * V_ox_tank / (4.0 * math.pi)) ** (1.0/3.0)

# Length function for propellant tanks (returns mm)
def prop_length_mm(R):
    # Fuel tank: flat-ended cylinder
    L_fuel = V_fuel_tank / (math.pi * R**2)          # meters
    # Oxidizer tank: capsule (cylinder + two hemispheres)
    V_cyl_ox = V_ox_tank - (4.0/3.0) * math.pi * R**3
    L_cyl_ox = V_cyl_ox / (math.pi * R**2)           # meters
    L_ox = L_cyl_ox + 2.0 * R                         # meters
    return (L_fuel + L_ox) * 1000.0                  # mm

# Find root such that prop_length_mm(R) = L_available_prop
def objective(R):
    return prop_length_mm(R) - L_available_prop

try:
    # Bracket must have opposite signs; check at bounds
    f_low = objective(0.05)
    f_high = objective(Rmax_ox)
    if f_low * f_high > 0:
        raise ValueError("Function does not change sign in the bracket.")
    R_prop = opt.brentq(objective, 0.05, Rmax_ox)
except (ValueError, RuntimeError) as e:
    raise RuntimeError("Could not find a valid radius. "
                       "The 1400mm limit may be too short for 320kg of propellant.") from e

# ----------------------------------------------------------------------
# 5. Final Dimensions & Thicknesses (mm)
# ----------------------------------------------------------------------
# Fuel tank (flat-ended cylinder)
L_tot_fuel = (V_fuel_tank / (math.pi * R_prop**2)) * 1000.0
t_fuel_wall = (P_op * R_prop) / sigma_allow * 1000.0
t_fuel_flat = R_prop * math.sqrt(0.45 * P_op / sigma_allow) * 1000.0

# Oxidizer tank (capsule)
L_cyl_ox = (V_ox_tank - (4.0/3.0) * math.pi * R_prop**3) / (math.pi * R_prop**2)
L_tot_ox = (L_cyl_ox + 2.0 * R_prop) * 1000.0
t_ox_wall = (P_op * R_prop) / sigma_allow * 1000.0
t_ox_hemi = (P_op * R_prop) / (2.0 * sigma_allow) * 1000.0

# Helium tank (sphere)
t_He_sphere = (P_storage * R_He) / (2.0 * sigma_allow) * 1000.0

# ----------------------------------------------------------------------
# 6. Display Results
# ----------------------------------------------------------------------
print("--- Tank Dimensions (Decoupled Helium) ---")
print(f"Propellant Tank Radius: {R_prop*1000:.2f} mm")
print(f"Helium Tank Radius:     {R_He*1000:.2f} mm")

print("\nFuel Tank (Flat Cylinder):")
print(f"  Total Length:     {L_tot_fuel:.1f} mm")
print(f"  Wall Thickness:   {t_fuel_wall:.3f} mm")
print(f"  Flat End Thickness: {t_fuel_flat:.3f} mm")

print("\nOxidizer Tank (Capsule):")
print(f"  Total Length:     {L_tot_ox:.1f} mm")
print(f"  Wall Thickness:   {t_ox_wall:.3f} mm")
print(f"  Hemi End Thickness: {t_ox_hemi:.3f} mm")

print("\nPressurizer Tank (Sphere):")
print(f"  Total Length:     {L_tot_He:.1f} mm")
print(f"  Wall Thickness:   {t_He_sphere:.3f} mm")

print(f"\nTotal System Length: {L_tot_fuel + L_tot_ox + L_tot_He:.1f} mm / {limit_mm} mm")

print(f"Fuel Tank Mass:       {m_fuel_liq:6.2f} kg (Flat Cylinder)")
print(f"Oxidizer Tank Mass:   {m_ox_liq:6.2f} kg (Capsule)")
print(f"Pressurizer Mass:     {m_He_needed:6.2f} kg (Sphere)")