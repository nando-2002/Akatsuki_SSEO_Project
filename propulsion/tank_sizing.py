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
rho_al = 2700.0              # kg/m^3

f_fuel = 0.95
f_ox   = 0.95
limit_mm = 1400.0            # total length limit in mm

I_sp_ome = 319.0                 # seconds (assumed specific impulse)
I_sp_rcs = 290.0                 # seconds (assumed specific impulse for RCS)

delta_v_ome = 1010.0                # m/s (delta-V requirement)
delta_v_rcs = 200.0                # m/s (delta-V requirement for RCS)
g0 = 9.81                      # m/s^2 (standard gravity)

MR_ome = math.exp(delta_v_ome / (I_sp_ome * g0))  # Mass ratio
MR_rcs = math.exp(delta_v_rcs / (I_sp_rcs * g0))  # Mass ratio for RCS
m_total_dry = 320.0 * 1.20  # kg (320kg propellant + 20% margin)
m_rcs = MR_rcs * m_total_dry - m_total_dry
m_ome = MR_ome * m_total_dry * MR_rcs - m_total_dry*MR_rcs

m_prop = (m_ome + m_rcs)*(1.0 + 0.02 + 0.05)  # kg (propellant mass with 2% margin)
m_total_wet = m_total_dry + m_prop + m_rcs  # kg (total wet mass including RCS propellant)
# ----------------------------------------------------------------------
# 2. Propellant & Pressurant Masses / Volumes
# ----------------------------------------------------------------------
#m_total_dry = 320.0 * 1.03  # kg (320kg propellant + 3% margin)
O_over_F = 0.8

m_fuel_liq = m_prop / (1.0 + O_over_F)
m_ox_liq   = m_prop - m_fuel_liq

V_fuel_tank = (m_fuel_liq / rho_fuel) / f_fuel* 1.1  # m^3 (10% margin on fuel tank volume)
V_ox_tank   = (m_ox_liq   / rho_ox)   / f_ox*    1.1  # m^3 (10% margin on oxidizer tank volume)

# Ullage volume (gas needed at operating pressure)
V_ullage_total = (V_fuel_tank - (m_fuel_liq / rho_fuel)) + \
                 (V_ox_tank   - (m_ox_liq   / rho_ox))

m_He_needed = (P_op * V_ullage_total) / (R_gas_He * T_op)*1.2  # kg (20% margin on pressurant mass)
V_He_tank = (m_He_needed * R_gas_He * T_storage) / P_storage

V_comb = V_fuel_tank + V_He_tank          # total volume of combined tank

# ----------------------------------------------------------------------
# 3. Determine Feasible Combined Tank Radius
# ----------------------------------------------------------------------
# Minimum oxidizer tank length = diameter of a sphere of volume V_ox_tank
R_ox_min_sphere = (3.0 * V_ox_tank / (4.0 * math.pi)) ** (1.0/3.0)
L_ox_min = 2.0 * R_ox_min_sphere

# Maximum allowed length for combined tank (to leave enough space for oxidizer)
L_comb_max = limit_mm / 1000.0 - L_ox_min
if L_comb_max <= 0:
    raise ValueError("Total length limit is too short to fit both tanks.")

# Minimum radius of combined tank to satisfy length limit
R_com_min = math.sqrt(V_comb / (math.pi * L_comb_max))

# Choose a radius (e.g., 10% larger than the minimum, adjust as needed)
R_com = R_com_min * 1.10

# Actual lengths
L_comb = V_comb / (math.pi * R_com**2)               # combined tank length
L_remain = limit_mm / 1000.0 - L_comb               # leftover length for oxidizer

# ----------------------------------------------------------------------
# 4. Sizing the Oxidizer Tank (Capsule) for the Given Length
# ----------------------------------------------------------------------
def oxidizer_radius(L_ox, V_ox):
    """
    For a capsule tank with total length L_ox (m) and volume V_ox (m³),
    return the required radius (m). Assumes L_ox >= 2*(3V_ox/(4π))^(1/3).
    """
    R_max = L_ox / 2.0
    # Check feasibility
    R_min_sphere = (3.0 * V_ox / (4.0 * math.pi)) ** (1.0/3.0)
    if L_ox < 2.0 * R_min_sphere - 1e-9:
        return None
    if abs(L_ox - 2.0 * R_min_sphere) < 1e-9:
        return R_min_sphere

    # Cubic: 2R³ - 3 L_ox R² + 3 V_ox/π = 0
    def f(R):
        return 2.0 * R**3 - 3.0 * L_ox * R**2 + 3.0 * V_ox / math.pi

    try:
        return opt.brentq(f, 1e-9, R_max)
    except ValueError:
        # If f(R_max) is zero, the root is at the upper bound
        if abs(f(R_max)) < 1e-9:
            return R_max
        else:
            raise

R_ox = oxidizer_radius(L_remain, V_ox_tank)
if R_ox is None:
    raise RuntimeError("Could not fit oxidizer tank in the available length.")

L_ox = L_remain   # by construction, the oxidizer length equals the leftover length

# ----------------------------------------------------------------------
# 5. Thicknesses (All in meters)
# ----------------------------------------------------------------------
# Combined tank (cylindrical, two flat ends, one bulkhead)
t_wall_com = (P_storage * R_com) / sigma_allow
t_end_fuel = R_com * math.sqrt(0.45 * P_op / sigma_allow)
t_end_He   = R_com * math.sqrt(0.45 * P_storage / sigma_allow)
t_bulkhead = R_com * math.sqrt(0.45 * (P_storage - P_op) / sigma_allow)

# Oxidizer tank (capsule: cylinder + two hemispheres)
t_wall_ox = (P_op * R_ox) / sigma_allow
t_hemi_ox = (P_op * R_ox) / (2.0 * sigma_allow)

# ----------------------------------------------------------------------
# 6. Tank Shell Masses (Aluminum)
# ----------------------------------------------------------------------
# Combined tank volumes
L_fuel = V_fuel_tank / (math.pi * R_com**2)
L_He   = V_He_tank   / (math.pi * R_com**2)

vol_cyl_com = 2.0 * math.pi * R_com * (L_fuel + L_He) * t_wall_com
vol_end_fuel = math.pi * R_com**2 * t_end_fuel
vol_end_He   = math.pi * R_com**2 * t_end_He
vol_bulkhead = math.pi * R_com**2 * t_bulkhead
mass_com_tank = rho_al * (vol_cyl_com + vol_end_fuel + vol_end_He + vol_bulkhead)

# Oxidizer tank volumes
L_cyl_ox = (V_ox_tank - (4.0/3.0) * math.pi * R_ox**3) / (math.pi * R_ox**2)
vol_cyl_ox = 2.0 * math.pi * R_ox * L_cyl_ox * t_wall_ox
vol_hemi_ox = 4.0 * math.pi * R_ox**2 * t_hemi_ox
mass_ox_tank = rho_al * (vol_cyl_ox + vol_hemi_ox)

# ----------------------------------------------------------------------
# 7. Propulsion System Mass Budget
# ----------------------------------------------------------------------
total_propellant = m_fuel_liq + m_ox_liq
total_pressurant = m_He_needed
total_tank_structure = mass_com_tank + mass_ox_tank
total_propulsion_system = total_propellant + total_pressurant + total_tank_structure

target_dry_mass = 320.0  *1.2  # kg + 20% margin
remaining_dry_mass = target_dry_mass - total_propulsion_system

# ----------------------------------------------------------------------
# 8. Print Results
# ----------------------------------------------------------------------
print("\n--- Combined Tank (Fuel + Helium) ---")
print(f"Radius:          {R_com*1000:.2f} mm")
print(f"Fuel length:     {L_fuel*1000:.1f} mm")
print(f"Helium length:   {L_He*1000:.1f} mm")
print(f"Total length:    {L_comb*1000:.1f} mm")
print(f"\nOxidizer Tank (Capsule):")
print(f"Radius:          {R_ox*1000:.2f} mm")
print(f"Total length:    {L_ox*1000:.1f} mm")
print(f"Total system length: {(L_comb+L_ox)*1000:.1f} mm / {limit_mm} mm")

print("\n--- Thicknesses ---")
print(f"Combined tank wall:        {t_wall_com*1000:.3f} mm")
print(f"Combined tank fuel end:    {t_end_fuel*1000:.3f} mm")
print(f"Combined tank He end:      {t_end_He*1000:.3f} mm")
print(f"Bulkhead thickness:        {t_bulkhead*1000:.3f} mm")
print(f"Oxidizer tank wall:        {t_wall_ox*1000:.3f} mm")
print(f"Oxidizer hemispherical ends: {t_hemi_ox*1000:.3f} mm")

print("\n--- Masses (with margin) ---")
print(f"Fuel propellant:        {m_fuel_liq:6.2f} kg")
print(f"Oxidizer propellant:    {m_ox_liq:6.2f} kg")
print(f"Pressurant (He):        {m_He_needed:6.2f} kg")
print(f"Combined tank shell:    {mass_com_tank:6.2f} kg")
print(f"Oxidizer tank shell:    {mass_ox_tank:6.2f} kg")
print(f"Total tank structure:   {total_tank_structure:6.2f} kg")
print(f"Total propulsion system: {total_propulsion_system:6.2f} kg")

print(f"\nTarget dry mass (given): {target_dry_mass:6.2f} kg")
if remaining_dry_mass >= 0:
    print(f"Remaining for spacecraft (excluding propulsion): {remaining_dry_mass:6.2f} kg")
else:
    print(f"WARNING: Propulsion system exceeds target dry mass by {-remaining_dry_mass:6.2f} kg")