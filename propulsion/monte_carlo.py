import math
import scipy.optimize as opt
import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# 0. Margins (Nominal values – will be varied in Monte Carlo)
# ----------------------------------------------------------------------
nominal_margins = {
    'dry_mass_margin': 0.10,
    'propellant_margin': 0.01,
    'ullage_margin': 0.10,
    'fuel_volume_margin': 0.05,
    'oxidizer_volume_margin': 0.05,
    'pressurant_margin': 0.10
}

# Range for each margin: (min, max) – can be adjusted
margin_ranges = {
    'dry_mass_margin': (0.0, 0.2),
    'propellant_margin': (0.0, 0.02),
    'ullage_margin': (0.0, 0.2),
    'fuel_volume_margin': (0.0, 0.1),
    'oxidizer_volume_margin': (0.0, 0.1),
    'pressurant_margin': (0.0, 0.2)
}

# ----------------------------------------------------------------------
# 1. Constants & Material Properties (fixed)
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
I_sp_rcs = 302                 # seconds (assumed specific impulse for RCS)

delta_v_ome = 1010.0                # m/s (delta-V requirement)
delta_v_rcs = 200.0                # m/s (delta-V requirement for RCS)
g0 = 9.81                      # m/s^2 (standard gravity)

# ----------------------------------------------------------------------
# Function to run sizing for given margins
# ----------------------------------------------------------------------
def run_sizing(margins):
    """
    Perform tank sizing with given margins dictionary.
    Returns a dictionary of results or None if infeasible.
    """
    try:
        # Extract margins
        dry_mass_margin = margins['dry_mass_margin']
        propellant_margin = margins['propellant_margin']
        ullage_margin = margins['ullage_margin']
        fuel_volume_margin = margins['fuel_volume_margin']
        oxidizer_volume_margin = margins['oxidizer_volume_margin']
        pressurant_margin = margins['pressurant_margin']

        # Mass ratios
        MR_ome = math.exp(delta_v_ome / (I_sp_ome * g0))
        MR_rcs = math.exp(delta_v_rcs / (I_sp_rcs * g0))

        m_total_dry = 320.0 * (1.0 + dry_mass_margin)  # kg (dry mass with margin)
        m_rcs = MR_rcs * m_total_dry - m_total_dry
        m_ome = MR_ome * m_total_dry * MR_rcs - m_total_dry * MR_rcs

        m_prop = (m_ome + m_rcs) * (1.0 + propellant_margin + ullage_margin)  # total propellant (OME+RCS) with margins
        m_total_wet = m_total_dry + m_prop + m_rcs  # not used later, but kept for completeness

        #O_over_F = 0.8
        O_over_F = rho_ox / rho_fuel  # convert to mass ratio
        m_fuel_liq = m_prop / (1.0 + O_over_F)
        m_ox_liq   = m_prop - m_fuel_liq

        V_fuel_tank = (m_fuel_liq / rho_fuel) / f_fuel * (1.0 + fuel_volume_margin)
        V_ox_tank   = (m_ox_liq   / rho_ox)   / f_ox   * (1.0 + oxidizer_volume_margin)

        # Ullage volume (gas needed at operating pressure)
        V_ullage_total = (V_fuel_tank - (m_fuel_liq / rho_fuel)) + \
                         (V_ox_tank   - (m_ox_liq   / rho_ox))

        m_He_needed = (P_op * V_ullage_total) / (R_gas_He * T_op) * (1.0 + pressurant_margin)
        V_He_tank = (m_He_needed * R_gas_He * T_storage) / P_storage

        V_comb = V_fuel_tank + V_He_tank

        # Geometry sizing
        R_ox_min_sphere = (3.0 * V_ox_tank / (4.0 * math.pi)) ** (1.0/3.0)
        L_ox_min = 2.0 * R_ox_min_sphere
        L_comb_max = limit_mm / 1000.0 - L_ox_min
        if L_comb_max <= 0:
            return None

        R_com_min = math.sqrt(V_comb / (math.pi * L_comb_max))
        R_com = R_com_min * 1.10
        L_comb = V_comb / (math.pi * R_com**2)
        L_remain = limit_mm / 1000.0 - L_comb

        # Oxidizer radius
        def oxidizer_radius(L_ox, V_ox):
            R_max = L_ox / 2.0
            R_min_sphere = (3.0 * V_ox / (4.0 * math.pi)) ** (1.0/3.0)
            if L_ox < 2.0 * R_min_sphere - 1e-9:
                return None
            if abs(L_ox - 2.0 * R_min_sphere) < 1e-9:
                return R_min_sphere
            def f(R):
                return 2.0 * R**3 - 3.0 * L_ox * R**2 + 3.0 * V_ox / math.pi
            try:
                return opt.brentq(f, 1e-9, R_max)
            except ValueError:
                if abs(f(R_max)) < 1e-9:
                    return R_max
                else:
                    raise
        R_ox = oxidizer_radius(L_remain, V_ox_tank)
        if R_ox is None:
            return None

        L_ox = L_remain

        # Thicknesses
        t_wall_com = (P_storage * R_com) / sigma_allow
        t_end_fuel = R_com * math.sqrt(0.45 * P_op / sigma_allow)
        t_end_He   = R_com * math.sqrt(0.45 * P_storage / sigma_allow)
        t_bulkhead = R_com * math.sqrt(0.45 * (P_storage - P_op) / sigma_allow)

        t_wall_ox = (P_op * R_ox) / sigma_allow
        t_hemi_ox = (P_op * R_ox) / (2.0 * sigma_allow)

        # Tank masses
        L_fuel = V_fuel_tank / (math.pi * R_com**2)
        L_He   = V_He_tank   / (math.pi * R_com**2)

        vol_cyl_com = 2.0 * math.pi * R_com * (L_fuel + L_He) * t_wall_com
        vol_end_fuel = math.pi * R_com**2 * t_end_fuel
        vol_end_He   = math.pi * R_com**2 * t_end_He
        vol_bulkhead = math.pi * R_com**2 * t_bulkhead
        mass_com_tank = rho_al * (vol_cyl_com + vol_end_fuel + vol_end_He + vol_bulkhead)

        L_cyl_ox = (V_ox_tank - (4.0/3.0) * math.pi * R_ox**3) / (math.pi * R_ox**2)
        vol_cyl_ox = 2.0 * math.pi * R_ox * L_cyl_ox * t_wall_ox
        vol_hemi_ox = 4.0 * math.pi * R_ox**2 * t_hemi_ox
        mass_ox_tank = rho_al * (vol_cyl_ox + vol_hemi_ox)

        total_propellant = m_fuel_liq + m_ox_liq
        total_pressurant = m_He_needed
        total_tank_structure = mass_com_tank + mass_ox_tank
        total_propulsion_system = total_propellant + total_pressurant + total_tank_structure

        target_dry_mass = 320.0 * (1.0 + dry_mass_margin)
        remaining_dry_mass = target_dry_mass - total_propulsion_system

        return {
            'total_propulsion': total_propulsion_system,
            'total_tank_structure': total_tank_structure,
            'total_propellant': total_propellant,
            'remaining_dry_mass': remaining_dry_mass,
            'dry_mass_margin': dry_mass_margin,
            'propellant_margin': propellant_margin,
            'ullage_margin': ullage_margin,
            'fuel_volume_margin': fuel_volume_margin,
            'oxidizer_volume_margin': oxidizer_volume_margin,
            'pressurant_margin': pressurant_margin
        }
    except Exception:
        return None

# ----------------------------------------------------------------------
# Monte Carlo Simulation
# ----------------------------------------------------------------------
n_samples = 20000
results = []

for _ in range(n_samples):
    # Sample each margin uniformly within its range
    margins = {}
    for key, (low, high) in margin_ranges.items():
        margins[key] = np.random.uniform(low, high)
    res = run_sizing(margins)
    if res is not None:
        results.append(res)

print(f"Completed {len(results)} feasible designs out of {n_samples} samples.")

# Convert to arrays for plotting
if len(results) > 0:
    total_propulsion = np.array([r['total_propulsion'] for r in results])
    total_tank_struct = np.array([r['total_tank_structure'] for r in results])
    total_propellant = np.array([r['total_propellant'] for r in results])
    remaining_dry = np.array([r['remaining_dry_mass'] for r in results])

    # Scatter plots
    plt.figure(figsize=(15, 5))

    plt.subplot(1, 3, 1)
    plt.scatter(total_tank_struct, total_propulsion, alpha=0.5, s=5)
    plt.xlabel('Total Tank Structure Mass (kg)')
    plt.ylabel('Total Propulsion System Mass (kg)')
    plt.title('Propulsion vs Tank Structure Mass')
    plt.grid(True)

    plt.subplot(1, 3, 2)
    plt.scatter(total_propellant, total_propulsion, alpha=0.5, s=5)
    plt.xlabel('Total Propellant Mass (kg)')
    plt.ylabel('Total Propulsion System Mass (kg)')
    plt.title('Propulsion vs Propellant Mass')
    plt.grid(True)

    plt.subplot(1, 3, 3)
    plt.scatter(total_propulsion, remaining_dry, alpha=0.5, s=5)
    plt.xlabel('Total Propulsion System Mass (kg)')
    plt.ylabel('Remaining Dry Mass (kg)')
    plt.title('Remaining Dry Mass vs Propulsion Mass')
    plt.grid(True)

    plt.tight_layout()
    plt.show()

    # Histograms of key margins (optional)
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    for i, key in enumerate(margin_ranges.keys()):
        values = [r[key] for r in results]
        axes[i].hist(values, bins=30, alpha=0.7, edgecolor='k')
        axes[i].set_title(key)
        axes[i].set_xlabel('Margin')
        axes[i].set_ylabel('Frequency')
    plt.tight_layout()
    plt.show()

else:
    print("No feasible designs found. Try adjusting the margin ranges or length limit.")