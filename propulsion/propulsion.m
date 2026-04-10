%% Akatsuki Propulsion System Sizing & Visualization (Final Integrated Version)
clear; clc; close all;

% ----------------------------------------------------------------------
% 0. Margins (Mission-specific requirements)
% ----------------------------------------------------------------------
margins.dry_mass_margin = 0.20;      % 20% margin on dry mass
margins.propellant_margin = 0.05;    % 5% margin on propellants (OME & RCS)
margins.ullage_margin = 0.10;        % 10% volume margin for tanks
margins.pressurant_margin = 0.20;    % 20% margin on Helium mass

% ----------------------------------------------------------------------
% 1. Constants & Material Properties
% ----------------------------------------------------------------------
rho_fuel = 1010.0;           % kg/m^3 (N2H4)
rho_ox   = 1433.0;           % kg/m^3 (MON-3)
R_gas_He = 2077.3;           % J/(kg·K)
gamma_He = 1.67;
P_storage = 227e5;           % Pa (227 bar)
T_storage = 273.15 + 59.5;   % K
T_op      = 273.15 + 59.5;   % K

% RCS Pressures
op_p_i = 27.6;               % bar
op_p_f = 9;                  % bar
p_loss = 0.4;                % bar
P_switch = (op_p_i + p_loss) * 1e5; % Pa (28 bar)
P_f_rcs  = (op_p_f + p_loss) * 1e5; % Pa
B = P_switch / P_f_rcs;

% Material: Aluminum 6061-T6
sigma_yield = 240e6;         % Pa
SF          = 2.0;           % Safety Factor
limit_mm    = 1400.0;        % Total length limit in mm (1.4m)

% Baseline Mission Inputs (From Akatsuki historical data)
baseline_prop_mass = 196.3;  % kg
baseline_wet_mass = 517.6;   % kg
baseline_dry_mass = baseline_wet_mass - baseline_prop_mass;
I_sp_ome = 319.0;            % seconds
delta_v_ome = 1010.0;        % m/s
g0 = 9.81;

% ----------------------------------------------------------------------
% 2. Sizing Logic (Mass & Volume)
% ----------------------------------------------------------------------
MR_ome = exp(delta_v_ome / (I_sp_ome * g0));

% Reverse engineer the baseline RCS budget
M_final_after_VOI_baseline = baseline_wet_mass / MR_ome; 
baseline_prop_RCS = M_final_after_VOI_baseline - baseline_dry_mass; 

% Apply margins to Dry Mass and RCS Propellant
m_total_dry = baseline_dry_mass * (1.0 + margins.dry_mass_margin); 
m_prop_rcs_real = baseline_prop_RCS * (1.0 + margins.propellant_margin);

% Calculate the Real Wet Mass required to complete the OME burn 
m_post_ome_real = m_total_dry + m_prop_rcs_real; 
m_wet_real = m_post_ome_real * MR_ome;

% Calculate OME propellant needed and apply margin
m_prop_ome_real = (m_wet_real - m_post_ome_real) * (1.0 + margins.propellant_margin);

% Split OME propellant into Fuel and Oxidizer (1.8 ratio logic)
OF_ratio = 0.8;  % convert to mass ratio
m_fuel_liq_ome = m_prop_ome_real / (1+OF_ratio);
m_ox_liq = m_prop_ome_real - m_fuel_liq_ome;

% Total liquid fuel (OME fuel + all RCS propellant)
m_fuel_liq_total = m_fuel_liq_ome + m_prop_rcs_real;

% Intermediate Liquid Volumes
v_fu_liq_ome = m_fuel_liq_ome / rho_fuel;
v_ox_liq = m_ox_liq / rho_ox;
v_prop_liq_rcs = m_prop_rcs_real / rho_fuel;

% ----------------------------------------------------------------------
% 3. Pressurant and Tank Physical Volumes
% ----------------------------------------------------------------------
% Gas Mass Calculations
M_gas_real_ome = (P_switch * (v_fu_liq_ome + v_ox_liq)) / (R_gas_He * T_op) * (gamma_He / (1 - P_switch/P_storage));
M_gas_real_rcs = (P_switch * v_prop_liq_rcs) / (R_gas_He * T_op);
m_He_needed = (M_gas_real_ome + M_gas_real_rcs) * (1.0 + margins.pressurant_margin);

% Tank Volumes (Applying 10% Ullage)
V_tank_ox = v_ox_liq * (1.0 + margins.ullage_margin);
V_tank_fu = (v_fu_liq_ome + v_prop_liq_rcs + (v_prop_liq_rcs/(B-1))) * (1.0 + margins.ullage_margin); 
V_tank_he = (m_He_needed * R_gas_He * T_storage) / P_storage;

% ----------------------------------------------------------------------
% 4. Geometry (Radii & Thicknesses)
% ----------------------------------------------------------------------
% CASE 1: All Spherical
R_ox_sph = (3 * V_tank_ox / (4*pi))^(1/3);
R_fuel_sph = (3 * V_tank_fu / (4*pi))^(1/3);
R_He_sph = (3 * V_tank_he / (4*pi))^(1/3);
L_T1 = 2*(R_ox_sph + R_fuel_sph + R_He_sph);

% CASE 2: Capsule (Cylinder + Hemisphere) for Fuel
R_ox_cap = R_ox_sph;
R_He_cap = R_He_sph;
R_fuel_cap = (9 * V_tank_fu / (16*pi))^(1/3);
L_fuel_cyl = (4/9) * R_fuel_cap; 
L_T2 = 2*R_ox_cap + 2*R_He_cap + 2*R_fuel_cap + L_fuel_cyl;

% Thicknesses (in mm for the table)
% Using thin-walled sphere formula: t = P*R / (2*sigma*SF)
t_ox_mm = ((P_storage * R_ox_sph) / (2 * sigma_yield * SF)) * 1000;
t_he_mm = ((P_storage * R_He_sph) / (2 * sigma_yield * SF)) * 1000;
t_fuel_sph_mm = ((P_storage * R_fuel_sph) / (2 * sigma_yield * SF)) * 1000;
t_fuel_cap_mm = ((P_storage * R_fuel_cap) / (2 * sigma_yield * SF)) * 1000;

% --- Physical Constants ---
rho_Al = 2700; % Density of Aluminum 6061-T6 [kg/m^3]

% --- 1. SPHERICAL OXIDIZER TANK ---
% Using: M = rho * (4/3) * pi * ((r+t)^3 - r^3)
r_ox = 0.24;              % Internal radius [m]
t_ox = t_ox_mm / 1000;    % Wall thickness [m]
mass_ox = rho_Al * (4/3) * pi * ((r_ox + t_ox)^3 - r_ox^3);

% --- 2. SPHERICAL HELIUM TANK ---
r_he = 0.23;              % Internal radius [m]
t_he = t_he_mm / 1000;    % Wall thickness [m]
mass_he = rho_Al * (4/3) * pi * ((r_he + t_he)^3 - r_he^3);

% --- 3. CAPSULE FUEL TANK (Cylinder + Hemispheres) ---
% The two ends form a sphere; the middle is a hollow cylinder.
r_f = 0.33;               % Internal radius [m]
L_f = L_fuel_cyl;         % Length of the cylindrical section [m]
t_f = t_fuel_cap_mm / 1000; % Wall thickness [m]

% Mass of the two hemispherical ends (Exact Sphere Volume)
mass_fuel_ends = rho_Al * (4/3) * pi * ((r_f + t_f)^3 - r_f^3);

% Mass of the cylindrical section (Exact Hollow Cylinder Volume)
% Volume = PI * (R_outer^2 - R_inner^2) * L
mass_fuel_cyl = rho_Al * pi * ((r_f + t_f)^2 - r_f^2) * L_f;

mass_fuel_total = mass_fuel_ends + mass_fuel_cyl;

% --- 4. RESULTS DISPLAY ---
fprintf('==========================================\n');
fprintf('       TANK MASS CALCULATION (EXACT)      \n');
fprintf('==========================================\n');
fprintf('Oxidizer Tank (Spherical): %8.2f kg\n', mass_ox);
fprintf('Helium Tank   (Spherical): %8.2f kg\n', mass_he);
fprintf('Fuel Tank     (Capsule):   %8.2f kg\n', mass_fuel_total);
fprintf('------------------------------------------\n');
fprintf('TOTAL STRUCTURAL MASS:     %8.2f kg\n', mass_ox + mass_he + mass_fuel_total);
fprintf('==========================================\n');

% =========================================================================
% --- 7. CONSOLE OUTPUT FORMATTING ---
% =========================================================================
fprintf('\n========================================================\n');
fprintf('        PROPULSION SYSTEM MASS AND VOLUME BUDGET        \n');
fprintf('========================================================\n\n');

fprintf('--- 1. SYSTEM LEVEL MASSES ---\n');
fprintf('Dry Mass (Sized)           : %8.2f kg  [+%d%% Margin]\n', m_total_dry, margins.dry_mass_margin*100);
fprintf('Wet Mass (Launch)          : %8.2f kg\n', m_wet_real);
fprintf('Total Propellant Mass      : %8.2f kg\n\n', m_prop_ome_real + m_prop_rcs_real);

fprintf('--- 2. PROPELLANT BREAKDOWN ---\n');
fprintf('RCS Propellant (N2H4)      : %8.2f kg  [+%d%% Margin]\n', m_prop_rcs_real, margins.propellant_margin*100);
fprintf('OME Propellant Total       : %8.2f kg  [+%d%% Margin]\n', m_prop_ome_real, margins.propellant_margin*100);
fprintf('  -> OME Fuel (N2H4)       : %8.2f kg\n', m_fuel_liq_ome);
fprintf('  -> OME Oxidizer (MON-3)  : %8.2f kg\n\n', m_ox_liq);

fprintf('--- 3. TANK SIZING (VOLUMES & RADII) ---\n');
fprintf('Oxidizer Tank Volume       : %8.4f m^3 [+%d%% Ullage]\n', V_tank_ox, margins.ullage_margin*100);
fprintf('Oxidizer Tank Radius       : %8.4f m\n', R_ox_sph);
fprintf('Fuel Tank Volume           : %8.4f m^3 [+%d%% Ullage]\n', V_tank_fu, margins.ullage_margin*100);
fprintf('Fuel Tank Radius (Sph)     : %8.4f m\n', R_fuel_sph);
fprintf('High-Pressure Gas Tank Vol : %8.4f m^3\n', V_tank_he);
fprintf('High-Pressure Gas Tank Rad : %8.4f m\n\n', R_He_sph);

fprintf('--- 4. PRESSURANT (HELIUM) BUDGET ---\n');
fprintf('Gas Mass (Regulated Phase) : %8.2f kg\n', M_gas_real_ome);
fprintf('Gas Mass (Blowdown Phase)  : %8.2f kg\n', M_gas_real_rcs);
fprintf('Total Pressurant Mass      : %8.4f kg  [+%d%% Margin]\n', m_He_needed, margins.pressurant_margin*100);
fprintf('--- 5. TANK WALL THICKNESSES ---\n');
fprintf('Oxidizer Tank Thickness    : %8.2f mm\n', t_ox_mm);
fprintf('Helium Tank Thickness      : %8.2f mm\n', t_he_mm);
fprintf('Fuel Tank Thickness        : %8.2f mm\n\n', t_fuel_cap_mm);
fprintf('========================================================\n\n');

 