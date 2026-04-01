%% Tank Sizing: Decoupled Helium Radius
% - Fuel tank: Flat-ended cylinder (Radius: R_prop)
% - Oxidizer tank: Capsule (Radius: R_prop)
% - Helium tank: Sphere (Radius: R_He)
% - Total length constraint <= 1400 mm

clear; clc;

%% 1. Constants & Material Properties
rho_fuel = 1021;        
rho_ox   = 1440;        
R_gas_He = 2077;        

P_storage = 22.7e6;     % 227 bar
P_op      = 1.45e6;     % 14.5 bar
T_storage = 273.15 + 59.5; 
T_op      = 273.15 + 59.5;

% Material: Aluminum 6061-T6
sigma_yield = 240e6;    
SF          = 2.0;      
sigma_allow = sigma_yield / SF;

f_fuel = 0.95; f_ox = 0.95;
limit_mm = 1400;

%% 2. Masses and Volumes
m_total = 320; 
O_over_F = 1.34;
m_fuel_liq = m_total / (1 + O_over_F);   
m_ox_liq   = m_total - m_fuel_liq;       

V_fuel_tank = (m_fuel_liq / rho_fuel) / f_fuel;
V_ox_tank   = (m_ox_liq / rho_ox) / f_ox;
V_ullage_total = (V_fuel_tank - (m_fuel_liq/rho_fuel)) + (V_ox_tank - (m_ox_liq/rho_ox));

m_He_needed = (P_op * V_ullage_total) / (R_gas_He * T_op);
V_He_tank = (m_He_needed * R_gas_He * T_storage) / P_storage;

%% 3. Size Helium Tank (As a Sphere)
R_He = (3 * V_He_tank / (4 * pi))^(1/3);
L_tot_He = 2 * R_He * 1000; % mm

%% 4. Solve for Propellant Tank Radius (R_prop)
% The remaining length available for Fuel and Ox
L_available_prop = limit_mm - L_tot_He;

% Max radius for Ox tank capsule
Rmax_ox = (3 * V_ox_tank / (4 * pi))^(1/3);

% Length function for Propellant tanks only
prop_length_mm = @(R) ( ...
    ( V_fuel_tank / (pi * R^2) ) + ...                         % Fuel (Cylinder)
    ( (V_ox_tank - (4/3)*pi*R^3) / (pi*R^2) + 2*R ) ...        % Ox (Capsule)
    ) * 1000;

try
    % We search for R_prop that fits in the remaining length
    R_prop = fzero(@(R) prop_length_mm(R) - L_available_prop, [0.05, Rmax_ox]);
catch
    error('Could not find a valid radius. The 1400mm limit may be too short for 320kg of propellant.');
end

%% 5. Final Dimensions & Thickness (mm)
% Fuel
L_tot_fuel = (V_fuel_tank / (pi * R_prop^2)) * 1000;
t_fuel_wall = (P_op * R_prop) / sigma_allow * 1000;
t_fuel_flat = R_prop * sqrt(0.45 * P_op / sigma_allow) * 1000; 

% Oxidizer
L_cyl_ox = (V_ox_tank - (4/3)*pi*R_prop^3) / (pi*R_prop^2);
L_tot_ox = (L_cyl_ox + 2*R_prop) * 1000;
t_ox_wall = (P_op * R_prop) / sigma_allow * 1000;
t_ox_hemi = (P_op * R_prop) / (2 * sigma_allow) * 1000;

% Helium (Sphere)
t_He_sphere = (P_storage * R_He) / (2 * sigma_allow) * 1000;

%% 6. Display Results
fprintf('--- Tank Dimensions (Decoupled Helium) ---\n');
fprintf('Propellant Tank Radius: %.2f mm\n', R_prop*1000);
fprintf('Helium Tank Radius:     %.2f mm\n', R_He*1000);

fprintf('\nFuel Tank (Flat Cylinder):\n');
fprintf('  Total Length:     %.1f mm\n', L_tot_fuel);
fprintf('  Wall Thickness:   %.3f mm\n', t_fuel_wall);
fprintf('  Flat End Thickness: %.3f mm\n', t_fuel_flat);

fprintf('\nOxidizer Tank (Capsule):\n');
fprintf('  Total Length:     %.1f mm\n', L_tot_ox);
fprintf('  Wall Thickness:   %.3f mm\n', t_ox_wall);
fprintf('  Hemi End Thickness: %.3f mm\n', t_ox_hemi);

fprintf('\nPressurizer Tank (Sphere):\n');
fprintf('  Total Length:     %.1f mm\n', L_tot_He);
fprintf('  Wall Thickness:   %.3f mm\n', t_He_sphere);

fprintf('\nTotal System Length: %.1f mm / %d mm\n', ...
    L_tot_fuel + L_tot_ox + L_tot_He, limit_mm);

fprintf('Fuel Tank Mass:       %6.2f kg (Flat Cylinder)\n', m_fuel_liq);
fprintf('Oxidizer Tank Mass:   %6.2f kg (Capsule)\n', m_ox_liq);
fprintf('Pressurizer Mass:     %6.2f kg (Sphere)\n', m_He_needed);