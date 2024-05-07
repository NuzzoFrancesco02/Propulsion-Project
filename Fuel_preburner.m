%% FUEL LINE

% Clear command window and workspace
clear
clc

% Define parameters for the fuel line
D_i = convlength(10.43,'in','m'); % Inner diameter (m)
L = convlength(4.37,'in','m'); % Length (m)
O_F = 0.869086; % Oxidizer-to-fuel ratio
%O_F = 0.1:0.01:1.5;
P = 4793; % Pressure (psia)

% Define properties for the fuel and oxidizer streams
H = struct('m_dot',[],'T',[],'P',[],'rho',[],'h',[],'Mm',[]);
O = H;
H.m_dot = 34.9; % Fuel mass flow rate (kg/s)
H.T = 148.2; % Fuel temperature (K)
H.P = 36.61; % Fuel pressure (Mpa)
H.rho = 42.42; % Fuel density (kg/m^3)
H.Mm = 2.016;
H.h = -2.503*H.Mm; % Fuel enthalpy (kJ/mol)

O.m_dot = 30.4; % Oxidizer mass flow rate (kg/s)
O.T = 115.4; % Oxidizer temperature (K)
O.P = 38.86; % Oxidizer pressure (Mpa)
O.rho = 1114.6; % Oxidizer density (kg/m^3)
O.Mm = 31.9988;
O.h = -0.3427*O.Mm; % Oxidizer enthalpy (kJ/mol)

% Calculate total mass flow rate and area
%m_tot = 60:0.1:70;
m_tot = H.m_dot + O.m_dot;
A = D_i^2/4*pi;
%A = 0.01:0.001:0.1; 
mass_area_ratio = m_tot ./ A;

% Perform calculations using a function called 'cea_preburner_input'
[sol, products] = cea_preburner_input(H, O, O_F, P, mass_area_ratio);
CEA_results_fuel = struct2table(sol) % Convert structure to table
struct2table(products); % Convert structure to table

% Calculate properties for the cooling stream
H_cool = struct('m_dot',[],'T',[],'P',[],'rho',[],'h',[]);
H_cool.m_dot = 2.7;
H_cool.T = 51.5;
H_cool.P = 41.07; % Pressure (Mpa)
H_cool.rho = 78.82;
H_cool.h = -7.539048;

% Update total mass flow rate and percentages
m_tot_new = m_tot + H_cool.m_dot;
perc_H_cool = H_cool.m_dot ./ m_tot_new;
perc_H = m_tot .* products.Mass_fraction(1) ./ m_tot_new;
perc_H20 = m_tot .* products.Mass_fraction(2) ./ m_tot_new;

% Define properties and store fuel information
vec_size = 1;
if isvector(mass_area_ratio)
    vec_size = length(mass_area_ratio);
elseif isvector(O_F)
    vec_size = length(O_F);
end
properties = struct('P_bar',sol.P,'T_before_main',[],'cp',sol.cp,'Mm',sol.Mm,'m_dot',m_tot_new*ones(length(vec_size),1),'gamma',sol.gamma);
Fuel = struct('properties',properties,'products',products);
%% PUMP FUEL

% Define parameters for the fuel pump
P_after_pump = 41.07; % Pressure after the pump (psia)
pump_eff_mech = 0.75; % Mechanical efficiency of the pump
pump_eff_iso = 0.98; % Isentropic efficiency of the pump
turb_eff_mech = 0.811; % Mechanical efficiency of the turbine
P_losses = P_after_pump - P; % Pressure losses
P_pre_pump = 1.72; % Pressure before the pump (psia)
T_pre_pump = 23.7; % Temperature before the pump (K)

h_pre_pump = -4401; % Enthalpy before the pump (23.7 K, 1.72 MPa)
%h_pre_pump = 48.692; % Enthalpy before the pump (23.7 K, 1.72 MPa)

h_post_pump = -3749; % Enthalpy after the pump (51.5 K, 41.07 MPa)
%h_post_pump = 697.76; % Enthalpy after the pump (51.5 K, 41.07 MPa)

H_m_dot = 70.3; % Mass flow rate of the fuel

% Calculate work done by the pump
W_pump = H_m_dot .* (h_post_pump - h_pre_pump);

% Calculate temperature after the turbine
T_pre_turb = sol.T;
T_after_turb = T_pre_turb - W_pump ./ (sol.cp .* m_tot_new' .* turb_eff_mech)

% Calculate temperature after isentropic process
T_after_iso = T_pre_turb + (T_after_turb - T_pre_turb) ./ pump_eff_iso;
Fuel.properties.T_before_main = T_after_turb;

% Calculate pressure at the end
P_end = sol.P .* (T_pre_turb ./ T_after_iso) .^ (sol.gamma ./ (1 - sol.gamma));
Fuel.properties.P_bar = P_end;
if length(O_F) > 1
    plot(O_F,T_pre_turb,'LineWidth',2)
    grid on;
    yline(986,'LineWidth',1.5,'LineStyle','-.','Color',"#A2142F");
    ylabel('T turbine'); xlabel('O/F')
    xline(0.869086,'LineWidth',2,'LineStyle',':','Color',"#EDB120")
    xticks(0.1:0.1:1.5)
    legend('T(O/F)','T = 986K','O/F = 0.8691')
elseif length(mass_area_ratio) > 1
    plot(mass_area_ratio,T_pre_turb,'LineWidth',2)
    grid on;
    yline(986,'LineWidth',2,'LineStyle','-.');
    ylabel('T turbine'); xlabel('m dot / A')
end

%% INJECTION PLATE FUEL LINE
% Definition of the cross-sectional area of the fuel preburner duct
A_preb_fuel = 0.0551;

% Calculation of the velocity of the fuel gas flow through the preburner
vel = sol.c * sol.Mach;
Q = vel * A_preb_fuel;

% Injection time and calculation of the required fuel volume in the preburner
t = 0.001434;
V_preb = Q * t;

% Calculation of the length of the preburner
L_preb = V_preb / A_preb_fuel;

% Definition of the diameter of a single injection hole, the area of the hole, and the discharge coefficient
%D_hole = 0.00173;
R_ox = 1.031/2*1e-3;
R_e_fuel = (4.93/2+0.58)*1e-3; R_i_fuel = 4.93/2*1e-3; %pag 25 NASA
A_1hole_ox = pi*R_ox^2;
A_1hole_fuel = pi*(R_e_fuel^2-R_i_fuel^2);
cd = 0.76;

% Calculation of the pressure difference between the system pressure and the pressure in the fuel and oxygen feed system
delta_P_fuel = -convpres(P, 'psi', 'Pa') + convpres(5310, 'psi', 'Pa');
delta_P_ox = -convpres(P, 'psi', 'Pa') + convpres(5636, 'psi', 'Pa');

% Calculation of the required area for the oxygen and fuel injection holes
A_hole_ox = O.m_dot / (cd * sqrt(2 * O.rho * delta_P_ox));
A_hole_fuel = H.m_dot / (cd * sqrt(2 * H.rho * delta_P_fuel));

% Calculation of the number of required injection holes and selection of the minimum between fuel and oxygen
N_ox = ceil(A_hole_ox / A_1hole_ox);
N_fuel = ceil(A_hole_fuel / A_1hole_fuel);
N = min(N_ox, N_fuel)

% Calculation of the new dimensions of the injection holes
A_1hole_ox = A_hole_ox / N;
A_1hole_fuel = A_hole_fuel / N;

R_1hole_ox = sqrt(A_1hole_ox / pi)
R_1hole_e_fuel = sqrt(A_1hole_fuel / pi + R_i_fuel^2)
% Calculation of the percentage of the total hole area compared to the preburner fuel area
A_perc = (A_hole_ox + A_hole_fuel) / A_preb_fuel * 100;

% Calculation of the velocities of the oxygen and fuel gases through the injection holes
u_ox = cd * sqrt(2 * delta_P_ox / O.rho);
u_fuel = cd * sqrt(2 * delta_P_fuel / H.rho);

% Calculation of the velocity ratio between fuel and oxygen
VR = u_fuel / u_ox;

% Calculation of a density-related parameter
J = (H.rho * u_fuel^2) / (O.rho * u_ox^2);

% Creation of a struct to store preburner fuel properties
preburner_fuel = struct('A', A_preb_fuel, 'V', V_preb, 'L', L_preb, ...
                        'N_holes', N * 2, 'D_hole_ox', D_1hole_ox, ...
                        'D_hole_fuel', D_1hole_fuel);
struct2table(preburner_fuel);
%% OX LINE

% Define parameters for the oxygen line
D_i = convlength(7.43,'in','m'); % Inner diameter (m)
L = convlength(4.25,'in','m'); % Length (m)
O_F = 0.6; % Oxidizer-to-fuel ratio
O_F = 0.1:0.1:10;
P = 4812; % Pressure (psia)

% Define properties for the oxygen and fuel streams
H = struct('m_dot',[],'T',[],'P',[],'rho',[],'h',[]);
O = H;
H.m_dot = 19.1; % Oxygen mass flow rate (kg/s)
H.T = 148.2; % Oxygen temperature (K)
H.P = 36.61; % Oxygen pressure (Mpa)
H.rho = 42.42; % Oxygen density (kg/m^3)
H.h = -2.503*2.016; % Oxygen enthalpy (kJ/mol)

O.m_dot = 11.3; % Fuel mass flow rate (kg/s)
O.T = 115.4; % Fuel temperature (K)
O.P = 39.53; % Fuel pressure (Mpa)
O.rho = 1116; % Fuel density (kg/m^3)
O.h = -0.3424*31.9988; % Fuel enthalpy (kJ/mol)

% Calculate total mass flow rate and area
m_tot = H.m_dot + O.m_dot;
A = D_i^2/4*pi;
mass_area_ratio = m_tot / A;

% Perform calculations using a function called 'cea_preburner_input'
[sol, products] = cea_preburner_input(H, O, O_F, P, mass_area_ratio);
CEA_results_ox = struct2table(sol) % Convert structure to table
struct2table(products); % Convert structure to table

% Define properties and store oxygen information
properties = struct('P_bar',sol.P,'T_before_main',[],'cp',sol.cp,'Mm',sol.Mm,'m_dot',m_tot*ones(length(O_F),1),'gamma',sol.gamma);
Ox = struct('properties',properties,'products',products);

%% PUMP OX

% Define parameters for the oxidizer pump
pump_eff_mech = 0.758; % Mechanical efficiency of the pump
boost_eff_mech = 0.718;
pump_eff_iso = 0.98; % Isentropic efficiency of the pump
turb_eff_mech = 0.746; % Mechanical efficiency of the turbine
P_after_boost = 27.89; % Pressure after the boost (psia)
boost_eff_mech = 0.718; % Mechanical efficiency of the boost
boost_eff_iso = 0.98; % Isentropic efficiency of the boost
P_pre_boost = 2.9; % Pressure before the boost (psia)
T_pre_boost = 93.7; % Temperature before the boost (K)
O_m_dot = 508.9; % Mass flow rate of the oxidizer
O_m_dot_boost = 41.73; % Mass flow rate of the boosted oxidizer

h_pre_pump = -397.35; % Enthalpy before the pump (93.7 K, 2.9 MPa)
%h_pre_pump = -4.0298;
h_post_pump = -363.65; % Enthalpy after the pump (104.3 K, 27.89 MPa)
%h_post_pump = -3.0415;

% Calculate work done by the pump and the boost
W_pump = O_m_dot * (h_post_pump - h_pre_pump);
h_pre_boost = -366.89; % Enthalpy before the boost (104.37 K, 26.96 MPa)
%h_pre_boost = -3.0543;
h_post_boost = -338.5; % Enthalpy after the boost (114.82 K, 48.06 MPa)
%h_post_boost = -2.1516;
W_boost = O_m_dot_boost .* (h_post_boost - h_pre_boost);

% Calculate temperature after the turbine
T_pre_turb = sol.T;
T_after_turb = T_pre_turb - (W_pump *pump_eff_mech + W_boost * boost_eff_mech) ./ (turb_eff_mech .* m_tot .* Ox.properties.cp)
Ox.properties.T_before_main = T_after_turb;

% Calculate temperature after isentropic process
T_after_iso = T_pre_turb + (T_after_turb - T_pre_turb) ./ pump_eff_iso;

% Calculate pressure at the end
P_end = sol.P .* (T_pre_turb ./ T_after_iso) .^ (sol.gamma ./ (1 - sol.gamma));
Ox.properties.P_bar = P_end;
if length(O_F) > 1
    plot(O_F,T_after_turb,'LineWidth',2)
end

%% INJECTION PLATE OX LINE
% Definition of the cross-sectional area of the oxidizer preburner duct
A_preb_ox = 0.028;

% Calculation of the velocity of the oxidizer gas flow through the preburner
vel = sol.c * sol.Mach;
Q = vel * A_preb_fuel;

% Injection time and calculation of the required oxidizer volume in the preburner
t = 0.001434;
V_preb = Q * t;

% Calculation of the length of the oxidizer preburner
L_preb = V_preb / A_preb_ox;

% Definition of the diameter of a single injection hole, the area of the hole, and the discharge coefficient
R_ox = 0.914*1e-3;
R_e_fuel = (4.7+0.648)*1e-3; R_i_fuel = 4.7*1e-3; %pag 25 NASA
A_1hole_ox = pi*R_ox^2;
A_1hole_fuel = pi*(R_e_fuel^2-R_i_fuel^2);
cd = 0.6;

% Calculation of the pressure difference between the system pressure and the pressure in the fuel and oxidizer feed system
delta_P_fuel = -convpres(P, 'psi', 'Pa') + convpres(5310, 'psi', 'Pa');
delta_P_ox = -convpres(P, 'psi', 'Pa') + convpres(5734, 'psi', 'Pa');

% Calculation of the required area for the oxidizer and fuel injection holes
A_hole_ox = O.m_dot / (cd * sqrt(2 * O.rho * delta_P_ox));
A_hole_fuel = H.m_dot / (cd * sqrt(2 * H.rho * delta_P_fuel));

% Calculation of the number of required injection holes and selection of the minimum between fuel and oxidizer
N_ox = ceil(A_hole_ox / A_1hole_ox);
N_fuel = ceil(A_hole_fuel / A_1hole_fuel);
N = min(N_ox, N_fuel)

% Calculation of the new dimensions of the injection holes
A_1hole_ox = A_hole_ox / N;
A_1hole_fuel = A_hole_fuel / N;
D_1hole_fuel = sqrt(4 * A_1hole_fuel / pi);
D_1hole_ox = sqrt(4 * A_1hole_ox / pi);

% Calculation of the percentage of the total hole area compared to the preburner fuel area
A_perc = (A_hole_ox + A_hole_fuel) / A_preb_fuel * 100;

% Calculation of the velocity ratio between fuel and oxidizer
VR = u_fuel / u_ox;

% Calculation of a density-related parameter
J = (H.rho * u_fuel^2) / (O.rho * u_ox^2);

% Creation of a struct to store preburner oxidizer properties
preburner_ox = struct('A', A_preb_fuel, 'V', V_preb, 'L', L_preb, ...
                        'N_holes', N * 2, 'D_hole_ox', D_1hole_ox, ...
                        'D_hole_fuel', D_1hole_fuel);
struct2table(preburner_ox);

%% Results

% Display results for the fuel line
fprintf('\t\t- - - - - - - - - - - - - - - - - - - - - - - - - -\n');
fprintf('\t\t- - - - - - - - - - - - - - - - - - - - - - - - - -\n');
Fuel_line_products = struct2table(Fuel.products)
Fuel_line_properties = struct2table(Fuel.properties)
fprintf('\t\t- - - - - - - - - - - - - - - - - - - - - - - - - -\n');
fprintf('\t\t- - - - - - - - - - - - - - - - - - - - - - - - - -\n');
% Display results for the oxidizer line
Ox_line_products = struct2table(Ox.products)
Ox_line_properties = struct2table(Ox.properties)
fprintf('\t\t- - - - - - - - - - - - - - - - - - - - - - - - - -\n');
fprintf('\t\t- - - - - - - - - - - - - - - - - - - - - - - - - -\n');
% Display geometric results
preburner_geometry_fuel = struct2table(preburner_fuel)
reburner_geometry_ox = struct2table(preburner_fuel)


%% MAIN
OX_m_dot_H2 = Ox.properties.m_dot*Ox.products.Mass_fraction(1)
OX_m_dot_H2O = Ox.properties.m_dot*Ox.products.Mass_fraction(2)

Fuel_m_dot_H2 = Fuel.properties.m_dot*Fuel.products.Mass_fraction(1)
Fuel_m_dot_H2O = Fuel.properties.m_dot*Fuel.products.Mass_fraction(2)

coolant_m_dot = 5+7.7;

O2_main = struct('m_dot',[],'T',[],'P',[],'rho',[],'h',[],'Mf',[]);
O2_main.m_dot = 380.6;
O2_main.T = 104.3;
O2_main.P = 24.03; % Pressure (Mpa)
O2_main.rho = 1128.4;
O2_main.h = -0.3685*31.9988;


fprintf('\t\t- - - - - - - - - - - - - - - - - - - - - - - - - -\n');
fprintf('\t\t- - - - - - - - - - - - - - - - - - - - - - - - - -\n');

m_tot_fuel = (OX_m_dot_H2 + Fuel_m_dot_H2 + coolant_m_dot);
m_tot_ox = O2_main.m_dot;
m_tot = m_tot_fuel+m_tot_ox;

mass_fraction_coolant = coolant_m_dot / m_tot;
mass_fraction_fuel_H2 = Fuel_m_dot_H2 / m_tot;
mass_fraction_ox_H2 = OX_m_dot_H2 / m_tot;
mass_fraction_fuel_H2O = Fuel_m_dot_H2O / m_tot;
mass_fraction_ox_H2O = OX_m_dot_H2O / m_tot;
mass_fraction_main_O2 = O2_main.m_dot / m_tot;

H2_cool = struct('Mf',mass_fraction_coolant,'T',255.4,'P',[]);
H2_fuel = struct('Mf',mass_fraction_fuel_H2,'T',Fuel.properties.T_before_main,'P',Fuel.properties.P_bar);
H2_ox = struct('Mf',mass_fraction_ox_H2,'T',Ox.properties.T_before_main,'P',Ox.properties.P_bar);
O2_main.Mf = mass_fraction_main_O2;

ae_at = 69;
[sol,prod] = cea_main_input(H2_cool,H2_fuel,H2_ox,O2_main,2871,ae_at);
%% plot
subplot(2,2,1)
plot(ae_at,sol.T,'LineWidth',2)
title('T - ae/at')
subplot(2,2,2)
plot(ae_at,sol.cf,'LineWidth',2)
title('cf - ae/at')
subplot(2,2,3)
plot(ae_at,sol.Isp_s,'LineWidth',2)
title('Isp_sea - ae/at')
subplot(2,2,4)
plot(ae_at,sol.Isp_v,'LineWidth',2)
title('Isp_vac - ae/at')

