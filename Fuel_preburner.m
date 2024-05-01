%% FUEL LINE
clear
clc

D_i = convlength(10.43,'in','m'); %m
L = convlength(4.37,'in','m'); %m
O_F = 0.86905;
%O_F = 0.1:0.1:10;
P = 4793; %psia


H = struct('m_dot',[],'T',[],'P',[],'rho',[],'h',[]);
O = H;
H.m_dot = 34.9; %kg/s
H.T = 148.2; %K 
H.P = 36.61; %Mpa
H.rho = 42.42; %kg/m^3
H.h = -2.503*2.016; %kJ/mol

O.m_dot = 30.4; %kg/s
O.T = 115.4; %K 
O.P = 38.86; %Mpa
O.rho = 1114.6; %kg/m^3
O.h = -0.3427*31.9988; %kJ/mol

m_tot = H.m_dot+O.m_dot;
A = D_i^2/4*pi;
mass_area_ratio = m_tot/A;
[sol,products] = cea_preburner_input(H,O,O_F,P,mass_area_ratio);
CEA_results_fuel = struct2table(sol)
struct2table(products);
% MIX
H_cool = struct('m_dot',[],'T',[],'P',[],'rho',[],'h',[]);
H_cool.m_dot = 2.7;
H_cool.T = 51.5;
H_cool.P = 41.07; %Mpa
H_cool.rho = 78.82;
H_cool.h = -7.539048;
m_tot_new = m_tot+H_cool.m_dot;
perc_H_cool = H_cool.m_dot/m_tot_new;
perc_H = m_tot*products.Mass_fraction(1)/m_tot_new;
perc_H20 = m_tot*products.Mass_fraction(2)/m_tot_new;
properties = struct('P_bar',sol.P,'T_before_main',[],'cp',sol.cp,'Mm',sol.Mm,'m_dot',m_tot_new,'gamma',sol.gamma);
Fuel = struct('properties',properties,'products',products);
%% PUMP FUEL
P_after_pump = 41.07; %psia; 
pump_eff_mech = 0.75;
pump_eff_iso = 0.98;
turb_eff_mech = 0.811;
P_losses = P_after_pump-P;
P_pre_pump = 1.72; %psia
T_pre_pump = 23.7; %K

h_pre_pump = -4401; % 23.7 K, 1.72 MPa
h_post_pump = -3749; % 51.5 K, 41.07 MPa
H_m_dot = 70.3;

W_pump = H_m_dot .* (h_post_pump-h_pre_pump);
T_pre_turb = sol.T;
T_after_turb = sol.T - W_pump./(sol.cp*m_tot_new*turb_eff_mech);

T_after_iso = T_pre_turb + (T_after_turb-T_pre_turb)/pump_eff_iso;
Fuel.properties.T_before_main = T_after_turb;
P_end = sol.P*(T_pre_turb/T_after_iso)^((sol.gamma)/(1-sol.gamma));
Fuel.properties.P_bar = P_end;
if length(O_F) > 1
    plot(O_F,T_after_turb,'LineWidth',2)
end
%% OX LINE
D_i = convlength(7.43,'in','m'); %m
L = convlength(4.25,'in','m'); %m
O_F = 0.6103;
P = 4812; %psia


H = struct('m_dot',[],'T',[],'P',[],'rho',[],'h',[]);
O = H;
H.m_dot = 19.5; %kg/s
H.T = 148.2; %K 
H.P = 36.61; %Mpa
H.rho = 42.42; %kg/m^3
H.h = -2.503*2.016; %kJ/mol

O.m_dot = 11.3; %kg/s
O.T = 115.4; %K 
O.P = 39.53; %Mpa
O.rho = 1116; %kg/m^3
O.h = -0.3424*31.9988; %kJ/mol

m_tot = H.m_dot+O.m_dot;
A = D_i^2/4*pi;
mass_area_ratio = m_tot/A;
[sol,products] = cea_preburner_input(H,O,O_F,P,mass_area_ratio);
CEA_results_ox = struct2table(sol)
struct2table(products);
properties = struct('P_bar',sol.P,'T_before_main',[],'cp',sol.cp,'Mm',sol.Mm,'m_dot',m_tot,'gamma',sol.gamma);
Ox = struct('properties',properties,'products',products);


%% PUMP OX
P_after_boost = 27.89; %psia; 
boost_eff_mech = 0.718;
boost_eff_iso = 0.98;
P_pre_boost = 2.9; %psia
T_pre_boost = 93.7; %K
O_m_dot = 508.9;
O_m_dot_boost = 41.73;
P_after_pump = 27.89; %psia; 
pump_eff_mech = 0.718;
pump_eff_iso = 0.98;
turb_eff_mech = 0.811;
P_losses = P_after_pump-P;
P_pre_pump = 2.9; %psia
T_pre_pump = 93.7; %K

h_pre_pump = -397.35; % 93.7 K, 2.9 MPa
h_post_pump = -363.65; % 104.3 K, 27.89 MPa
% turbopump
W_pump = O_m_dot * (h_post_pump-h_pre_pump);
% boost
h_pre_boost = -366.89; % 104.37 K, 26.96 MPa
h_post_boost = -338.5; % 114.82 K, 48.06 MPa
W_boost = O_m_dot_boost .* (h_post_boost-h_pre_boost);
T_pre_turb = sol.T;
T_after_turb = sol.T - (W_pump+W_boost)./(turb_eff_mech.*m_tot.*Ox.properties.cp);
Ox.properties.T_before_main = T_after_turb;

T_after_iso = T_pre_turb + (T_after_turb-T_pre_turb)/pump_eff_iso;
P_end = sol.P*(T_pre_turb/T_after_iso)^((sol.gamma)/(1-sol.gamma));
Fuel.properties.P_bar = P_end;
if length(O_F) > 1
    plot(O_F,T_after_turb,'LineWidth',2)
end

%% results
fprintf('\t\t- - - - - - - - - - - - - - - - - - - - - - - - - -\n');
fprintf('\t\t- - - - - - - - - - - - - - - - - - - - - - - - - -');
Fuel_line_products = struct2table(Fuel.products)
Fuel_line_properties = struct2table(Fuel.properties)
fprintf('\t\t- - - - - - - - - - - - - - - - - - - - - - - - - -\n');
fprintf('\t\t- - - - - - - - - - - - - - - - - - - - - - - - - -');
Ox_line_products = struct2table(Ox.products)
Ox_line_properties = struct2table(Ox.properties)

%%
%% FUEL LINE

% Clear command window and workspace
clear
clc

% Define parameters for the fuel line
D_i = convlength(10.43,'in','m'); % Inner diameter (m)
L = convlength(4.37,'in','m'); % Length (m)
O_F = 0.86905; % Oxidizer-to-fuel ratio
P = 4793; % Pressure (psia)

% Define properties for the fuel and oxidizer streams
H = struct('m_dot',[],'T',[],'P',[],'rho',[],'h',[]);
O = H;
H.m_dot = 34.9; % Fuel mass flow rate (kg/s)
H.T = 148.2; % Fuel temperature (K)
H.P = 36.61; % Fuel pressure (Mpa)
H.rho = 42.42; % Fuel density (kg/m^3)
H.h = -2.503*2.016; % Fuel enthalpy (kJ/mol)

O.m_dot = 30.4; % Oxidizer mass flow rate (kg/s)
O.T = 115.4; % Oxidizer temperature (K)
O.P = 38.86; % Oxidizer pressure (Mpa)
O.rho = 1114.6; % Oxidizer density (kg/m^3)
O.h = -0.3427*31.9988; % Oxidizer enthalpy (kJ/mol)

% Calculate total mass flow rate and area
m_tot = H.m_dot + O.m_dot;
A = D_i^2/4*pi;
mass_area_ratio = m_tot / A;

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
perc_H_cool = H_cool.m_dot / m_tot_new;
perc_H = m_tot * products.Mass_fraction(1) / m_tot_new;
perc_H20 = m_tot * products.Mass_fraction(2) / m_tot_new;

% Define properties and store fuel information
properties = struct('P_bar',sol.P,'T_before_main',[],'cp',sol.cp,'Mm',sol.Mm,'m_dot',m_tot_new*ones(length(O_F),1),'gamma',sol.gamma);
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
h_post_pump = -3749; % Enthalpy after the pump (51.5 K, 41.07 MPa)
H_m_dot = 70.3; % Mass flow rate of the fuel

% Calculate work done by the pump
W_pump = H_m_dot .* (h_post_pump - h_pre_pump);

% Calculate temperature after the turbine
T_pre_turb = sol.T;
T_after_turb = T_pre_turb - W_pump ./ (sol.cp * m_tot_new * turb_eff_mech);

% Calculate temperature after isentropic process
T_after_iso = T_pre_turb + (T_after_turb - T_pre_turb) ./ pump_eff_iso;
Fuel.properties.T_before_main = T_after_turb;

% Calculate pressure at the end
P_end = sol.P .* (T_pre_turb ./ T_after_iso) .^ (sol.gamma ./ (1 - sol.gamma));
Fuel.properties.P_bar = P_end;
if length(O_F) > 1
    plot(O_F,T_after_turb,'LineWidth',2)
end
%% OX LINE

% Define parameters for the oxygen line
D_i = convlength(7.43,'in','m'); % Inner diameter (m)
L = convlength(4.25,'in','m'); % Length (m)
%O_F = 0.6103; % Oxidizer-to-fuel ratio
O_F = 0.1:0.1:10;
P = 4812; % Pressure (psia)

% Define properties for the oxygen and fuel streams
H = struct('m_dot',[],'T',[],'P',[],'rho',[],'h',[]);
O = H;
H.m_dot = 19.5; % Oxygen mass flow rate (kg/s)
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
P_after_boost = 27.89; % Pressure after the boost (psia)
boost_eff_mech = 0.718; % Mechanical efficiency of the boost
boost_eff_iso = 0.98; % Isentropic efficiency of the boost
P_pre_boost = 2.9; % Pressure before the boost (psia)
T_pre_boost = 93.7; % Temperature before the boost (K)
O_m_dot = 508.9; % Mass flow rate of the oxidizer
O_m_dot_boost = 41.73; % Mass flow rate of the boosted oxidizer

h_pre_pump = -397.35; % Enthalpy before the pump (93.7 K, 2.9 MPa)
h_post_pump = -363.65; % Enthalpy after the pump (104.3 K, 27.89 MPa)

% Calculate work done by the pump and the boost
W_pump = O_m_dot * (h_post_pump - h_pre_pump);
h_pre_boost = -366.89; % Enthalpy before the boost (104.37 K, 26.96 MPa)
h_post_boost = -338.5; % Enthalpy after the boost (114.82 K, 48.06 MPa)
W_boost = O_m_dot_boost .* (h_post_boost - h_pre_boost);

% Calculate temperature after the turbine
T_pre_turb = sol.T;
T_after_turb = T_pre_turb - (W_pump + W_boost) ./ (turb_eff_mech .* m_tot .* Ox.properties.cp);
Ox.properties.T_before_main = T_after_turb;

% Calculate temperature after isentropic process
T_after_iso = T_pre_turb + (T_after_turb - T_pre_turb) ./ pump_eff_iso;

% Calculate pressure at the end
P_end = sol.P .* (T_pre_turb ./ T_after_iso) .^ (sol.gamma ./ (1 - sol.gamma));
Ox.properties.P_bar = P_end;
if length(O_F) > 1
    plot(O_F,T_after_turb,'LineWidth',2)
end
%% Results

% Display results for the fuel line
fprintf('\t\t- - - - - - - - - - - - - - - - - - - - - - - - - -\n');
fprintf('\t\t- - - - - - - - - - - - - - - - - - - - - - - - - -');
Fuel_line_products = struct2table(Fuel.products)
Fuel_line_properties = struct2table(Fuel.properties)
fprintf('\t\t- - - - - - - - - - - - - - - - - - - - - - - - - -\n');
fprintf('\t\t- - - - - - - - - - - - - - - - - - - - - - - - - -');

% Display results for the oxidizer line
Ox_line_products = struct2table(Ox.products)
Ox_line_properties = struct2table(Ox.properties)






