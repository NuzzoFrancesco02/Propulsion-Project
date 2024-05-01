%% 
clear
clc

g = 9.81; % m^2/s
O_F = 6.032; % rapporto di miscela
t = 480; % tempo di burning
Pcc = 3009 * 6894.76; %Pa   pressione camera di combustione
Tcc = 3776.68; % K  da CEA
Tsl = 376600 * 4.448222; %N    spinta al livello di mare
Tvac = 470800 * 4.448222; %N    spinta nel vuoto
gamma = 1.1441; % da CEA
Mm = 13.411; % da CEA
R = 8314/Mm; % da CEA
Pe = 0.56772 * 1e5; % Pa ???????
I_sl = 4558.1/g; % da CEA
I_vac = 4741.7/g; % da CEA
eps = 69; 

Mcc = 0.2; % ? ? ? ? ? 
alpha_conv = deg2rad(45); %°
alpha_div = deg2rad(15); %°
%%
GAMMA = sqrt(gamma*(2/(gamma+1))^((gamma+1)/(gamma-1)));
m_dot = Tsl/2/(g*I_sl);
c_star = sqrt(R*Tcc)/GAMMA;
c_t = I_sl*g/c_star;
u_e = sqrt(2*g/(g-1)*R*Tcc*(1-(Pe/Pcc)^((gamma-1)/(gamma))))

%% Dimensionamento ugello
A_t = Tsl / 2 / (Pcc * c_t);
D_t = sqrt(4 * A_t / pi);
A_e = eps*A_t;
D_e = sqrt(4 * A_e / pi);
u_c = Mcc * sqrt(gamma * R * Tcc);
rho_cc = Pcc/(R * Tcc);
A_cc = m_dot / (rho_cc * u_c);
D_cc = sqrt(4 * A_cc / pi);
L_conv = (D_cc - D_t) / (2 * atan(alpha_conv));
L_div = (D_e - D_t) / (2*atan(alpha_div));
L_noz = L_div + L_conv;

lambda = (1 + cos(alpha_div) ) / 2;
L_tot = 14*0.3048;
L_cc = L_tot-L_noz;
V_cc = A_cc*L_cc;
L_star = V_cc/L_cc
%% CALCOLO PORTATE
clear 
clc
alpha_stech = 6.032;

m_o2_l = 67 * 2.205; % kg/s    portata ox sinistra
m_o2_r = 25 * 2.205; % kg/s     portata ox destra
m_h2_l = 65 * 2.205; % kg/s     portata h2 sinistra
m_h2_r = 42 * 2.205; % kg/s     portata h2 destra
m_h2_r_burn = m_o2_r / alpha_stech;
m_h2_l_burn = m_o2_l / alpha_stech;
m_h2_refill_l = 17 * 2.205; 
m_h2_cons_l = m_h2_l - m_h2_l_burn
m_h2_cons_r = m_h2_r - m_h2_r_burn

% bilancio
p_min = 3095 * 6894.76;
cp_h2 = 14.304;
cp_pre_l = 8.1056;
cp_post_l = 8.6837;
T_tot = (m_h2_refill_l*cp_h2*(273.15)+(m_h2_l+m_o2_l)*1130.39*cp_pre_l)/((m_h2_refill_l+m_h2_l+m_o2_l)*cp_post_l)
%% propr reag
T_H2 = 150; %K
Mm_H2 = 2.008;
R_H2 = 8.314/Mm_H2;
h_H2 = 1.870079;
u_H2 = h_H2-R_H2*T_H2

