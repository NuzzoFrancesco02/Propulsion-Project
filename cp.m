%% cp properties
function [gamma_medio] = cp(T0,T1)

cp_file = fopen('fluid-3.txt','r');
str = fscanf(cp_file,'%c');
% Cerca la posizione del testo "COMB END"
posizione_mass_fraction = strfind(str, 'Temperature');
 
% Trova il testo relativo alla sezione "COMB END"
sezione_mass_fraction = extractAfter(str, posizione_mass_fraction(1));

% Divide la sezione "COMB END" in righe
righe_comb_end = strsplit(sezione_mass_fraction, '\n');
T = []; Cp = []; Cv = [];
for i = 2 : numel(righe_comb_end)-1
    riga = strsplit(righe_comb_end{i});
    T = [T; str2double(riga(1))];
    Cv = [Cv; str2double(riga(8))];
    Cp = [Cp; str2double(riga(9))];
end
gamma_medio = Cp./Cv;
if nargin == 1
    Cp = interp1(T,Cp,T0);
    gamma_medio = interp1(T,gamma_medio,T0);
else
    T_dis = T0:0.001:T1;
    %plot(T,Cp,'o','LineWidth',2); hold on
    Cp = interp1(T,Cp,T_dis);
    gamma_medio = interp1(T,gamma_medio,T_dis);
    %plot(T_dis,Cp,'LineWidth',2)
    %figure(); plot(T,gamma)
    gamma_medio = 1/(gamma_medio(end)-gamma_medio(1))*trapz(gamma_medio);
end

%%

