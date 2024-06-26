function [sol, products] = cea_main_input(H_cool, H_fuel, H_ox, O2_main,P,ae_at)
    % Initialize solution struct
    sol = struct('P',[],'T',[],'rho',[],'h',[],'Mm',[],'cp',[],'gamma',[],'Isp_s',[],'Isp_v',[],'cf',[],...
        'sound',[],'Mach',[]);
    
    % Loop through each oxidizer to fuel ratio
    for i = 1 : length(ae_at)
        % Open input file and read content
        case_name = 'main';
        file_input = fopen(['input_',case_name,'.txt'],'r');
        str = fscanf(file_input,'%c');
        pat = "%" + lettersPattern(1);
        data_types = extract(str,pat);
        
        % Change directory to CEA folder
        currentDir = pwd;
        cd("CEA/");
        
        % Create and open input file for CEA
        file_case = fopen(['input_',case_name,'.inp'],'w');
        fprintf(file_case,str,'rocket',P,ae_at(i),H_cool.Mf,H_cool.T,H_fuel.Mf,H_fuel.T,H_ox.Mf,H_fuel.T,...
            O2_main.Mf,O2_main.T,O2_main.h);
        
        filename = ['input_',case_name];  % Nome del file da passare a FCEA2

        % Execute FCEA2
        setenv('PATH', [getenv('PATH') ':/usr/local/bin:/usr/bin:/bin/']);
        fprintf(fopen("FCEA2_input.txt",'w'),filename);
        [~,~] = system('./FCEA2 < FCEA2_input.txt');
        delete('FCEA2_input.txt');
        fclose(file_case);
        unsetenv('PATH');
        
        % Open output file and read content
        file_output = fopen(['input_',case_name,'.out'],'r');
        str = fscanf(file_output,'%c');
        
        % Find "COMB END" section
        posizione_comb_end = strfind(str, 'AFTER POINT 2');
        sezione_comb_end = extractAfter(str, posizione_comb_end(1));
        posizione_comb_end = strfind(sezione_comb_end,'CHAMBER');
        sezione_comb_end = extractAfter(sezione_comb_end, posizione_comb_end(1));
        fine_sezione_mass_fraction = strfind(sezione_comb_end, 'MASS FRACTIONS');
        sezione_comb_end = sezione_comb_end(1:fine_sezione_mass_fraction(1));
        
        % Extract relevant data from "COMB END" section
        righe_comb_end = strsplit(sezione_comb_end, '\n');
        exit = [];
        for j = 2 : numel(righe_comb_end)
            riga = righe_comb_end{j};
            valori = strsplit(riga);
            if numel(valori) >= 3
                str_val = valori{end};
                if ~isnan(strfind(str_val(end-2:end),'-'))
                    str_val = replace(str_val,'-','e-')
                end
                exit = [exit; str2double(str_val)];
            end
        end
        
        % Store data in solution struct
        sol.P = [sol.P; exit(2)];
        sol.T = [sol.T; exit(3)];
        sol.rho = [sol.rho; exit(4)];
        sol.h = [sol.h; exit(5)];
        sol.Mm = [sol.Mm; exit(9)];
        sol.cp = [sol.cp; exit(10)];
        sol.gamma = [sol.gamma; exit(11)];
        sol.sound = [sol.sound; exit(12)];
        sol.Mach = [sol.Mach; exit(13)];
        sol.cf = [sol.cf; exit(17)];
        sol.Isp_s = [sol.Isp_s; exit(19)];
        sol.Isp_v = [sol.Isp_v; exit(18)];
        
        % Find "MASS FRACTIONS" section
        posizione_mass_fraction = strfind(str, 'FRACTIONS');
        sezione_mass_fraction = extractAfter(str, posizione_mass_fraction(1));
        righe_mass_fraction = strsplit(sezione_mass_fraction, '\n');
        fine_sezione_mass_fraction = strfind(sezione_mass_fraction, '* THERMODYNAMIC');
        sezione_mass_fraction = sezione_mass_fraction(1:fine_sezione_mass_fraction(1));
        righe_mass_fraction = strsplit(sezione_mass_fraction, '\n');
        
        % Extract mass fractions data
        number_products = numel(righe_mass_fraction) - 2;
        if i == 1
            products = struct('name',[],'Mass_fraction',[]);
        end
        var_M_fraction = zeros(number_products,1);
        var_name = {};
        flag = 1;
        if length(ae_at) == 1
            flag = 0;
        end
        for j = 1 : number_products + 1
            riga = righe_mass_fraction{j};
            valori = strsplit(riga);
            if numel(valori) >= 3
                var_name = [var_name; valori{2}];
                var_M_fraction(j-1) = str2double(valori{end-1});
            end
        end
        products.name = [products.name; var_name(:,end)];
        if length(ae_at) > 1
            products.name = [products.name; '----'];
            products.Mass_fraction = [products.Mass_fraction; ' '];
        end
        products.Mass_fraction = [products.Mass_fraction; var_M_fraction];
        if i ~= length(ae_at)
            delete(['input_',case_name,'.inp']); delete(['input_',case_name,'.out']);
        end
        cd(currentDir)
    end
    
    % Return solution and products
   
end
