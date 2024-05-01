function [sol, products] = cea_preburner_input(struct1, struct2, O_F, P, mass_area_ratio)
    % Initialize solution struct
    sol = struct('P',[],'T',[],'rho',[],'h',[],'Mm',[],'cp',[],'gamma',[]);
    
    % Loop through each oxidizer to fuel ratio
    for i = 1 : length(O_F)
        % Open input file and read content
        problem = 'rocket';
        file_input = fopen("input_rocket.txt",'r');
        str = fscanf(file_input,'%c');
        pat = "%" + lettersPattern(1);
        data_types = extract(str,pat);
        
        % Change directory to CEA folder
        currentDir = pwd;
        cd("CEA/");
        
        % Create and open input file for CEA
        file_case = fopen("input_problem.inp",'w');
        fprintf(file_case,str,problem,O_F(i),mass_area_ratio,P,struct1.T,struct1.h,struct2.T,struct2.h);
        
        filename = 'input_problem';  % Nome del file da passare a FCEA2

        % Execute FCEA2
        setenv('PATH', [getenv('PATH') ':/usr/local/bin:/usr/bin:/bin/']);
        fprintf(fopen("FCEA2_input.txt",'w'),'input_problem');
        [~,~] = system('./FCEA2 < FCEA2_input.txt');
        delete('FCEA2_input.txt');
        fclose(file_case);
        unsetenv('PATH');
        
        % Open output file and read content
        file_output = fopen('input_problem.out','r');
        str = fscanf(file_output,'%c');
        
        % Find "COMB END" section
        posizione_comb_end = strfind(str, 'COMB END');
        sezione_comb_end = extractAfter(str, posizione_comb_end(1));
        fine_sezione_mass_fraction = strfind(sezione_comb_end, 'PERFORMANCE PARAMETERS');
        sezione_comb_end = sezione_comb_end(1:fine_sezione_mass_fraction(1));
        
        % Extract relevant data from "COMB END" section
        righe_comb_end = strsplit(sezione_comb_end, '\n');
        colonna2 = [];
        for j = 2 : numel(righe_comb_end)
            riga = righe_comb_end{j};
            valori = strsplit(riga);
            if numel(valori) >= 3
                colonna2 = [colonna2; str2double(valori{end-1})];
            end
        end
        
        % Store data in solution struct
        sol.P = [sol.P; colonna2(2)];
        sol.T = [sol.T; colonna2(3)];
        sol.rho = [sol.rho; colonna2(4)];
        sol.h = [sol.h; colonna2(5)];
        sol.Mm = [sol.Mm; colonna2(9)];
        sol.cp = [sol.cp; colonna2(12)];
        sol.gamma = [sol.gamma; colonna2(13)];

        % Find "MASS FRACTIONS" section
        posizione_mass_fraction = strfind(str, 'MASS FRACTIONS');
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
        if length(O_F) == 1
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
        if length(O_F) > 1
            products.name = [products.name; '----'];
            products.Mass_fraction = [products.Mass_fraction; ' '];
        end
        products.Mass_fraction = [products.Mass_fraction; var_M_fraction];
        cd(currentDir)
    end
    
    % Return solution and products
    delete("input_problem.inp"); delete("input_problem.out");

end
