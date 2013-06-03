function [shear_heating,stress] = VSH(nz,t0,dt,nzlit,T,timestep,nzcrust,melt,z,viscosity)
% the function shear_heating calculates the volumetric heat production rate 
% due to dissipation of mechanical work; deformation is assumed to be 
% simple shear, so VSH = strain-rate*stress (Turcotte & Schubert 2002).  
%==========================================================================

data = xlsread('data','Q10:S19');
if isempty(data) == 0
    m = numel(data(:,1));
    year = 365.25*3600*24;                       % seconds per year [s]
    Ma = 1e6*year;                               % seconds per Ma [s]

    z = wrev(abs(z));                            % reshape z vector

    % initialize start and end timesteps of ductile events---------------------

    starting = zeros(m,1);
    ending = zeros(m,1);
    strain_rate = zeros(m,1);
    shear_heating = zeros(1,nz);                  % initialize VSH vector
    stress = zeros(1,nz);                         % initialize stress vector

% calculate shear heating--------------------------------------------------

for i = 1:m
    
    starting(i) = round(1+(t0-data(i,1))/(dt/Ma));
    ending(i) = round((t0-data(i,2))/(dt/Ma));
    strain_rate(i) = 10.^(data(i,3)); 
    
    % evaluate shear heating at the current timestep-----------------------
    if starting(i) <= timestep && ending(i) >= timestep
        
        % set the limits of ductile deformation, find nz points------------
        
        check = T-573.15;
        for n = 1:nz
            if check(n) < 0
                check(n) = NaN;
            end
        end
        
        z_number = sum(isnan(check));            % find depth at which 
                                                 % ductile deformation
                                                 % begins
        
        % evaluate stress through the crust in case of shearing------------
        
        for j = 1:nz 
            if j <= z_number
                stress(j) = 0.023409*z(j);          % Byerlee 
                shear_heating(j) = 0;
            elseif j > z_number && j <= nzcrust
                 % wet quartzite Gleason & Tullis (1995)        
                 if melt(j) <= 0.1  
                     % stress in solid crustal rocks
                     A = 1.1e-4;                    % constant [MPa^-n]
                     Q = 223*1e3;                   % enthalpy [J*mol^-1]
                     R = 8.314;                     % gas constant
                     stress(j) = (strain_rate(i)/(A*exp(-Q/(R*T(j))))).^(1/4);
                     shear_heating(j) = strain_rate(i).*(stress(j)*1e6); % W/m3
                 else   
                     % stress in molten crustal rocks
                     stress(j) = (strain_rate(i).*viscosity(j))./1e6; 
                     shear_heating(j) = strain_rate(i).*stress(j); % W/m3
                 end
            elseif j >= nzcrust && j <= nzlit
                if melt(j) <= 0.1
                    % stress in solid mantle rocks
                    % Olivine Hirth & Kohlstedt (2003)
                    A = 1.1*10^5;                    % constant [MPa^-n]
                    Q = 530*1e3;                     % enthalpy [J*mol^-1]
                    R = 8.314;                       % gas constant
                    stress(j) = (strain_rate(i)/(A*exp(-Q/(R*T(j))))).^(1/3.5);
                    shear_heating(j) = strain_rate(i)*(stress(j)*1e6); % W/m3
                else
                    % stress in molten mantle rocks
                    stress(j) = (strain_rate(i)*viscosity(j))./1e6;
                    shear_heating(j) = strain_rate(i).*stress(j); % W/m3
                end
            else
            end 
        end
    end 
end
else
    shear_heating = zeros(1,nz);
    stress = zeros(1,nz);
end
end
 