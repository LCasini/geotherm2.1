function[thermo_Tlow,thermo_Thigh,thermo_Plow,thermo_Phigh,thermo_n] = thermobarometry(Ma,nt,dt,t0,tfin)
% the function thermobarometry is used to display geologic age-depth-T
% constraints provided by the user from filling thermobarometry.xls file
%==========================================================================
data = xlsread('thermobarometry','A12:G100');     % load user-defined data
n = numel(data(:,1));
thermo_n = n;

% set thermobarometric constraints to be displayed-------------------------
%
% initialize vectors, all n-by-1 vectors where n is the number of
% thermobarometric constraints
Tav = data(:,1);
Tsd = data(:,2);
Pav = data(:,3);
Psd = data(:,4);
age_av = data(:,5);
age_sd = data(:,6);


age_max = age_av+age_sd;
age_min = age_av-age_sd;

% loop through the thermobarometric points---------------------------------
for i = 1:n
    if age_max(i) > t0
    age_max(i) = t0;
    end
    if age_min(i) < tfin
    age_min(i) = tfin;
    end
end
%--------------------------------------------------------------------------
%
% initialize the number of timestep of thermobarometric n-points-----------

init_dtn = zeros(1,n);
fin_dtn = zeros(1,n);

% initialize matrices, all n-by-nt-----------------------------------------

thermo_Tlow = zeros(n,nt);
thermo_Thigh = zeros(n,nt);
thermo_Plow = zeros(n,nt);
thermo_Phigh = zeros(n,nt);

Tlow = Tav-Tsd;              % temperature [K]
Thigh = Tav+Tsd;             % temperature [K]
Plow = (Pav-Psd)*36;         % pressure [km]
Phigh = (Pav+Psd)*36;        % pressure [km]

% find initial timestep of constraints

for i = 1:n
    % loop through the n-points
    % ensures the constraints are within experimental timestep loop--------
    init_dtn(i) = round(1+(t0-age_max(i))/(dt/Ma));
    fin_dtn(i) = round((t0-age_min(i))/(dt/Ma));            
    
    for n = 1:nt
        thermo_Tlow(i,init_dtn(i):fin_dtn(i)) = Tlow(i);
        thermo_Thigh(i,init_dtn(i):fin_dtn(i)) = Thigh(i);
        thermo_Plow(i,init_dtn(i):fin_dtn(i)) = Plow(i);
        thermo_Phigh(i,init_dtn(i):fin_dtn(i)) = Phigh(i);
    end
    
end

end