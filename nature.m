function [t0,tfin,L,dz,H,Hsd,Hsd_minus,Moho,lc,nzcrust,density,variation] = nature(nt,z0,nzlit,nzsublit,heat_rate,heat_sd,heat_sd_minus)
% the function nature calculates the rate of vertical displacement in
% the 1D model lithosphere. Rates are derived from user-defined model
% parameters loaded from an external excel file (data.xls). IMPORTANT: 
% the data.xls file must be in the same folder as the matlab scripts. 
%--------------------------------------------------------------------------
% NOTE that before runninng numerical simulations, the mainscript prompt
% the user to fill in the data.m file with specific user-defined
% values.
%--------------------------------------------------------------------------
% outputs:
% t0 = scalar, stores the start time of experiment [Ma]
% tfin = scalar, stores the end of experiment [Ma]
% dt = scalar, stores the dimension of timestep [s]
% z = 5-by-nt matrix, stores layer thickness over the nt timesteps [m]
% L = 1-by-nt vector, stores lithosphere thickness over nt timesteps [m]
% dz = 2-by-nt matrix, stores grid spacing in lithosphere, dz(1,n) and
%      sub-lithospheric mantle, dz(2,n) [m]
% H = nz-by-nt matrix, stores the heat production rate at each gridpoints
%     over the nt timesteps [W m-1K-1]
% density = nz-by-nt matrix, stores density at each gridpoints over the nt
%           timesteps [kg m-3]
%
%==========================================================================
year = 365.25*3600*24;                % seconds per year [s]
Ma = 1e6*year;                        % seconds per Ma [s]
dat = xlsread('data','B9:C18');

bb = find(isnan(dat(:,1)) == 0);

data = zeros(numel(bb),2);

for i = 1:numel(bb)
    data(i,:) = dat(bb(i),:);
end

n = numel(data(:,1));              % evaluate the number of user-defined 
                                   % checkpoints
nz = nzlit+nzsublit;

densityvalues = [2500 2700 2900 3500 3500];

%-----CHECK FOR CONSISTENCY OF DATA-SET------------------------------------
                                 
check0 = zeros(numel(data(:,1)),numel(data(1,:)));

for i = 1:numel(data(:,1))
    for j = 1:numel(data(1,:))
        check0(i,j) = isnan(data(i,j)); 
    end
end

% locate non-numeric values in the user-defined dataset

[row,col] = find(check0);         % find the index of the column that  
                                  % returns error
                         
    if 1 == isempty(row) && isempty(col) 
        
        if n < 2
            msgbox({'At least 2 checkpoints required!';...
                'Please, modify crust_time.xls file';...
                'and save changes befor going on with ';...
                'GEOTHERM'})
            
            data = xlsread('data','B9:C18');
        end
        
        pause
        
    else msgbox({'data are not consistent! ERROR: ';...
            'some value is missing, or written';...
            'in invalid numeric format. Please';...
            ', check the dataset consistency  ';...
            'and save changes before going    ';...
            'on with GEOTHERM.m               '})
        
        data = xlsread('data','B9:C18');
    end
    
%-----EVALUATE USER-DEFINED CONSTRAINTS------------------------------------
%data = sort(data,'descend');

t0 = data(1,1);           % define the initial time, start of experiment
tfin = data(n,1);         % define the final time, end of experiment  
rates = zeros(1,n-1);     % initialize the displacement/rate vector
dtn = zeros(1,n-1);       % initialize the number of timesteps within each 
                          % time interval
% ensures the number of dt is equal to predefined nt
dtn(n-1) = nt-sum(dtn(1:n-2));
                         
dt = ((t0-tfin)*Ma)/nt;   % define the time step
dh = zeros(1,nt);         % initialize the finite displacement vector 
                          % [m/timestep]
dh(1) = 0;                

% user-defined constraints loop--------------------------------------------

for i = 1:n-1
    rates(i) = ((data(i,2)-data(i+1,2))*1e3)/((data(i,1)-data(i+1,1))*Ma);
    dtn(i) = round(((data(i,1)-data(i+1,1))*Ma)/dt);
    dh(sum(dtn(1:i))-dtn(i)+1:sum(dtn(1:i))) = dt*rates(i)*ones(1,dtn(i));    
end

%-----UPDATE LAYER THICKNESS AND HEAT PRODUCTION RATES VECTOR--------------
%
variation = cumsum(dh);                % bulk lithosphere thickness change 
L = zeros(1,nt);                       % lithosphere thickness
Moho = zeros(1,nt);                    % Moho depth
Moho(1) = sum(z0(1:3));                % Moho depth at t=t0

lc = Moho(1)-z0(3);                    % lower crust depth at t=t0
L(1) = sum(z0(1:4));                   % lithosphere thickness at time = t0 

% spacing of numerical grid------------------------------------------------

dz = zeros(2,nt);                      % initialize matrix
dz(1,1) = L(1)/(nzlit-1);              % grid spacing, lit. at time = t0
dz(2,1) = (1650*1e3-L(1))/nzsublit;    % grid spacing, sublit. at time = t0

% layer thickness----------------------------------------------------------

Z = zeros(5,nt);                       % layer thicknes
variation_bulk = repmat(variation,4,1);
z_bulk = repmat(cumsum(z0(1:4)),1,nt);
Z(1:4,:) = z_bulk-variation_bulk;
Z(:,1) = z0;                           % layer thickness, time = t0

% number of gridpoints in each layer---------------------------------------

dzn = zeros(5,nt);                     % initialize matrix
dzn(1:4,1) = fix(Z(1:4,1)/dz(1,1));    % gridpoints in lit., time = t0
dzn(5,1) = nz-(sum(dzn(1:4,1)));       % gridpoints in lit., time = t0
%--------------------------------------
H = zeros(nz,nt);
Hsd = zeros(nz,nt);
Hsd_minus = zeros(nz,nt);
density = zeros(nz,nt);

for n = 2:nt
    L(n) = L(1)-variation(n);          % update lithosphere thickness [m]
    Moho(n) = Moho(1)-variation(n);    % constant sub-continental mantle 
    dz(1,n) = L(n)/nzlit-1;            % update grid spacing
    dz(2,n) = (1650*1e3-L(n))/nzsublit;
    % loop trough the z direction------------------------------------------
    for i = 1:4
        if Z(i,n) < 0
            Z(i,n) = 0;
        else Z(i+1,n) = z0(i+1);
        end
    end
    Z(5,n) = 1650*1e3-L(n);
    %----------------------------------------------------------------------
    for m = 1:5
        dzn(m,n) = fix(Z(m,n)/dz(1,n));              % number of gridpoints    
        dzn(5,n) = nz-(sum(dzn(1:4,n)));    
        % heat production rate over time
        H(sum(dzn(1:m,n))-dzn(m,n)+1:sum(dzn(1:m,n)),n) = ...
            heat_rate(m)*ones(1,dzn(m,n)); 
        H(sum(dzn(1:m,1))-dzn(m,1)+1:sum(dzn(1:m,1)),1) = ...
            heat_rate(m)*ones(1,dzn(m,1));
        Hsd(sum(dzn(1:m,n))-dzn(m,n)+1:sum(dzn(1:m,n)),n) = ...
            heat_sd(m)*ones(1,dzn(m,n));
        Hsd(sum(dzn(1:m,1))-dzn(m,1)+1:sum(dzn(1:m,1)),1) = ...
            heat_sd(m)*ones(1,dzn(m,1));
        Hsd_minus(sum(dzn(1:m,n))-dzn(m,n)+1:sum(dzn(1:m,n)),n) = ...
            heat_sd_minus(m)*ones(1,dzn(m,n));
        Hsd_minus(sum(dzn(1:m,1))-dzn(m,1)+1:sum(dzn(1:m,1)),1) = ...
            heat_sd_minus(m)*ones(1,dzn(m,1));
        density(sum(dzn(1:m,n))-dzn(m,n)+1:sum(dzn(1:m,n)),n) = ...
            densityvalues(m)*ones(1,dzn(m,n));
        density(sum(dzn(1:m,1))-dzn(m,1)+1:sum(dzn(1:m,1)),1) = ...
            densityvalues(m)*ones(1,dzn(m,1));
    end                                      
end
nzcrust = round((Moho./dz(1,n))+1);
end
