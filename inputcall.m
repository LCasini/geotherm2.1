function [z0,heat_rate,heat_sd, heat_sd_minus] = inputcall
% the function prompt the user to select appropriate inputs to the model.
% Note that layer thickness must be specified in [m], and volumetric heat
% production rates must be set in [uW/m3].
%--------------------------------------------------------------------------
% outputs are:
% z0 = 5-by-1 vector, stores the thickness of layers
% L0 = scalar, stores the bulk thickness of lithosphere
% dzn0 = 1-by-5 vector, stores the number of gridpoints in each layer
% dz0 = 1-by-2 vector, stores the spacing of gridpoints in lithosphere 
% dz0(1,1), and sub-lithospheric mantle dz0(1,2)
% heat_rate = 5-by-1 vector,stores the specific heat production rates of 
% layers
%==========================================================================
data = xlsread('data','G10:J13');        % load user-defined data

% define model geometry at time = t0

z0 = 1e3*[data(:,1);1650-sum(data(:,1))];% layer thickness [m]                                            
heat_rate = [data(:,2);data(4,2)];       % specific heat production rate      
                                         % of the 5 layers [W/m3]
heat_sd = [data(:,3);data(4,3)];         % standard deviation to heat-
                                         % production rate, plus
heat_sd_minus = [data(:,4);data(4,4)];   % standard deviation to heat-
                                         % production rate, minus                                       
% check the dataset--------------------------------------------------------

check0 = xlsread('data','C9:C9');

if check0 == sum(z0(1:3))/1e3;
else msgbox({'WARNING: layers are not consistent';...
        'check bulk initial crustal thickness'});
end
end
