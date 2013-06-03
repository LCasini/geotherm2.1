function viscosity = ni(melt,nz,nzlit,nzcrust)
% the function ni calculates the effective viscosity of solid (melt < 0.1)
% and partially molten rocks (melt > 0.1), according to Ranalli (1995) and 
% Bittner & Schmeling (1995) respectively. 
% The input arguments 'melt',T, nz are passed to the function 
% ni.m by  the mainscript GEOTHERM
%==========================================================================

% the function calculates the effective viscosity
viscosity = zeros(1,nz); % initialize the ni vector
 
%-----loop through the crust-----------------------------------------------

% calculate viscosity

for i = 1:nz
    if i <= nzcrust 
        if melt(i) <= 0.1 
            % viscosity for solid crustal rocks [Pa s-1]                      
            viscosity(i) = 1e-21;
        else
            % viscosity of molten crustal rocks [Pa s-1]  
            viscosity(i) = (5.0*10^14)*exp(2.5+...
                (1-melt(i)).*((1-melt(i))./melt(i)).^0.48);
        end
    elseif i > nzcrust && i <= nzlit
        if melt(i) <= 0.1
            % viscosity for solid mantle rocks [Pa s-1]   
            viscosity(i) = 1e-19;
        else
            % viscosity of molten mantle rocks [Pa s-1]  
            viscosity(i) = (10^13)*exp(2.5+...
                (1-melt(i)).*((1-melt(i))./melt(i)).^0.48);
        end
    else 
        viscosity(i) = 1e-20;
    end
end
end
