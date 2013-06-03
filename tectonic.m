function [startime,depth,width,Ttect] = tectonic(t0,dt,nzlit,nzcrust,dz)
% the function tectonic.m allow the user to simulate different tectonic
% settings, including: i) impinging of mantle plume beneath the
% lithosphere, ii) breakoff of subcontinental mantle
%==========================================================================
data = xlsread('data','N10:N11');
year = 365.25*3600*24;                               % seconds per year [s]
Ma = 1e6*year;                                       % seconds per Ma [s]



breakoff = 1+nzcrust;

% evaluate the specified tectonic setting----------------------------------

if sum(data) > 0
    
    if (data(1,1)) > 0
        
        T = 1530+273.15; 
        startime = round((t0-data(1,1))*Ma/dt);
        d = nzlit+1;
        w = d+round(80000/dz(2,startime)); 
        
    elseif (data(2,1)) > 0
        
        T = 1350+273.15;
        startime = round((t0-data(2,1))/(dt/Ma));
        d = breakoff(startime);
        w = d+round(12000/dz(1,startime));
    
    end
    
    Ttect = T;
    depth = d;
    width = w;
else
    depth = 1; width = 1; Ttect = 273.15; startime = t0;
end
end













