function bestfit = compare(thermo_Tav,thermo_Pav,Tsurf,Tmantle,T,z,nz,dz)
% the function compare is used to evaluate the degree of fitting between
% natural data provided by thermobarometric/age constraints and the
% calculated geotherm at the current timestep
%==========================================================================

check = numel(thermo_Tav)-sum(isnan(thermo_Tav/0));

% compare natural data (if available) and numerical results

if check > 0
    
    % initialize vectors--------------------------------------------------- 
    
    Tdata = zeros(1,check+2);
    Tdata(1) = Tsurf;                          % enforce BC
    Tdata(check+2) = Tmantle;                  % enforce BC
    Pdata = zeros(1,check+2);
    Pdata(1) = 0;                              % enforce BC
    Pdata(check+2) = 1650*1e3;                 % enforce BC
    
    % find P,T values of constraints at the current timestep---------------
    
    [~,~,T_values] = find(thermo_Tav);         % temperature [K]
    [~,~,P_values] = find(thermo_Pav*1e3);     % pressure [m]
    
    Pmax = round(max(P_values)/dz);            % find the lower and upper  
    Pmin = round(min(P_values)/dz);            % points of the evaluation
                                               % interval
    Tdata(2:check+1) = T_values;
    Pdata(2:check+1) = P_values;
    
    % interpolate natural data---------------------------------------------
    
    [p,~,mu] = polyfit(Pdata,Tdata,2);         % find the coefficients of 
                                               % second-order polynomial 
    Tnatural = zeros(nz,1);

    for i = 1:nz
        Tnatural(i) = p(1).*((z(i)-mu(1))/mu(2)).^2 +...
            p(2)*((z(i)-mu(1))/mu(2)) + p(3);
    end
    
    % find the relative error
    
    bestfit = median(100*(Tnatural(Pmin:Pmax)-T(Pmin:Pmax))./T(Pmin:Pmax));
    
else                                         
    bestfit = [];
end

end