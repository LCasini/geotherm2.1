function [Cp,k] = kappa(T,nz)
% the function kappa calculates the value of thermal diffusivity
% for variable temperatures, following equations 2,3 and 4,5 in 
% Whittington et al(2009)-an average rock density of 2700 kgm-3 is assumed
%==========================================================================

Cp = zeros(1,nz);                      % initialize specific heat capacity

k = zeros(1,nz);                       % initialize thermal diffusivity

% output are the values of thermal diffusivity,specific heat Cp
%--------------------------------------------------------------------------
for i = 1:nz
    if T(i) < 846
        k(i) = (567.3./T(i)-0.062)*1e-6;     % diffusivity [m2 s-1]
        Cp(i) = (199.5+0.0857.*T(i)-5.0*...  % J kg-1 K-1
            10^6.*(T(i).^-2))*4.5089;
        % togliere Cp come funzione di T e mettere il valore nominale ---------------
    else 
        k(i) = (0.732-0.000135.*T(i))*1e-6;  % diffusivity [m2 s-1]
        Cp(i) = (229.32+0.0323.*T(i)-47.9*...% J kg-1 K-1
            10^-6.*(T(i).^-2))*4.5089;
    end
end

end
