function [melt,pers,perl] = M(T,z,Moho,lc)
% the function M calculates the proportion of melt/solid 
% according to (Burg & Gerya, 2005):
% M = 0 for T < Tsolidus,
% M = (T-Tsolidus)/(Tliquidus-Tsolidus) for Tsolidus<T<Tliquidus,
% M = 1 for T> Tliquidus
%==========================================================================

gs =  1398*(z.^-0.04);                          % granite solidus
gs(1) = 1080;
gl = 1767*(z.^-0.05);                           % granite liquidus  
gl(1) = 1187;

tons = zeros(1,numel(z));

for j = 1:numel(z)
    % tonalite solidus
    if z(j) <=49000
        tons(j) = -3e-8*(z(j).^2)+0.004*z(j)+1076;
    elseif z(j) > 49000 && z(j) <= 74375
        tons(j) = -4.17e-4*z(j)+1220;
    elseif z(j) > 74375 && z(j) <= 94500
        tons(j) = 0.005*z(j)+765.9;
    else
        tons(j) = -0.006*z(j)+1883;
    end
    if tons(j) < 273.15
        tons(j) = 273.15;
    end
end

tonl = tons+120; % rough tonalite liquidus

pers = zeros(1,numel(z));
perl = zeros(1,numel(z));

for j  = 1:numel(z)
    % peridotite solidus/liquidus
    if z(j) < 525000
        pers(j) = -4.7e-009*(z(j).^2)+0.004*z(j)+1400; 
        perl(j) = -9e-010*(z(j).^2)+0.0015*z(j)+2100;  
    else 
        pers(j) = 3.7e-15*(z(j).^3)-9.7e-9*(z(j).^2)+0.0087*z(j)-110;
        perl(j) = 2.6e-15*(z(j).^3)-6.9e-9*(z(j).^2)+0.0064*z(j)+680;
    end
end

% the function calculates the melt proportion within the vector z

melt = zeros(1,numel(z)); % initialize the melt_ratio vector

for i = 1:numel(z)
    if z(i) <= lc                         % loop through upper+middle crust
        if T(i) <= gs(i) 
            melt(i) = 0;
        elseif T(i) > gs(i) && T(i) <= gl(i)
            melt(i) = (T(i)-gs(i))./(gl(i)-gs(i));
        else melt(i) = 1;
        end
    elseif z(i) > lc && z(i) <= Moho       % loop through lower crust
        if T(i) <= tons(i) 
            melt(i) = 0;
        elseif T(i) > tons(i) && T(i) <= tonl(i)
            melt(i) = (T(i)-tons(i))./(tonl(i)-tons(i));
        else melt(i) = 1;
        end
    else
        if T(i) < pers(i)
            melt(i) = 0;
        elseif T(i) > pers(i) && T(i) <= perl(i)
            melt(i) = (T(i)-pers(i))./(perl(i)-pers(i));
        else melt(i) = 1;
        end
    end
end

end
 