function change = lowercrust(melt,nzcrust,nzuppercrust,nz)
% the function calculates the change of lower crustal thickness as number
% of grid points
%==========================================================================
 
partial_change = zeros(1,nz);  

for j = nzuppercrust+1:nzcrust
    % evaluate threshold for melt collection (0.04)
    if melt(j) >= 0.04 
        partial_change(j) = melt(j) - 0.02;
    end
    % evaluate threshold for refractory lower crust
    if partial_change(j) > 0.35
        partial_change(j) = 0.35;
    end
end
                                                     
change = round(sum(partial_change));                                                    
                                                     
end