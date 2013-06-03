function [grt_S_z,grt_S_Temp,grt_L_z,grt_L_Temp,per_z,per_Temp,and_z,and_Temp,ky_z,ky_Temp,sill_z,sill_Temp] = petrology(con)
% the function allow the user to choose and plot relevant petrologic curves
% to be displayed in the main figure
%==========================================================================

% load curve limits--------------------------------------------------------

grt_S_z = con(:,2)./1e3;
grt_S_Temp = con(:,1);
grt_L_z = con(:,4)./1e3;
grt_L_Temp = con(:,3);
per_z = con(:,6)./1e3;
per_Temp = con(:,5);
and_z = con(:,8)./1e3;
and_Temp = con(:,7);
ky_z = con(:,10)./1e3;
ky_Temp = con(:,9);
sill_z = con(:,12)./1e3;
sill_Temp = con(:,11);

end

