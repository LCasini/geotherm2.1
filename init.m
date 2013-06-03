function T = init(Tsurf,Tmantle,H,nz,nzcrust,dz,nzlit)
% the function setup the initial temperature profile based on constant
% thermal conductivity, density and basal heat flow at 1650 km [Wm-2] 
%--------------------------------------------------------------------------
T = zeros(1,nz);                             % preallocate T vector
h0 = H'/1e3;                                 % heat plroduction rate [W/m3]

h0(1:nzlit) = h0(1:nzlit)*dz(1,1);
h0(nzlit+1:nz) = h0(nzlit+1:nz)*dz(2,1);

q0 = (Tmantle/165e4)*5500;                   % basal heat flow, 1650 km 
kappa = zeros(1,nz);                         % static thermal conductivity
   kappa(1:nzcrust) = 3.5;                   % average crust
   kappa(nzcrust+1:nz) = 4:(nz-nzcrust):5.5; % average mantle
        
q = wrev(cumsum(wrev(h0))+q0);               % incremental heat flow

% evaluate initial temperature profile

for i = 1:nz
    T(i) = (q(i)/kappa(i))-h0(i)/(2*kappa(i));
end

T = cumsum(T)+Tsurf;

end
