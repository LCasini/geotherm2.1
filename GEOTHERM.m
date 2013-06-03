% GEOTHERM
% Leonardo Casini 2013, University of Sassari, Italy.
%==========================================================================
% GEOTHERM calculates the shape of transient 1D geotherms in a deforming
% crustal section. The model crust is composed of three layers-crust and 
% two-layers mantle. Each layer is characterized by specific heat 
% production rate, density, and thickness. Dirichlet boundary conditions  
% are used at both ends of the profile (0 - 1650 km). 
% 
% The upper node is keep at constant 273.15 K temperature (surface),
% whereas lower node is keep at 2323.15 K - this gives approximately 1300°C 
% at 150 km depth in a stable crust, using a linear 0.5°C/km mantle adiabat
%==========================================================================
clear all
close all
msgbox({'Fill required fields in tables 1,2,3 and 4 of data.xls ';... 
    'before going on with simulation using GEOTHERM.m. ';...
    'Please, follow instructions in the spreadsheet'});
                                      % warning message
pause

%-------SET PHYSICAL PARAMETERS (default values)---------------------------
Tsurf = 273.15;                       % surface temperature [K]
Tmantle = 2050+273.15;                % temperature at 1650 km depth [K]
Tplume = 1500+273.15;                 % temperature of mantle plume [K]
year = 365.25*3600*24;                % seconds per year [s]
Ma = 1e6*year;                        % seconds per Ma [s]

%-------SET NUMERICAL PARAMETERS-------------------------------------------
nzlit = 101;                          % gridpoints, lithosphere
nzsublit = 50;                        % gridpoints, sub-lithosphere
nz = nzlit+nzsublit;                  % number of gridpoints in z-direction
                                      % only the first 101 points
                                      % (lithosphere above LAB) are
                                      % displayed in the plot
nt = 101;                             % number of time steps to compute

Newnt = input('Set number of time steps (default is 100): ');
if Newnt > 0                          % the user has set the number of time 
                                      % step to compute
    nt = Newnt+1;
else                                  % the user hint 'enter' 
end

pause(0.5);

%-------USER-DEFINED PARAMETERS--------------------------------------------

[z0,heat_rate,heat_sd,heat_sd_minus] = inputcall;    % load user-defined 
                                                     % inputs from data.xls 
                                                     % file (initial
                                                     % configuration)

[t0,tfin,L,dz,H,Hsd,Hsd_minus,Moho,lc,nzcrust,...    % evaluate geological 
    density,variation] = nature(nt,z0,nzlit,...      % history
    nzsublit,heat_rate,heat_sd,heat_sd_minus);     

nzuppercrust = round((lc./dz(1,1))+1);               % grid number at the
                                                     % interface middle/
                                                     % lower crust

dt = ((t0-tfin)/nt)*Ma;                              % calc. time step [s]

%-------setup geodynamic setting and deformation---------------------------
 
[startime,depth,width,Ttect] = tectonic(t0,dt,nzlit,nzcrust,dz);

granite = 2.37;                                      % model heat 
                                                     % production rate for
                                                     % granite
                                                     
%-------read and evaluate petrologic curves to be displayed----------------

con = xlsread('petrology','A5:L15');

[grt_S_z,grt_S_Temp,grt_L_z,grt_L_Temp,per_z,per_Temp,...
    and_z,and_Temp,ky_z,ky_Temp,sill_z,sill_Temp] = petrology(con);  

%-------read thermobarometric constraints to be displayed------------------

[thermo_Tlow,thermo_Thigh,thermo_Plow,thermo_Phigh,thermo_n] =...
    thermobarometry(Ma,nt,dt,t0,tfin);

thermo_Tav = thermo_Tlow+(thermo_Thigh-thermo_Tlow)/2; 

thermo_Pav = thermo_Plow+(thermo_Phigh-thermo_Plow)/2;

%==========================================================================
% setup initial conditions

time = t0;                                         % setup initial time

z = [linspace(-1650*1e3,-L(1),nzsublit),...        % setup depth vector
    linspace(-L(1)+dz(1),0,nzlit)];
z = wrev(abs(z));

% Initial and Dirichlet boundary conditions
T = init(Tsurf,Tmantle,H(:,1),nz,nzcrust(1),dz(:,1),nzlit);
Tminus = init(Tsurf,Tmantle,H(:,1)+Hsd(:,1),nz,nzcrust(1),dz(:,1),nzlit);
T(1) = Tsurf;                                      % upper BC
T(nz) = Tmantle;                                   % lower BC

T = T';                                            % make T a column vector

%-------initialize movie structure-----------------------------------------

nFrames = nt;                          % initialize the number of movie 
                                       % frames
mov(1:nFrames) = struct('cdata',...    % preallocate movie structure
    [],'colormap', []);                                      

%-------CALCULATES NEW TEMPERATURES----------------------------------------

timestep = 0;                                % setup predefined numbering
v_bestfit  = NaN*ones(1,nt);
t_bestfit = NaN*ones(1,nt);
v_Sq = NaN*ones(1,nt);t_Sq = NaN*ones(1,nt);
t_Mq = NaN*ones(1,nt);v_Mq = NaN*ones(1,nt);

volume_change = 0;
melt_ext = 0;                                % initialize melt extracted,
                                             % lower crust
for n = 1:nt
    
    % time step loop
    timestep = timestep+1;
    z = [linspace(-1650*1e3,-L(n),nzsublit),linspace(-L(n)+dz(1),0,nzlit)];
    z = wrev(abs(z));                        % reshape the z-vector
    
    % update melt proportion, viscosity, stress and shear heating
    
    [melt,pers,perl] = M(T,z,Moho(n),lc);    % update melt proportion
    
    viscosity = ni(melt,nz,nzlit,...         % update dynamic viscosity
        nzcrust(n));
    
    [shear_heating,stress] = VSH(...         % update shear heating
        nz,t0,dt,nzlit,T,timestep...
        ,nzcrust(n),melt,z,viscosity);
    
    %======================================================================
    % update lower crust thickness/composition if partial melting occur
        
    change = lowercrust(melt,nzcrust(n),nzuppercrust,nz); 
    volume_change = volume_change+change*dz(1,n);
   
    if volume_change > 0.35*z0(3)
        volume_change = 0.35*z0(3);
    end 
    
    lc = (Moho(1)-z0(3)-variation(n))+volume_change;
    nzuppercrust = round((lc./dz(1,n))+1);
    
    time = time-(dt/Ma);                      % update time [Ma]
    [Cp,k] = kappa(T,nz);                     % update diffusivity and heat
                                              % capacity
    %======================================================================
    
    % calculates new temperatures
    
    number = 1:nz;                            % Setup numbering       
    s = zeros(1,nz);                          % initialize constant
    b = zeros(1,nz);
    f = zeros(1,nz);
    b(2:nz-1) = k(1:nz-2)+k(2:nz-1);          % backward approximation 
    f(2:nz-1) = k(2:nz-1)+k(3:nz);            % forward approximation
    
    s(1:nzlit) = dt/(2*dz(1,n)^2);            % constant 
    s(nzlit+1:nz) = dt/(2*dz(2,n)^2);         % constant
    
    cost = 1 + s.*b + s.*f;
    
    % construct the A matrix at the current timestep
    
    A = zeros(nz);        
    A(1,1) = 1; A(nz,nz) = 1; 
    
    for j = 2:nz-1
            i = number(j);
            A(i,j-1) = -s(j)*f(j);            
            A(i,j) = cost(j);               
            A(i,j+1) = -s(j)*b(j);    
    end
    
    rhs = T + dt*((H(:,n).*1e-6 + shear_heating')./(density(:,n).*Cp'));
    rhs_minus = T + dt*(((H(:,n)-Hsd_minus(:,n)).*1e-6 +...
        shear_heating')./(density(:,n).*Cp'));
    rhs_plus = T + dt*(((H(:,n)+Hsd(:,n)).*1e-6 +...
        shear_heating')./(density(:,n).*Cp'));
    Tnew_minus = A\rhs_minus;
    Tnew_plus = A\rhs_plus;
    Tnew_minus(1) = Tsurf; Tnew_plus(1) = Tsurf;
    Tnew_minus(nz) = Tmantle; Tnew_plus(nz) = Tmantle;
    
    % solve for new temperatures
    
    Tnew = A\rhs;                             % calculate temperature at 
                                              % the current timestep
    Tnew(1) = Tsurf;                          % enforce boundary condition
    Tnew(nz) = Tmantle;                       % lower BC
    
    % update tectonic setting----------------------------------------------
    
    if timestep == startime                 
        Tnew(depth:width) = Ttect;
    else
    end
    
    if volume_change < 0.35*z0(3)
        
        emplacement = find(Tnew >= 773.15);   % find the nz point which  
                                              % is the roof of intrusions
        Tnew(emplacement(1):emplacement(1)+change-1) = 680+273.15;
        H(emplacement(1):emplacement(1)...
            +change-1,n) = granite; 
        H(nzuppercrust+1:nzuppercrust+1+...
            change-1,n) = H(nzuppercrust,n);  
        
        density(emplacement(1):emplacement(1)...
            +change-1,n) = 2500;
    else
        melt(nzuppercrust+1:nzcrust(n)) = 0.02;
    end
    
    T  = Tnew;                                % update temperature [K] 
    Tminus = Tnew_minus;
    Tplus = Tnew_plus;
    
    % update heat flow-----------------------------------------------------
    
    Sq = sum(H(1:nzcrust(n)).*dz(1,n)./1e3)+...
        (T(nzcrust(n)+1)*1e3)./z(nzcrust(n)+1);
    
    Mq = (T(nzcrust(n)+1)*1e3)./z(nzcrust(n)+1); 
       
    
    % compare experimental and numerical data==============================   

    bestfit = compare(thermo_Tav(:,n),thermo_Pav(:,n),...
        Tsurf,Tmantle,T,z,nz,dz(1,n));
  
    % updat geometry-------------------------------------------------------
    
    LAB = find(T >= pers');                   % update astenoshpere 
    lithosphere_Temp = [273.16 T(nz)];        % update lithosphere T
    
    if isempty(LAB) == 0
        lithosphere_z = [z(LAB(1))...         % update lithosphere z
            z(LAB(1))]./1e3; 
    else
        lithosphere_z = [250 250];
    end
    
    Moho_z = [Moho(n)/1e3 Moho(n)/1e3];       % update Moho thickness
    Moho_Temp = [273.16 T(nz)];               % update Moho T
    lc_z = [lc/1e3 lc/1e3];                   % update lower crut thickness
    lc_Temp = [273.16 T(nz)];                 % update lower crust T
                                   
%========PLOT THE SOLUTIONS================================================

if (mod(n,1) == 0)                            % plot solution every 2
                                              % timestep
    figure(1), clf
    subplot(10,13,[15:19 28:32 41:45 54:58 67:71 80:84 93:97 106:110]);
    hold on
    
    % add Moho, lower crust,LAB and shade the section----------------------
    
    area(Moho_Temp,[1650*1e3 1650*1e3],'FaceColor',...
        [1 0.74 0.23],'EdgeColor',[1 0.74 0.23]);
    
    plot(lithosphere_Temp,lithosphere_z,'--w',...
        'LineWidth',3);
    
    area(Moho_Temp,Moho_z,'FaceColor',...
        [0.61 0.51 0.89],'EdgeColor',[0.61 0.51 0.89]);
    
    plot(Moho_Temp,Moho_z,'--w',...
        'LineWidth',2);
    
    text(300,Moho_z(1),'Moho','FontSize',8,'FontWeight','bold',...
        'Fontname','Helvetica','BackgroundColor','w');
    
    text(300,lithosphere_z(1),'LAB','FontSize',9,'FontWeight',...
        'bold','Fontname','Helvetica','BackgroundColor','w');
    
    area(lc_Temp,lc_z,'FaceColor',...
        [0.45 0.30 0.85],'EdgeColor',[0.45 0.30 0.85]);
    
    plot(lc_Temp,lc_z,'--w');
    
    % add petrologic boundaries-------------------------------------------- 
    
    plot(grt_S_Temp,grt_S_z,'k','LineWidth',1);
    plot(grt_L_Temp,grt_L_z,':k','LineWidth',1);
    plot(per_Temp,per_z,'Color',[0 0 0.62],'LineWidth',2);
    plot(and_Temp,and_z,'k','LineWidth',1);
    plot(ky_Temp,ky_z,'k','LineWidth',1);
    plot(sill_Temp,sill_z,'k','LineWidth',1);
    
    % display thermobarometric constraints at the current timestep---------
    
    set (gca, 'Xdir', 'normal','Ydir','reverse');
    set (gca,'XAxisLocation','Top');
    xlabel ('temperature [K]','FontWeight','bold','FontSize',8);
    ylabel ('depth [km]','FontWeight','bold','FontSize',8);
    set(gca,'Ylim',[0 150]);
    set(gca,'XLim',[273.15,1700]);
    title(['time is = ',num2str(time),' Ma']);
    set(gca,'YTick',[0 10 20 30 40 50 100 150],'TickDir','Out')
    plot(T,z/1e3,'w','LineWidth',2);
    
    % evaluate standard deviation
    
    sd_plus = Hsd(:,n)./H(:,n);
    sd_minus = Hsd_minus(:,n)./H(:,n);
    plot(Tplus,z/1e3,':w');
    plot(Tminus,z/1e3,':w');
    
    for j = 1:thermo_n 
        plot([thermo_Tlow(j,n),thermo_Thigh(j,n)],[thermo_Pav(j,n),...
            thermo_Pav(j,n)],'Color','k','LineWidth',1);
        plot([thermo_Tav(j,n),thermo_Tav(j,n)],[thermo_Plow(j,n),...
            thermo_Phigh(j,n)],'Color','k','LineWidth',1);
        plot(thermo_Tav(j,n),thermo_Pav(j,n),'wo','MarkerSize',3,...
            'MarkerEdgeColor','k','MarkerFaceColor','w',...
            'LineWidth',1);
    end
    
    % display melt proportion at current timestep--------------------------
    
    subplot(10,13,[21:22 34:35]);
    hold on
    plot(melt,z./1e3,'r','LineWidth',2);
    set (gca, 'Xdir', 'normal','FontSize',6,'Ydir','reverse','FontSize',6);
    set (gca,'XAxisLocation','Top');
    ylabel ('depth [km]','FontWeight','bold','FontSize',6);
    xlabel ('','FontSize',4);
    set(gca,'Ylim',[0 60]);
    set(gca,'XLim',[0,1]);
    title(['melt, ',num2str(time),' Ma'],'FontSize',6);
    plot([0 1],Moho_z,':k');
    plot([0 1],[z(nzuppercrust)./1e3 z(nzuppercrust)./1e3],'--k');
    set(gca,'YTick',[0 20 40 60],'TickDir','Out')
    
    % display deviatoric stress at current timestep------------------------
    
    subplot(10,13,[24:25 37:38]);
    hold on
    plot(stress,z./1e3,'b','LineWidth',2);
    set (gca, 'Xdir','normal','FontSize',6,'Ydir','reverse','FontSize',6);
    set (gca,'XAxisLocation','Top');
    ylabel ('depth [km]','FontWeight','bold','FontSize',6);
    xlabel ('','FontSize',4);
    set(gca,'Ylim',[0 60]);
    set(gca,'XLim',[1,500]);
    title(['stress, ',num2str(time),' Ma'],'FontSize',6);
    plot([1 500],Moho_z,':k');
    plot([1 500],[z(nzuppercrust)./1e3 z(nzuppercrust)./1e3],'--k')
    text(515,z(nzuppercrust)./1e3,'Lower Crust','FontSize',6,...
        'FontWeight','bold','Fontname','Helvetica');
    text(515,Moho_z(1),'Moho','FontSize',6,'FontWeight','bold',...
        'Fontname','Helvetica');
    set(gca,'YTick',[0 20 40 60],'TickDir','Out')
    
    % display median error in the crust the current timestep---------------
    if isempty(bestfit) == 0
        v_bestfit(n) = bestfit;
        t_bestfit(n) = time;
    else
        v_bestfit(n) = 0;
        t_bestfit(n) = time;
    end
    subplot(10,13,[60:64 73:77]);
    hold on
    plot(t_bestfit,v_bestfit,'b','LineWidth',2);
    plot([t0 tfin],[0 0],':k');
    text(tfin-1,0,'best fit','FontSize',6,'FontWeight','bold',...
        'Fontname','Helvetica');
    text(tfin-1,-25,'T underestimated','FontSize',6,'FontWeight','bold',...
        'Fontname','Helvetica');
    text(tfin-1,+25,'T overestimated','FontSize',6,'FontWeight','bold',...
        'Fontname','Helvetica');
    set (gca, 'Xdir','reverse','FontSize',6,'Ydir','normal','FontSize',6);
    title(['fit at = ',num2str(time),' Ma'],'FontSize',8);
    ylabel ('error [%]','FontWeight','bold','FontSize',6);
    set(gca,'Ylim',[-50 50]);
    set(gca,'XLim',[tfin,t0]);
    set(gca,'TickDir','Out');
    
    % display Moho and surface heat flow at the current timestep
    
    v_Sq(n) = Sq;
    t_Sq(n) = time;
    v_Mq(n) = Mq;
    t_Mq(n) = time;
    subplot(10,13,[99:103 112:116]);
    hold on
    plot(t_Sq,v_Sq,'r','LineWidth',2);
    plot(t_Mq,v_Mq,'b','LineWidth',2);
    set (gca, 'Xdir','reverse','FontSize',6,'Ydir','normal','FontSize',6);
    title(['heat flow at = ',num2str(time),' Ma'],'FontSize',8);
    ylabel ('heat flow [mW/m2]','FontWeight','bold','FontSize',6);
    set(gca,'Ylim',[0 200]);
    set(gca,'XLim',[tfin,t0]);
    text(298,190,'red = Surface','FontSize',6,'FontWeight','bold',...
        'Fontname','Helvetica');
    text(298,160,'blue = Moho','FontSize',6,'FontWeight','bold',...
        'Fontname','Helvetica');
    set(gca,'TickDir','Out');
    drawnow
    
    mov(n) = getframe(gcf);     % capture movie frames
end  

end

% export simulation to video (AVI format, default)
filename = input('Save results as experiment number: '); 
movie2avi(mov, ['experiment ',num2str(filename),'.avi'],'compression',...
    'none');




