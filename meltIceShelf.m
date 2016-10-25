%% matlab function to compute the volumetric height change of a melting floating ice shelf
%% input:
%% T0                              % in situ deg Celsius of the seawater before melting
%% S0                              % Salinity PSU of the seawater before melting
%%  meltArea                    % m^2 total area where the shelf is diluted 
%%  meltDepth                    % m depth of affected water column
%%  lon,lat                          %approximate position of where the ice is going to be melted
%%  Icethickness                   % m thickness of the ice shelf
%%  Icearea                      % m^2 %area of the to be melted shelf
%% output:
%% deltaH                        %m volumetric height difference due to melting
%% Tmean                         % deg C mean temperature of the sea waterafter melting
%% Smean                         % PSU mean salinity of the sea waterafter melting
%% W_ihmean                      % remaining mass fraction of ice versus water ( if not zero this measn not all ice has melted, and the sea water reached freezing temperature)\
%%
%% Author R. Rietbroek, Oct 2016
function  [deltaH,Tmean,Smean,W_ihmean] = meltIceShelf(T0,S0,meltArea,meltDepth,lon,lat, Icethickness,Icearea,Tice)


%some constants used here
rho_w=1027;                         % kg/m^3 average density of sea water
rho_ice=917;                        % kg/m^3 average density of ice
g=9.81;                             % m/s^2 mean gravity 
h2dbar=g*rho_w*1e-4;                %db/m conversion factor from meter of sea water to decibar
ndp=50;                             %discretize depth steps in ndp 

%% set up sea pressure column (since mass is conserved this remains constant throughout the melting process)
p=[0:(meltDepth*h2dbar)/ndp:meltDepth*h2dbar];

SA0=gsw_SA_from_SP(S0,meltDepth/2,lon,lat); % g/kg initial absolut salinity of the sea water computed from 34.5 PSU found in literature (use half the water depth to the representative depth)
%% set salinity to be constant for all depths
SAin=SA0*ones(1,length(p));

CTin=gsw_CT_from_t(SAin,T0*ones(1,length(p)),p);



%%set up sea-ice mixture by making a homogenous slush puppy (constant ice/water mass ratio over depth)
icemass=Icearea*Icethickness*rho_ice; %kg of ice
totalmass=meltArea*meltDepth*rho_w; %kg/m2 (the mass of the ice follows from Archimedes)

w_Ih_in=icemass/(totalmass)*ones(1,length(p));
t_Ih_in=Tice*ones(1,length(p));

[SA_final, CT_final, w_Ih_final] = gsw_melting_ice_into_seawater(SAin,CTin,p,w_Ih_in,t_Ih_in);


%compute the difference in volumetric height change before and after melting by integartion of the specific volume

spvolin=gsw_specvol(SAin,CTin,p);
spvolfinal=gsw_specvol(SA_final,CT_final,p);

%integrate density over depth with the trapezium rule (rho_w serves as a reference value for the sea water density)

sthin=trapz(p/h2dbar,spvolin*rho_w);
sthfinal=trapz(p/h2dbar,spvolfinal*rho_w);

%change in volumetric height
deltaH=sthfinal-sthin;


%compute in situ temperature
t = gsw_t_from_CT(SA_final,CT_final,p);

%final average temperature
Tmean=mean(t);

%final mean ice water ratio (should be zero for the entire sheet to be melted)
W_ihmean=mean(w_Ih_final);

if W_ihmean > 0
fprintf(1,'Note: not all ice has melted, with this configuration\n'); 
end

%final mean salinity (PSU)

[SP, in_ocean] = gsw_SP_from_SA(SA_final,p,lon,lat);
Smean=mean(SP);



