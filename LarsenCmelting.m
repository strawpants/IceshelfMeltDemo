%% A back-on-the envelope calcuation on computing local sea level change when a chunk of Larsen C would melt
%% Uses the Gibb's sea water toolbox version 3.05
%% Author: R. Rietbroek, 25 Oct 2016
clear all
close all

%% set up initial parameters
rho_w=1027;                         % kg/m^3 average density of sea water
rho_ice=917;                        % kg/m^3 average density of ice
g=9.81;                             % m/s^2 mean gravity 
h2dbar=g*rho_w*1e-4;                %db/m conversion factor from meter of sea water to decibar


%iceshelf and melting basin characteristics
LarsenCarea=6000e6;
meltarea=3*LarsenCarea;
Hoce=2000;

%approximate lon and latitude of point where the ice sheet is melted
lon=-60;
lat=-68;

%% intial salinity of sea water
S0=34.5;

Hice=200;                           % m ice shelf thickness
Tice=-10;                           % Deg Celsius insitu Temperature of the ice





%make a temperature dependent plot
Tin=[-2:0.1:2 2:22];
dH=[];
Tout=[];
Sout=[];
W_ihout=[];
for i=1:length(Tin)
    [dH(i),Tout(i),Sout(i),W_ihout(i)] = meltIceShelf(Tin(i),S0,meltarea,Hoce,lon,lat,Hice,LarsenCarea,Tice);
    
end


%make some plots

fw=16;
lw=3;
figure(1)
clf
[hAx,hLine1,hLine2]=plotyy(Tin,dH*1e2,Tin,Tout-Tin);
%  [hAx,hLine1,hLine2]=plotyy([Tin(ivecmelted),Tin(ivecpartialmelt)],[dH(ivecmelted)*1e2,dH(ivecpartialmelt)*1e2],Tin,Tout-Tin);

%create a background path to delineate the partial metling section
idx=find(eq(W_ihout,0));
xmelt=Tin(idx(1));
x1=xlim;
y1=ylim;

xp=[x1(1) x1(1)  xmelt xmelt];
yp=[y1(1) y1(2) y1(2) y1(1)];
col=[99 184 255]/255;
a=patch(xp,yp,col);
uistack(a,'down')
uistack(a,'down')

%some additonal info
text((xmelt+x1(1))/2,y1(2)*0.9,sprintf('Partial\nmelting'),'FontSize',fw,'HorizontalAlignment','center');
text((xmelt+x1(2))/2,y1(2)*0.9,sprintf('Complete\nmelting'),'FontSize',fw,'HorizontalAlignment','center');

grid on
set(hLine1,'LineWidth',lw,'LineStyle','--');
set(hLine2,'LineWidth',lw);
 
xlabel('Initial sea water temperature [^{\circ}C]','FontSize',fw)
ylabel(hAx(1),'Sea surface height change [cm]','FontSize',fw) % left y-axis
ylabel(hAx(2),'Sea water Temp. change [^{\circ}C]','FontSize',fw) % right y-axis
set(hAx(1),'FontSize',fw)
set(hAx(2),'FontSize',fw)
print -dpng LarsenCSealevel.png

%  figure(2)
%  clf
%  plot(Tin,100*(S0-Sout)/S0,'b-','LineWidth',3)
%  xlabel('Initial temperature [^{\circ}C]')
%  ylabel('Salinity decrease [ %]')
%  hold on
%  
%  figure(3)
%  clf
%  plot(Tin,W_ihout*100,'b-','LineWidth',3)
%  xlabel('Initial temperature [^{\circ}C]')
%  ylabel('remaining ice fraction [%]')
%  hold on


