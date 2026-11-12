% This script generates an example plot of inferred density, N2, M2,
% velocity, w. For Thompson et al. (Comm. Earth Env. 2026), Figure 4, use
% SWOT ID '155'.
%
% AFT, Revised version, Feb, 2026
% Last updated: May 2, 2026


%%% BEFORE RUNNING THIS SCRIPT, RUN THE CreateDensityFields.m SCRIPT FOR
%%% THE DESIRED SWOT SWATH. This will generate the density, density
%%% gradient, velocity and velocity gradient fields.


% Load the appropriate SWOT file

SWOT_swath_no = input('Enter SWOT swath no.:  ');
filename = (['SWOT_L3_LR_SSH_Expert_' SWOT_swath_no '_v1.0.2_fcs.nc'])

% Load the SSH data from a single SWOT swath
cd /Users/andy/Work/SWOT/Benoit/SWOTMaps/Expert/
lat = ncread(filename,'latitude'); 
lon = ncread(filename,'longitude');
time = ncread(filename,'time');
time = time/60/60/24 + datenum(2000,1,1,0,0,0) - datenum(2024,1,1,0,0,0);

ssha = ncread(filename,'ssha');
mdt = ncread(filename,'mdt');
ugeo = ncread(filename,'ugos');
vgeo = ncread(filename,'vgos');
ugeo_fit13 = ncread(filename,'ugos_fit_13');
vgeo_fit13 = ncread(filename,'vgos_fit_13');
ucgos_fit13 = ncread(filename,'ucgos_fit_13');
vcgos_fit13 = ncread(filename,'vcgos_fit_13');

speed = sqrt( ugeo.^2 + vgeo.^2 );
adt = ssha + mdt;

ff = sw_f(-54);

% Find the appropriate dates
date_indices = find(time(1,:)>=0 & time(1,:)<=31);
'Here are the dates available:'
time(1,date_indices)
if length(date_indices)>1
    ll = input('Choose the date entry:  ')
    choose_time = find(time(1,:)==time(1,date_indices(ll)));
    date_index = date_indices(ll);
else
    choose_time = find(time(1,:)==time(1,date_indices));
    date_index = date_indices;
end

adt_index = adt(:,:,date_index);


% Get the appropriate hydrographic data

cd /Users/andy/Work/SWOT/Benoit/GliderData/

lon_buoy = ncread(['FOCUS_data_trunc_swath_' SWOT_swath_no '_' num2str(date_index) '.nc'],'lon_buoy');
lat_buoy = ncread(['FOCUS_data_trunc_swath_' SWOT_swath_no '_' num2str(date_index) '.nc'],'lat_buoy');
buoyancy_fig = ncread(['FOCUS_data_trunc_swath_' SWOT_swath_no '_' num2str(date_index) '.nc'],'buoyancy');
w_fig = ncread(['FOCUS_data_trunc_swath_' SWOT_swath_no '_' num2str(date_index) '.nc'],'w_vel');
vort_fig = ncread(['FOCUS_data_trunc_swath_' SWOT_swath_no '_' num2str(date_index) '.nc'],'vorticity');


figure(186); close
hh = figure(186); set(hh,'Position',[ 59           1        1078         865 ])
    
choose_sec = 55;
    
pl1 = subplot(4,2,1)
pcolor(lon,lat,adt_index); shading flat; colorbar; 
hold on
hax = haxby(64); 
set(gca,'colormap',hax)
set(gca,'YTick',[-57:-53],'YTickLabel',['57';'56';'55';'54';'53'])
plot(lon(choose_sec,:),lat(choose_sec,:),'k-','LineWidth',2)
axis([151 156 -57.7 -52.5])
caxis([-1.2 0.3])
set(gca,'FontSize',16)
%xlabel('Longitude (°E)','FontSize',16)
ylabel('Latitude (°S)','FontSize',16)
t1 = text(157.65,-56.0682,'ADT (m)','FontSize',16,'Rotation',90)

pl3 = subplot(4,2,3)
pcolor(lon_buoy,lat_buoy,squeeze(buoyancy_fig(:,:,200))); shading flat; colorbar; hold on
axis([151 156 -57.7 -52.5])
plot(lon(choose_sec,:),lat(choose_sec,:),'k-','LineWidth',2)
ccc = cmocean('delta');
set(gca,'YTick',[-57:-53],'YTickLabel',['57';'56';'55';'54';'53'])
set(gca,'Colormap',ccc)
caxis([-.004 .004]);
set(gca,'FontSize',16)
%xlabel('Longitude (°E)','FontSize',16)
ylabel('Latitude (°S)','FontSize',16)
t2 = text(157.65,-56.9399,'Buoyancy (m s^-^2)','FontSize',16,'Rotation',90)

pl5 = subplot(4,2,5)
pcolor(lon_buoy,lat_buoy,cross_uc_fit13_buoy); shading flat; colorbar; hold on
axis([151 156 -57.7 -52.5])
plot(lon(choose_sec,:),lat(choose_sec,:),'k-','LineWidth',2)
ccc = cmocean('balance');
set(gca,'YTick',[-57:-53],'YTickLabel',['57';'56';'55';'54';'53'])
set(gca,'Colormap',ccc)
caxis([-1.5 1.5]);
set(gca,'FontSize',16)
%xlabel('Longitude (°E)','FontSize',16)
ylabel('Latitude (°S)','FontSize',16)
t3 = text(157.65,-57.9919,'Cross-track velocity (m s^-^1)','FontSize',16,'Rotation',90)


pl7 = subplot(4,2,7)
pcolor(lon_buoy,lat_buoy,squeeze(w_fig(:,:,200))*60*60*24); shading flat; colorbar; hold on
axis([151 156 -57.7 -52.5])
plot(lon(choose_sec,:),lat(choose_sec,:),'k-','LineWidth',2)
ccc = cmocean('balance');
set(gca,'YTick',[-57:-53],'YTickLabel',['57';'56';'55';'54';'53'])
set(gca,'Colormap',ccc)
caxis([-200 200]);
set(gca,'FontSize',16)
xlabel('Longitude (°E)','FontSize',16)
ylabel('Latitude (°S)','FontSize',16)
t4 = text(157.65,-57.8416,'Vertical velocity (m day^-^1)','FontSize',16,'Rotation',90)

pl2 = subplot(4,2,2)
kk = find(N2_sm<0); N2_sm(kk) = NaN;
pcolor(lat_buoy(choose_sec,:),depthints,squeeze(log10(N2(choose_sec,:,:)))'); shading flat; colorbar; axis ij; hold on
contour(lat_buoy(choose_sec,:),depthints,squeeze((buoyancy(choose_sec,:,:)))',[-.005:.001:.005],'Color',[1 1 1])
ccc = cmocean('solar');
set(gca,'Colormap',ccc)   
axis([-57.7 -52.7 0 950])
caxis([-6 -4]);
set(gca,'FontSize',16)
set(gca,'XTick',[-57:-53],'XTickLabel',['57';'56';'55';'54';'53'])
%xlabel('Latitude (°S)','FontSize',16)
ylabel('Depth (m)','FontSize',16)
t5 = text(-51.9,741,'log_1_0(N^2) (s^-^2)','FontSize',16,'Rotation',90)

pl4 = subplot(4,2,4)
pcolor(lat_buoy(choose_sec,:),depthints,squeeze((M2_along_sm(choose_sec,:,:)))'); shading flat; colorbar; axis ij; hold on
contour(lat_buoy(choose_sec,:),depthints,squeeze((buoyancy(choose_sec,:,:)))',[-.005:.001:.005],'Color',[.5 .5 .5])
ccc = cmocean('balance');
set(gca,'Colormap',ccc)   
axis([-57.7 -52.7 0 950])
caxis([-1.2e-7 1.2e-7]);
set(gca,'FontSize',16)
set(gca,'XTick',[-57:-53],'XTickLabel',['57';'56';'55';'54';'53'])
%xlabel('Latitude (°S)','FontSize',16)
ylabel('Depth (m)','FontSize',16)
t6 = text(-51.9,631,'M^2 (s^-^2)','FontSize',16,'Rotation',90)

pl6 = subplot(4,2,6)
pcolor(lat_buoy(choose_sec,:),depthints,squeeze((full_across_sm(choose_sec,:,:)))'); shading flat; colorbar; axis ij; hold on
contour(lat_buoy(choose_sec,:),depthints,squeeze((vort_fig(choose_sec,:,:)))'/ff,[-1:.1:-.1],'Color',[.5 .5 .5],'LineStyle','--')
contour(lat_buoy(choose_sec,:),depthints,squeeze((vort_fig(choose_sec,:,:)))'/ff,[.1:.1:1],'Color',[.5 .5 .5],'LineStyle','-')
contour(lat_buoy(choose_sec,:),depthints,squeeze((vort_fig(choose_sec,:,:)))'/ff,[0 0],'Color',[0 0 0],'LineStyle','-')
ccc = cmocean('balance');
set(gca,'Colormap',ccc)   
axis([-57.7 -52.7 0 950])
caxis([-1.5 1.5]);
set(gca,'FontSize',16)
set(gca,'XTick',[-57:-53],'XTickLabel',['57';'56';'55';'54';'53'])
%xlabel('Latitude (°S)','FontSize',16)
ylabel('Depth (m)','FontSize',16)
t7 = text(-51.9,1005,'Cross-track velocity (m s^-^1)','FontSize',16,'Rotation',90)


pl8 = subplot(4,2,8)
pcolor(lat_buoy(choose_sec,:),depth_trunc,squeeze((w_fig(choose_sec,:,:))*60*60*24)'); shading flat; colorbar; axis ij; hold on
contour(lat_buoy(choose_sec,:),depthints,squeeze((buoyancy(choose_sec,:,:)))',[-.005:.001:.005],'Color',[.5 .5 .5])
ccc = cmocean('balance');
set(gca,'Colormap',ccc)   
axis([-57.7 -52.7 0 950])
caxis([-200 200]);
set(gca,'FontSize',16)
xlabel('Latitude (°S)','FontSize',16)
set(gca,'XTick',[-57:-53],'XTickLabel',['57';'56';'55';'54';'53'])
ylabel('Depth (m)','FontSize',16)
t8 = text(-51.9,977,'Vertical velocity (m day^-^1)','FontSize',16,'Rotation',90)


set(pl1,'Position',[0.07 0.78 0.22 0.19])
set(pl3,'Position',[0.07 0.54 0.22 0.19])
set(pl5,'Position',[0.07 0.30 0.22 0.19])
set(pl7,'Position',[0.07 0.06 0.22 0.19])

set(pl2,'Position',[0.46 0.78 0.45 0.19])
set(pl4,'Position',[0.46 0.54 0.45 0.19])
set(pl6,'Position',[0.46 0.30 0.45 0.19])
set(pl8,'Position',[0.46 0.06 0.45 0.19])




