% Overview figure for Thompson et al. Nat. Comm. draft
%
% Last updated: May 2, 2026

% Add paths used for plotting
addpath /Users/andy/Work/Matlab/cmocean_v1.4/cmocean/
addpath /Users/andy/Work/Matlab/
addpath /Users/andy/Work/Matlab/colmaps/

yn = input('Do you need to load the data (y==1, n==0)?  ');

if yn==1

%%% First load some data %%%

% Load EKE
%ncdisp('means_SSH_1993_2017.nc')
loneke = ncread('means_SSH_1993_2017.nc','longitude');
lateke = ncread('means_SSH_1993_2017.nc','latitude');
eke = ncread('means_SSH_1993_2017.nc','EKE');

% Load bathymetry
%ncdisp('GEBCO2023_0.02.nc')
lonbathy = ncread('GEBCO2023_0.02.nc','lon');
latbathy = ncread('GEBCO2023_0.02.nc','lat');
bathy = ncread('GEBCO2023_0.02.nc','elevation');

% Load DUACS altimetry

cd /Users/andy/Work/SWOT/Benoit/SSHData/
sshtime = ncread('cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.25deg_P1D_1731273896002.nc','time');
sshlon = ncread('cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.25deg_P1D_1731273896002.nc','longitude');
sshlat = ncread('cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.25deg_P1D_1731273896002.nc','latitude');
sshadt = ncread('cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.25deg_P1D_1731273896002.nc','adt');
sshsla = ncread('cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.25deg_P1D_1731273896002.nc','sla');

sshtime = (sshtime/60/60/24) + datenum(1970,1,1,0,0,0) - datenum(2024,1,1,0,0,0);
cd /Users/andy/Work/Manuscripts/Andy_FOCUS/

% Load select SWOT data
% January 9
cd /Users/andy/Work/SWOT/Benoit/SWOTMaps/Expert/
lat_Jan9 = ncread('SWOT_L3_LR_SSH_Expert_155_v1.0.2_fcs.nc','latitude');
lon_Jan9 = ncread('SWOT_L3_LR_SSH_Expert_155_v1.0.2_fcs.nc','longitude');
time_Jan9 = ncread('SWOT_L3_LR_SSH_Expert_155_v1.0.2_fcs.nc','time');
time_Jan9 = time_Jan9/60/60/24 + datenum(2000,1,1,0,0,0) - datenum(2024,1,1,0,0,0);
ssha_Jan9 = ncread('SWOT_L3_LR_SSH_Expert_155_v1.0.2_fcs.nc','ssha');
mdt_Jan9 = ncread('SWOT_L3_LR_SSH_Expert_155_v1.0.2_fcs.nc','mdt');
ugeo_Jan9 = ncread('SWOT_L3_LR_SSH_Expert_155_v1.0.2_fcs.nc','ugos');
vgeo_Jan9 = ncread('SWOT_L3_LR_SSH_Expert_155_v1.0.2_fcs.nc','vgos');
speed_Jan9 = sqrt( ugeo_Jan9.^2 + vgeo_Jan9.^2 );
adt_Jan9 = ssha_Jan9 + mdt_Jan9;
cd /Users/andy/Work/Manuscripts/Andy_FOCUS/

% January 19
cd /Users/andy/Work/SWOT/Benoit/SWOTMaps/Expert/
lat_Jan19 = ncread('SWOT_L3_LR_SSH_Expert_433_v1.0.2_fcs.nc','latitude');
lon_Jan19 = ncread('SWOT_L3_LR_SSH_Expert_433_v1.0.2_fcs.nc','longitude');
time_Jan19 = ncread('SWOT_L3_LR_SSH_Expert_433_v1.0.2_fcs.nc','time');
time_Jan19 = time_Jan19/60/60/24 + datenum(2000,1,1,0,0,0) - datenum(2024,1,1,0,0,0);
ssha_Jan19 = ncread('SWOT_L3_LR_SSH_Expert_433_v1.0.2_fcs.nc','ssha');
mdt_Jan19 = ncread('SWOT_L3_LR_SSH_Expert_433_v1.0.2_fcs.nc','mdt');
ugeo_Jan19 = ncread('SWOT_L3_LR_SSH_Expert_433_v1.0.2_fcs.nc','ugos');
vgeo_Jan19 = ncread('SWOT_L3_LR_SSH_Expert_433_v1.0.2_fcs.nc','vgos');
speed_Jan19 = sqrt( ugeo_Jan19.^2 + vgeo_Jan19.^2 );
adt_Jan19 = ssha_Jan19 + mdt_Jan19;
cd /Users/andy/Work/Manuscripts/Andy_FOCUS/

% Load the glider positions
load /Users/andy/Work/SWOT/Benoit/GliderData/SG681/possies681.mat
load /Users/andy/Work/SWOT/Benoit/GliderData/SG659/possies659.mat

end


%%% Now make the figure %%%

figure(77); close
hh = figure(77); set(hh,'Position',[ 9         123        1274         743])

pl1 = subplot(2,2,1)

kk = find(loneke>=130 & loneke<=170);
ll = find(lateke<=-40 & lateke>=-65);
contourf(loneke(kk),lateke(ll),log10(eke(kk,ll))',50,'Color','none'); shading flat
axis xy; colorbar
caxis([-2.5 -1])
ccc = cmocean('solar'); 
set(gca,'ColorMap',ccc)
axis([142 163 -59 -47])

hold on

kk = find(lonbathy>=140 & lonbathy<=165);
ll = find(latbathy<=-45 & latbathy>=-60);
%contour(lonbathy(kk),latbathy(ll),bathy(kk,ll)',[-6000:1000:-500],'Color',[0.8 0.8 0.8]);
contour(lonbathy(kk),latbathy(ll),bathy(kk,ll)',[-6000:1000:-500],'Color',[0.2 0.2 0.2],'LineWidth',1);
hold on

% Plot the SWOT position
plot(lon_Jan19(5,:),lat_Jan19(5,:),'Color','w','LineWidth',2);
leftedge = polyfit(lon_Jan19(5,518:618),lat_Jan19(5,518:618),2); xxle = linspace(lon_Jan19(5,end),158,101); yyle = leftedge(1)*xxle.^2 + leftedge(2)*xxle + leftedge(3); plot(xxle,yyle,'Color','w','LineWidth',2)
plot(lon_Jan19(27,:),lat_Jan19(27,:),'Color','w','LineWidth',2);
leftedge = polyfit(lon_Jan19(27,518:618),lat_Jan19(27,518:618),2); xxle = linspace(lon_Jan19(27,end),158,101); yyle = leftedge(1)*xxle.^2 + leftedge(2)*xxle + leftedge(3); plot(xxle,yyle,'Color','w','LineWidth',2)
plot(lon_Jan19(41,:),lat_Jan19(41,:),'Color','w','LineWidth',2);
leftedge = polyfit(lon_Jan19(41,518:618),lat_Jan19(41,518:618),2); xxle = linspace(lon_Jan19(41,end),158,101); yyle = leftedge(1)*xxle.^2 + leftedge(2)*xxle + leftedge(3); plot(xxle,yyle,'Color','w','LineWidth',2)
plot(lon_Jan19(66,:),lat_Jan19(66,:),'Color','w','LineWidth',2);
leftedge = polyfit(lon_Jan19(66,518:618),lat_Jan19(66,518:618),2); xxle = linspace(lon_Jan19(66,end),159,101); yyle = leftedge(1)*xxle.^2 + leftedge(2)*xxle + leftedge(3); plot(xxle,yyle,'Color','w','LineWidth',2)
line([150 157 157 150 150],[-55 -55 -52 -52 -55],'Color','w','LineWidth',2,'LineStyle','--')


set(gca,'YTick',[-58:2:-46],'YTickLabel',['58';'56';'54';'52';'50';'48';'46'])
xlabel('Longitude (°E)','FontSize',18)
ylabel('Latitude (°S)','FontSize',18)
set(gca,'FontSize',16);


pl2 = subplot(2,2,2)

% Plot the position of the PF
kk = find(sshlon>=140 & sshlon<=165);
ll = find(sshlat<=-45 & sshlat>=-60);
mm = find(sshtime==9);

contourf(sshlon(kk),sshlat(ll),sshadt(kk,ll,mm)',[-1.2:.05:1],'Color','none'); hold on
contour(sshlon(kk),sshlat(ll),sshadt(kk,ll,mm)',[-1.2:.1:1],'Color',[0.4 0.4 0.4]); hold on
contour(sshlon(kk),sshlat(ll),sshadt(kk,ll,mm)',[-.49 -.49],'Color',[1 1 1],'LineWidth',2); hold on
axis([145 160 -58 -48])
ccc = haxby(64);
set(gca,'ColorMap',ccc)
colorbar

line([150 157 157 150 150],[-55 -55 -52 -52 -55],'Color','k','LineWidth',2,'LineStyle','--')

plot(lon_659(5:5:end),lat_659(5:5:end),'s-','MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor',[0 171 142]/256,'Color','k')
plot(lon_681(5:5:end),lat_681(5:5:end),'o-','MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor',[138 147 255]/256,'Color','k')

plot(lon_659(3),lat_659(3),'s','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','w')
plot(lon_659(end),lat_659(end),'s','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','r')

plot(lon_681(3),lat_681(3),'o','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','w')

set(gca,'YTick',[-58:2:-48],'YTickLabel',['58';'56';'54';'52';'50';'48'])
xlabel('Longitude (°E)','FontSize',18)
ylabel('Latitude (°S)','FontSize',18)
set(gca,'FontSize',16);



pl3 = subplot(2,2,3)
% January 9

mm = find(sshtime==9);
choose_time = 8;

imagesc(sshlon,sshlat,sshsla(:,:,mm)'); shading flat; axis xy
hold on
pcolor(lon_Jan9,lat_Jan9,ssha_Jan9(:,:,choose_time)); shading flat; axis xy
axis([150 157 -55 -52])
ccc = cmocean('balance'); 
set(gca,'ColorMap',ccc)
caxis([-1 0.5])
set(gca,'FontSize',16); colorbar
set(gca,'YTick',[-55:.5:-52],'YTickLabel',['55.0';'54.5';'54.0';'53.5';'53.0';'52.5';'52.0'])
xlabel('Longitude (°E)','FontSize',18)
ylabel('Latitude (°S)','FontSize',18)


pl4 = subplot(2,2,4)
% January 19

mm = find(sshtime==19);
choose_time = 9;

imagesc(sshlon,sshlat,sshsla(:,:,mm)'); shading flat; axis xy
hold on
pcolor(lon_Jan19,lat_Jan19,ssha_Jan19(:,:,choose_time)); shading flat; axis xy
axis([150 157 -55 -52])
ccc = cmocean('balance'); 
set(gca,'ColorMap',ccc)
caxis([-1 0.5])
set(gca,'FontSize',16); colorbar
set(gca,'YTick',[-55:.5:-52],'YTickLabel',['55.0';'54.5';'54.0';'53.5';'53.0';'52.5';'52.0'])
xlabel('Longitude (°E)','FontSize',18)
ylabel('Latitude (°S)','FontSize',18)

set(pl1,'Position',[0.05 0.445 0.48 0.54])
set(pl2,'Position',[0.64 0.445 0.30 0.54])
set(pl3,'Position',[0.09 0.064 0.33 0.30])
set(pl4,'Position',[0.57 0.064 0.33 0.30])


