% Supplementary figure on smoothing sensitivity
%
% AFT, September 23, 2025


cd /Users/andy/Work/SWOT/Benoit/GliderData/

% load the vertical velocities with various smoothings and density
% relationhips

wGS5DS1 = ncread('FOCUS_data_trunc_swath_155_9_GS5_SSHDens1.nc','w_vel');
wGS8DS1 = ncread('FOCUS_data_trunc_swath_155_9_GS8_SSHDens1.nc','w_vel');
wGS10DS1 = ncread('FOCUS_data_trunc_swath_155_9_GS10_SSHDens1.nc','w_vel');
wGS15DS1 = ncread('FOCUS_data_trunc_swath_155_9_GS15_SSHDens1.nc','w_vel');

wGS5DS2 = ncread('FOCUS_data_trunc_swath_155_9_GS5_SSHDens2.nc','w_vel');
wGS8DS2 = ncread('FOCUS_data_trunc_swath_155_9_GS8_SSHDens2.nc','w_vel');
wGS10DS2 = ncread('FOCUS_data_trunc_swath_155_9_GS10_SSHDens2.nc','w_vel');
wGS15DS2 = ncread('FOCUS_data_trunc_swath_155_9_GS15_SSHDens2.nc','w_vel');


% Load different SSH-Dens relationships
load sshbindens_681_final_curve_090.mat


% Set the ssh bins
sshmin = -1.10;
sshmax = 0.20;
sshbinedges = sshmin:.02:sshmax;
midsshbin = (sshbinedges(2:end)+sshbinedges(1:end-1))/2;

% Set depth and density surfaces
siggy = 26.6:.02:27.7;
depthints = 11:2:989;



figure(88); close
hh = figure(88); set(hh,'Position',[ 75     7   932   859])

pl1 = subplot(3,3,1)
imagesc(midsshbin,depthints,sshbindens_681_final_curve_098); hold on
contour(midsshbin,depthints,sshbindens_681_final_curve_098,[26.6:.1:27.7],'Color',[.3 .3 .3],'LineWidth',.5)
colorbar
caxis([26.3 27.8])
axis([-1.1 0.2 0 960])
set(gca,'FontSize',14)
xlabel('Absolute Dynamic Topography (m)','FontSize',14)
ylabel('Depth (m)','FontSize',14)
t1 = text( -1.0725,  -50.2797 ,'(a)','FontSize',14)

pl4 = subplot(3,3,4)
imagesc(midsshbin,depthints,sshbindens_681_final_curve_098_v2); hold on
contour(midsshbin,depthints,sshbindens_681_final_curve_098,[26.6:.1:27.7],'Color',[.3 .3 .3],'LineWidth',.5)
colorbar
caxis([26.3 27.8])
axis([-1.1 0.2 0 960])
set(gca,'FontSize',14)
xlabel('Absolute Dynamic Topography (m)','FontSize',14)
ylabel('Depth (m)','FontSize',14)
t4 = text( -1.0725,  -50.2797 ,'(b)','FontSize',14)

pl7 = subplot(3,3,7)
imagesc(midsshbin,depthints,sshbindens_681_final_curve_098 - sshbindens_681_final_curve_098_v2);
colorbar
axis([-1.1 0.2 0 950])
ccc = cmocean('balance');
caxis([-.08 .08])
set(gca,'Colormap',ccc)
set(gca,'FontSize',14)
xlabel('Absolute Dynamic Topography (m)','FontSize',14)
ylabel('Depth (m)','FontSize',14)
t7 = text( -1.0725,  -50.2797 ,'(c)','FontSize',14)

pl2 = subplot(3,3,2)
pcolor(lon_buoy,lat_buoy,wGS5DS1(:,:,50)*60*60*24); shading flat
caxis([-200 200])
set(gca,'YTick',[-58:-53],'YTickLabel',['58';'57';'56';'55';'54';'53'])
set(gca,'Colormap',ccc)
set(gca,'FontSize',14)
xlabel('Longitude (°E)','FontSize',14)
ylabel('Latitude (°S)','FontSize',14)
t2a = text(149.3439,  -52.3388 ,'Scale = 8 km','FontSize',14)
t2b = text(154.4834,  -58.1175 ,'(d)','FontSize',14)
title('SSH-\sigma (panel a)','FontSize',14)

pl3 = subplot(3,3,3)
pcolor(lon_buoy,lat_buoy,wGS5DS2(:,:,50)*60*60*24); shading flat
colorbar
caxis([-200 200])
set(gca,'YTick',[-58:-53],'YTickLabel',['58';'57';'56';'55';'54';'53'])
set(gca,'Colormap',ccc)
set(gca,'FontSize',14)
xlabel('Longitude (°E)','FontSize',14)
t3a = text(149.3439,  -52.3388 ,'Scale = 8 km','FontSize',14)
t3b = text(154.4834,  -58.1175 ,'(g)','FontSize',14)
%ylabel('Latitude (°S)','FontSize',14)
title('SSH-\sigma (panel b)','FontSize',14)


pl5 = subplot(3,3,5)
pcolor(lon_buoy,lat_buoy,wGS8DS1(:,:,50)*60*60*24); shading flat
caxis([-200 200])
set(gca,'YTick',[-58:-53],'YTickLabel',['58';'57';'56';'55';'54';'53'])
set(gca,'Colormap',ccc)
set(gca,'FontSize',14)
xlabel('Longitude (°E)','FontSize',14)
ylabel('Latitude (°S)','FontSize',14)
t5a = text(149.3439,  -52.3388 ,'Scale = 15 km','FontSize',14)
t5b = text(154.4834,  -58.1175 ,'(e)','FontSize',14)



pl6 = subplot(3,3,6)
pcolor(lon_buoy,lat_buoy,wGS8DS2(:,:,50)*60*60*24); shading flat
colorbar
caxis([-200 200])
set(gca,'YTick',[-58:-53],'YTickLabel',['58';'57';'56';'55';'54';'53'])
set(gca,'Colormap',ccc)
set(gca,'FontSize',14)
xlabel('Longitude (°E)','FontSize',14)
%ylabel('Latitude (°S)','FontSize',14)
t6a = text(149.3439,  -52.3388 ,'Scale = 15 km','FontSize',14)
t6b = text(154.4834,  -58.1175 ,'(h)','FontSize',14)


pl8 = subplot(3,3,8)
pcolor(lon_buoy,lat_buoy,wGS15DS1(:,:,50)*60*60*24); shading flat
caxis([-200 200])
set(gca,'YTick',[-58:-53],'YTickLabel',['58';'57';'56';'55';'54';'53'])
set(gca,'Colormap',ccc)
set(gca,'FontSize',14)
xlabel('Longitude (°E)','FontSize',14)
ylabel('Latitude (°S)','FontSize',14)
t8a = text(149.3439,  -52.3388 ,'Scale = 30 km','FontSize',14)
t8b = text(154.4834,  -58.1175 ,'(f)','FontSize',14)


pl9 = subplot(3,3,9)
pcolor(lon_buoy,lat_buoy,wGS15DS2(:,:,50)*60*60*24); shading flat
colorbar
caxis([-200 200])
set(gca,'YTick',[-58:-53],'YTickLabel',['58';'57';'56';'55';'54';'53'])
set(gca,'Colormap',ccc)
set(gca,'FontSize',14)
xlabel('Longitude (°E)','FontSize',14)
%ylabel('Latitude (°S)','FontSize',14)
t9a = text(149.3439,  -52.3388 ,'Scale = 30 km','FontSize',14)
t9b = text(154.4834,  -58.1175 ,'(i)','FontSize',14)



set(pl1,'Position',[0.06 0.72 0.25 0.25])
set(pl2,'Position',[0.44 0.72 0.22 0.25])
set(pl3,'Position',[0.70 0.72 0.22 0.25])
set(pl4,'Position',[0.06 0.39 0.25 0.25])
set(pl5,'Position',[0.44 0.39 0.22 0.25])
set(pl6,'Position',[0.70 0.39 0.22 0.25])
set(pl7,'Position',[0.06 0.06 0.25 0.25])
set(pl8,'Position',[0.44 0.06 0.22 0.25])
set(pl9,'Position',[0.70 0.06 0.22 0.25])

cd /Users/andy/Work/Manuscripts/Andy_FOCUS/


figure(89); close
hh = figure(89); set(hh,'Position',[273    70   427   796])

wedges = -300:1:300;
midwedges = (wedges(2:end)+wedges(1:end-1))/2;

nn_1 = hist(wGS5DS1(:)*60*60*24,midwedges);
nn_2 = hist(wGS10DS1(:)*60*60*24,midwedges);
nn_3 = hist(wGS15DS1(:)*60*60*24,midwedges);

nn_4 = hist(wGS8DS1(:)*60*60*24,midwedges);
nn_5 = hist(wGS8DS2(:)*60*60*24,midwedges);

pl1 = subplot(2,1,1)
semilogy(midwedges,nn_1,'LineWidth',2,'Color',[0 .38 .66]); hold on
semilogy(midwedges,nn_2,'LineWidth',2,'Color',[.52 0 .20])
semilogy(midwedges,nn_3,'LineWidth',2,'Color',[.9 .58 0])
axis([-200 200 10 1e6])
set(gca,'FontSize',14)
xlabel('Vertical velocity (m day^-^1)','FontSize',14)
ylabel('Count','FontSize',14)
ll = legend('8 km','15 km','30 km')
set(ll,'FontSize',14)
t1 = text(-192,545300,'(j)'); set(t1,'FontSize',14)


pl2 = subplot(2,1,2)
semilogy(midwedges,nn_4,'LineWidth',2,'Color','k','LineStyle','-'); hold on
semilogy(midwedges,nn_5,'LineWidth',2,'Color','k','LineStyle','--')
axis([-200 200 10 1e6])
set(gca,'FontSize',14)
xlabel('Vertical velocity (m day^-^1)','FontSize',14)
ylabel('Count','FontSize',14)
ll = legend('SSH-\sigma (a)','SSH-\sigma (b)')
set(ll,'FontSize',14)
t1 = text(-192,545300,'(k)'); set(t1,'FontSize',14)



figure(90); close
hh = figure(90)

wholder_5_100 = wGS5DS1(:,:,16);
wholder_10_100 = wGS10DS1(:,:,16);
wholder_15_100 = wGS15DS1(:,:,16);
wholder_5_500 = wGS5DS1(:,:,82);
wholder_10_500 = wGS10DS1(:,:,82);
wholder_15_500 = wGS15DS1(:,:,82);

nn_1 = hist(wholder_5_100(:)*60*60*24,midwedges);
nn_2 = hist(wholder_10_100(:)*60*60*24,midwedges);
nn_3 = hist(wholder_15_100(:)*60*60*24,midwedges);
nn_4 = hist(wholder_5_500(:)*60*60*24,midwedges);
nn_5 = hist(wholder_10_500(:)*60*60*24,midwedges);
nn_6 = hist(wholder_15_500(:)*60*60*24,midwedges);


semilogy(midwedges,nn_1,'LineWidth',2,'Color','b','LineStyle','-'); hold on
semilogy(midwedges,nn_2,'LineWidth',2,'Color','r','LineStyle','-'); hold on
semilogy(midwedges,nn_3,'LineWidth',2,'Color','g','LineStyle','-'); hold on
semilogy(midwedges,nn_4,'LineWidth',2,'Color','y','LineStyle','-'); hold on
semilogy(midwedges,nn_5,'LineWidth',2,'Color','k','LineStyle','-'); hold on
semilogy(midwedges,nn_6,'LineWidth',2,'Color','m','LineStyle','-'); hold on


