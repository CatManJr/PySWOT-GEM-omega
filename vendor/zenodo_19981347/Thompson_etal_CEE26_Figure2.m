% This figure shows the relationship between SSH and dynamic topography
% calcuated along the path of the glider
%
% Last updated: May 2, 2026

% Add path for plotting and for the Gibbs Seawater Toolbox
addpath /Users/andy/Work/Matlab/colmaps/
addpath /Users/andy/Work/Matlab/cmocean_v1.4/cmocean/
addpath /Users/andy/Work/Matlab/seawater_ver3_3.1/
addpath /Users/andy/Work/Matlab/gsw_matlab_v3_06_11/ 
addpath /Users/andy/Work/Matlab/gsw_matlab_v3_06_11/library/
addpath /Users/andy/Work/Matlab/gsw_matlab_v3_06_11/thermodynamics_from_t/


% The glider data can be uploaded using the Zenodod site provided in
% Thompson et al. (Comm. Earth Env., 2026).
yn = input('Do you need to load the glider data (y==1, n==0)?  ');
if yn==1
    
    %%%%%% Load individual glider profiles %%%%%%%

    %%% Load SG681 data %%%
    
    % Make a matrix that holds lon, lat, time, depth, temperature, salinity
    % potential density
    
    cd /Users/andy/Work/SWOT/Benoit/GliderData/SG681/
    sg681mat = NaN*ones(557,7,50000);
    'Load data from SG681'
    
    for i=[2:9]

        holder =  ncread(['p681000' num2str(i) '.nc'],'longitude'); nn = length(holder);
        sg681mat(i,1,1:nn) = ncread(['p681000' num2str(i) '.nc'],'longitude');
        sg681mat(i,2,1:nn) = ncread(['p681000' num2str(i) '.nc'],'latitude');
        sg681mat(i,3,1:nn) = ncread(['p681000' num2str(i) '.nc'],'ctd_time');
        sg681mat(i,4,1:nn) = ncread(['p681000' num2str(i) '.nc'],'ctd_depth'); 
        sg681mat(i,5,1:nn) = ncread(['p681000' num2str(i) '.nc'],'theta');
        sg681mat(i,6,1:nn) = ncread(['p681000' num2str(i) '.nc'],'salinity'); 
        sg681mat(i,7,1:nn) = ncread(['p681000' num2str(i) '.nc'],'sigma_theta');
        i
        
    end
    for i=[10:99]
    
        holder =  ncread(['p68100' num2str(i) '.nc'],'longitude'); nn = length(holder);
        sg681mat(i,1,1:nn) = ncread(['p68100' num2str(i) '.nc'],'longitude');
        sg681mat(i,2,1:nn) = ncread(['p68100' num2str(i) '.nc'],'latitude');
        sg681mat(i,3,1:nn) = ncread(['p68100' num2str(i) '.nc'],'ctd_time');
        sg681mat(i,4,1:nn) = ncread(['p68100' num2str(i) '.nc'],'ctd_depth'); 
        sg681mat(i,5,1:nn) = ncread(['p68100' num2str(i) '.nc'],'theta');
        sg681mat(i,6,1:nn) = ncread(['p68100' num2str(i) '.nc'],'salinity'); 
        sg681mat(i,7,1:nn) = ncread(['p68100' num2str(i) '.nc'],'sigma_theta');
        i
        
    end
    for i=[100:198 200:393]    
        holder =  ncread(['p6810' num2str(i) '.nc'],'longitude'); nn = length(holder);
        sg681mat(i,1,1:nn) = ncread(['p6810' num2str(i) '.nc'],'longitude');
        sg681mat(i,2,1:nn) = ncread(['p6810' num2str(i) '.nc'],'latitude');
        sg681mat(i,3,1:nn) = ncread(['p6810' num2str(i) '.nc'],'ctd_time');
        sg681mat(i,4,1:nn) = ncread(['p6810' num2str(i) '.nc'],'ctd_depth'); 
        sg681mat(i,5,1:nn) = ncread(['p6810' num2str(i) '.nc'],'theta');
        sg681mat(i,6,1:nn) = ncread(['p6810' num2str(i) '.nc'],'salinity'); 
        sg681mat(i,7,1:nn) = ncread(['p6810' num2str(i) '.nc'],'sigma_theta');
        i
        
    end
     
    sg681_lon = NaN*ones(1,786);
    sg681_lat = NaN*ones(1,786);
    sg681_time = NaN*ones(1,786);
    
    % Determine mean longitude, latitude and time if downcast and upcast
    for i=2:393
        
        kk = find( squeeze(sg681mat(i,4,:)) == max(squeeze(sg681mat(i,4,:))) );
        mm = 2*i-1;
        sg681_lon(mm) = nanmean(squeeze(sg681mat(i,1,1:kk)));
        sg681_lat(mm) = nanmean(squeeze(sg681mat(i,2,1:kk)));
        sg681_time(mm) = nanmean(squeeze(sg681mat(i,3,1:kk)));
        
        mm = 2*i;
        sg681_lon(mm) = nanmean(squeeze(sg681mat(i,1,kk+1:end)));
        sg681_lat(mm) = nanmean(squeeze(sg681mat(i,2,kk+1:end)));
        sg681_time(mm) = nanmean(squeeze(sg681mat(i,3,kk+1:end)));
        
        i
    end
    sg681_time = sg681_time/60/60/24 + datenum(1970,1,1,0,0,0) - datenum(2024,1,1,0,0,0);

    
    sg681_temp = NaN*ones(786,490);
    sg681_salt = NaN*ones(786,490);
    sg681_potden = NaN*ones(786,490);
    
    % Grid the data at 2-m intervals to start
    for i=2:393
        
        kk = find( squeeze(sg681mat(i,4,:)) == max(squeeze(sg681mat(i,4,:))) ,1);
        mm = 2*i-1;
        sg681_theta(mm,1:490) = NaN;
        for j=1:490
            ll = find(squeeze(sg681mat(i,4,1:kk))>(2*j+8) & squeeze(sg681mat(i,4,1:kk))<(2*j+10) );
            sg681_theta(mm,j) = nanmean(squeeze(sg681mat(i,5,ll)));
            sg681_salt(mm,j) = nanmean(squeeze(sg681mat(i,6,ll)));
            sg681_potden(mm,j) = nanmean(squeeze(sg681mat(i,7,ll)));
        end
        
        mm = 2*i;
        sg681_theta(mm,1:490) = NaN;
        for j=1:490
            ll = find(squeeze(sg681mat(i,4,kk+1:end))>(2*j+8) & squeeze(sg681mat(i,4,kk+1:end))<(2*j+10) );
            sg681_theta(mm,j) = nanmean(squeeze(sg681mat(i,5,kk+ll)));
            sg681_salt(mm,j) = nanmean(squeeze(sg681mat(i,6,kk+ll)));
            sg681_potden(mm,j) = nanmean(squeeze(sg681mat(i,7,kk+ll)));
        end
        
        i
    end
    depthints = 11:2:989;
    
    %%% Error in dive 294 -- average dives 293 and 295 %%%
    sg681_theta(294,:) = (sg681_theta(293,:)+sg681_theta(295,:))/2;
    sg681_salt(294,:) = (sg681_salt(293,:)+sg681_salt(295,:))/2;
    sg681_potden(294,:) = (sg681_potden(293,:)+sg681_potden(295,:))/2;
    
    
    'Do some smoothing on SG681'
    sg681_potden_sm = NaN*ones(size(sg681_potden));
    kk = find(~isnan(sg681_time));
    for i=1:length(depthints)
        sg681_potden_sm(kk,i) = smooth(sg681_potden(kk,i));
        sg681_time_sm = sg681_time(kk);
        i
    end
    
    %%% Load SG659 data %%%
    
    % Make a matrix that holds lon, lat, time, depth, temperature, salinity
    % potential density
    
    cd /Users/andy/Work/SWOT/Benoit/GliderData/SG659/
    sg659mat = NaN*ones(594,7,50000);
    'Load data from SG681'

    
    %for i=[29:77 90:99]
    for i=[3:9]

        holder =  ncread(['p659000' num2str(i) '.nc'],'longitude'); nn = length(holder);
        sg659mat(i,1,1:nn) = ncread(['p659000' num2str(i) '.nc'],'longitude');
        sg659mat(i,2,1:nn) = ncread(['p659000' num2str(i) '.nc'],'latitude');
        sg659mat(i,3,1:nn) = ncread(['p659000' num2str(i) '.nc'],'ctd_time');
        sg659mat(i,4,1:nn) = ncread(['p659000' num2str(i) '.nc'],'ctd_depth'); 
        sg659mat(i,5,1:nn) = ncread(['p659000' num2str(i) '.nc'],'theta');
        sg659mat(i,6,1:nn) = ncread(['p659000' num2str(i) '.nc'],'salinity'); 
        sg659mat(i,7,1:nn) = ncread(['p659000' num2str(i) '.nc'],'sigma_theta');
        i
        
    end
    for i=[10:99]
    
        holder =  ncread(['p65900' num2str(i) '.nc'],'longitude'); nn = length(holder);
        sg659mat(i,1,1:nn) = ncread(['p65900' num2str(i) '.nc'],'longitude');
        sg659mat(i,2,1:nn) = ncread(['p65900' num2str(i) '.nc'],'latitude');
        sg659mat(i,3,1:nn) = ncread(['p65900' num2str(i) '.nc'],'ctd_time');
        sg659mat(i,4,1:nn) = ncread(['p65900' num2str(i) '.nc'],'ctd_depth'); 
        sg659mat(i,5,1:nn) = ncread(['p65900' num2str(i) '.nc'],'theta');
        sg659mat(i,6,1:nn) = ncread(['p65900' num2str(i) '.nc'],'salinity'); 
        sg659mat(i,7,1:nn) = ncread(['p65900' num2str(i) '.nc'],'sigma_theta');
        i
        
    end
    for i=[100:576 578:594]
        holder =  ncread(['p6590' num2str(i) '.nc'],'longitude'); nn = length(holder);
        sg659mat(i,1,1:nn) = ncread(['p6590' num2str(i) '.nc'],'longitude');
        sg659mat(i,2,1:nn) = ncread(['p6590' num2str(i) '.nc'],'latitude');
        sg659mat(i,3,1:nn) = ncread(['p6590' num2str(i) '.nc'],'ctd_time');
        sg659mat(i,4,1:nn) = ncread(['p6590' num2str(i) '.nc'],'ctd_depth'); 
        sg659mat(i,5,1:nn) = ncread(['p6590' num2str(i) '.nc'],'theta');
        sg659mat(i,6,1:nn) = ncread(['p6590' num2str(i) '.nc'],'salinity'); 
        sg659mat(i,7,1:nn) = ncread(['p6590' num2str(i) '.nc'],'sigma_theta');
        i
        
    end
     
    sg659_lon = NaN*ones(1,1188);
    sg659_lat = NaN*ones(1,1188);
    sg659_time = NaN*ones(1,1188);
    
    % Determine mean longitude, latitude and time if downcast and upcast
    for i=3:594
        
        kk = find( squeeze(sg659mat(i,4,:)) == max(squeeze(sg659mat(i,4,:))) );
        mm = 2*i-1;
        sg659_lon(mm) = nanmean(squeeze(sg659mat(i,1,1:kk)));
        sg659_lat(mm) = nanmean(squeeze(sg659mat(i,2,1:kk)));
        sg659_time(mm) = nanmean(squeeze(sg659mat(i,3,1:kk)));
        
        mm = 2*i;
        sg659_lon(mm) = nanmean(squeeze(sg659mat(i,1,kk+1:end)));
        sg659_lat(mm) = nanmean(squeeze(sg659mat(i,2,kk+1:end)));
        sg659_time(mm) = nanmean(squeeze(sg659mat(i,3,kk+1:end)));
        
        i
    end
    sg659_time = sg659_time/60/60/24 + datenum(1970,1,1,0,0,0) - datenum(2024,1,1,0,0,0);

    
    sg659_temp = NaN*ones(1188,490);
    sg659_salt = NaN*ones(1188,490);
    sg659_potden = NaN*ones(1188,490);
    
    % Grid the data at 2-m intervals to start
    for i=2:594
        
        kk = find( squeeze(sg659mat(i,4,:)) == max(squeeze(sg659mat(i,4,:))) ,1);
        mm = 2*i-1;
        sg659_theta(mm,1:490) = NaN;
        for j=1:490
            ll = find(squeeze(sg659mat(i,4,1:kk))>(2*j+8) & squeeze(sg659mat(i,4,1:kk))<(2*j+10) );
            sg659_theta(mm,j) = nanmean(squeeze(sg659mat(i,5,ll)));
            sg659_salt(mm,j) = nanmean(squeeze(sg659mat(i,6,ll)));
            sg659_potden(mm,j) = nanmean(squeeze(sg659mat(i,7,ll)));
        end
        
        mm = 2*i;
        sg659_theta(mm,1:490) = NaN;
        for j=1:490
            ll = find(squeeze(sg659mat(i,4,kk+1:end))>(2*j+8) & squeeze(sg659mat(i,4,kk+1:end))<(2*j+10) );
            sg659_theta(mm,j) = nanmean(squeeze(sg659mat(i,5,kk+ll)));
            sg659_salt(mm,j) = nanmean(squeeze(sg659mat(i,6,kk+ll)));
            sg659_potden(mm,j) = nanmean(squeeze(sg659mat(i,7,kk+ll)));
        end
        
        i
    end
             
end

%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     %%%%%%%% Load gridded SG659 data %%%%%%%%
% 
%     'Load gridded SG659 data'
%     cd /Users/andy/Work/SWOT/Benoit/GliderData/SG659/
% 
%     temp659 = ncread('SG659_FOCUS_20240807.nc','ConsTemperature_degC');
%     salt659 = ncread('SG659_FOCUS_20240807.nc','AbsSalinity_gperkg');
%     oxy659 = ncread('SG659_FOCUS_20240807.nc','DissolvedOxygen_umolperkg');
%     bbp659 = ncread('SG659_FOCUS_20240807.nc','Backscatter700nm_permeter');
% 
%     lon659 = ncread('SG659_FOCUS_20240807.nc','Longitude_degE');
%     lat659 = ncread('SG659_FOCUS_20240807.nc','Latitude_degN');
%     press659 = ncread('SG659_FOCUS_20240807.nc','Pressure_dbar');
%     dist659 = ncread('SG659_FOCUS_20240807.nc','HorzDist_m');
%     time659 = ncread('SG659_FOCUS_20240807.nc','Time_JulDay');
% 
%     sigma0659 = gsw_sigma0(salt659,temp659);
    
    
    %%%%%%%%%% Map SG681 position to DUACS SSH %%%%%%%%%%%%

    ' Interpolate DUACS SSH to SG681'
    cd /Users/andy/Work/SWOT/Benoit/SSHData/
    sshtime = ncread('cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.25deg_P1D_1731273896002.nc','time');
    sshlon = ncread('cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.25deg_P1D_1731273896002.nc','longitude');
    sshlat = ncread('cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.25deg_P1D_1731273896002.nc','latitude');
    sshadt = ncread('cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.25deg_P1D_1731273896002.nc','adt');
    sshsla = ncread('cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.25deg_P1D_1731273896002.nc','sla');

    sshtime = (sshtime/60/60/24) + datenum(1970,1,1,0,0,0) - datenum(2024,1,1,0,0,0);

    clear interp_adt_6* interp_sla_6*

    for i=1:length(sg681_time)
    
        if ~isnan(sg681_time(i))
        kk = find(sshtime>sg681_time(i),1);
        diff = sshtime(kk)-sg681_time(i);
        mergeadt = (1-diff)*sshadt(:,:,kk) + diff*sshadt(:,:,kk-1);
        mergesla = (1-diff)*sshsla(:,:,kk) + diff*sshsla(:,:,kk-1);
        interp_adt_681(i) = interp2(sshlat,sshlon,mergeadt,sg681_lat(i),sg681_lon(i));
        interp_sla_681(i) = interp2(sshlat,sshlon,mergesla,sg681_lat(i),sg681_lon(i));
        i
        else
            interp_adt_681(i) = NaN;
            interp_sla_681(i) = NaN;
        end

    end
    
     %%%%%%%%%% Map SG681 position to SWOT SSH on January 8 %%%%%%%%%%%%

    'Interpolate SWOT SSH to SG681 on January 8'
    cd /Users/andy/Work/SWOT/Benoit/SWOTMaps/Expert/

    lat = ncread('SWOT_L3_LR_SSH_Expert_127_v1.0.2_fcs.nc','latitude');
    lon = ncread('SWOT_L3_LR_SSH_Expert_127_v1.0.2_fcs.nc','longitude');
    time = ncread('SWOT_L3_LR_SSH_Expert_127_v1.0.2_fcs.nc','time');
    time = time/60/60/24 + datenum(2000,1,1,0,0,0) - datenum(2024,1,1,0,0,0);
    ssha = ncread('SWOT_L3_LR_SSH_Expert_127_v1.0.2_fcs.nc','ssha');
    mdt = ncread('SWOT_L3_LR_SSH_Expert_127_v1.0.2_fcs.nc','mdt');
    ugeo = ncread('SWOT_L3_LR_SSH_Expert_127_v1.0.2_fcs.nc','ugos');
    vgeo = ncread('SWOT_L3_LR_SSH_Expert_127_v1.0.2_fcs.nc','vgos');
    ugeo_fit13 = ncread('SWOT_L3_LR_SSH_Expert_127_v1.0.2_fcs.nc','ugos_fit_13');
    vgeo_fit13 = ncread('SWOT_L3_LR_SSH_Expert_127_v1.0.2_fcs.nc','vgos_fit_13');
    ucgos_fit13 = ncread('SWOT_L3_LR_SSH_Expert_127_v1.0.2_fcs.nc','ucgos_fit_13');
    vcgos_fit13 = ncread('SWOT_L3_LR_SSH_Expert_127_v1.0.2_fcs.nc','vcgos_fit_13');
    speed = sqrt( ugeo.^2 + vgeo.^2 );
    adt = ssha + mdt;
    cd ../../GliderData

    lonholder = lon(:);
    latholder = lat(:);
    sshholder = ssha(:,:,7)+mdt(:,:,7);
    interp_ssh_681_swot = NaN*ones(size(sg681_time));
    for i=514:535

        clear disty
        for j=1:length(lonholder)
            disty(j) = sw_dist([latholder(j) sg681_lat(i)],[lonholder(j) sg681_lon(i)],'km');
        end
        ll = find(disty == min(disty));
        interp_ssh_681_swot(i) = sshholder(ll);
        i
    end

    
    
    %%%%%%%%%% Map SG681 position to SWOT SSH on January 15 %%%%%%%%%%%%

    'Interpolate SWOT SSH to SG681 on January 15'
    cd /Users/andy/Work/SWOT/Benoit/SWOTMaps/Expert/

    lat = ncread('SWOT_L3_LR_SSH_Expert_310_v1.0.2_fcs.nc','latitude');
    lon = ncread('SWOT_L3_LR_SSH_Expert_310_v1.0.2_fcs.nc','longitude');
    time = ncread('SWOT_L3_LR_SSH_Expert_310_v1.0.2_fcs.nc','time');
    time = time/60/60/24 + datenum(2000,1,1,0,0,0) - datenum(2024,1,1,0,0,0);
    ssha = ncread('SWOT_L3_LR_SSH_Expert_310_v1.0.2_fcs.nc','ssha');
    mdt = ncread('SWOT_L3_LR_SSH_Expert_310_v1.0.2_fcs.nc','mdt');
    ugeo = ncread('SWOT_L3_LR_SSH_Expert_310_v1.0.2_fcs.nc','ugos');
    vgeo = ncread('SWOT_L3_LR_SSH_Expert_310_v1.0.2_fcs.nc','vgos');
    ugeo_fit13 = ncread('SWOT_L3_LR_SSH_Expert_310_v1.0.2_fcs.nc','ugos_fit_13');
    vgeo_fit13 = ncread('SWOT_L3_LR_SSH_Expert_310_v1.0.2_fcs.nc','vgos_fit_13');
    ucgos_fit13 = ncread('SWOT_L3_LR_SSH_Expert_310_v1.0.2_fcs.nc','ucgos_fit_13');
    vcgos_fit13 = ncread('SWOT_L3_LR_SSH_Expert_310_v1.0.2_fcs.nc','vcgos_fit_13');
    speed = sqrt( ugeo.^2 + vgeo.^2 );
    adt = ssha + mdt;
    cd ../../GliderData

    lonholder = lon(:);
    latholder = lat(:);
    sshholder = ssha(:,:,8)+mdt(:,:,8);
    for i=615:637

        clear disty
        for j=1:length(lonholder)
            disty(j) = sw_dist([latholder(j) sg681_lat(i)],[lonholder(j) sg681_lon(i)],'km');
        end
        ll = find(disty == min(disty));
        interp_ssh_681_swot(i) = sshholder(ll);
        i
    end
    
   
yn = input('Do you need to calculate the dynamic height for SG681 (y==1, n==0)?  ');
if yn==1
    
    %%%%%%%%%% Calculte the dynamic height for SG681 %%%%%%%%%%%%
    'Calculate dynamic height for SG681'
    g=9.81;
    dp = 2;
    dyn_height_681 = NaN*ones(1,length(sg681_time));
    dyn_height_681gsw = NaN*ones(length(sg681_time),490);
    depthints = 11:2:989;
    for i=1:length(sg681_time)
    
        %dyn_height_681(i) = -nansum(1./(1000+sg681_potden(i,10:450)))*dp*1020;
        kk = find(~isnan(sg681_salt(i,:))); 
        if max(depthints(kk)>900);
            dyn_height_681gsw_900(i,kk) = gsw_geo_strf_dyn_height(sg681_salt(i,kk),sg681_theta(i,kk),depthints(kk)',900);
        end
        %dyn_height_681sw(i,:) = sw_gpan(sg681_salt(i,:),sg681_theta(i,:),(11:2:989));
        i
    
    end
    kk = find(dyn_height_681gsw_900==0); dyn_height_681gsw_900(kk) = NaN;
    kk = find(dyn_height_681==0); dyn_height_681(kk) = NaN;
    
end
    
    
cd /Users/andy/Work/Manuscripts/Andy_FOCUS/

figure(37); close
hh = figure(37); set(hh,'Position',[ 63         186        1021         800])

pl1 = subplot(3,1,1)
plot(sg681_time,dyn_height_681gsw_900(:,1)*1.041/9.81-1.43,'b.-','MarkerSize',10); hold on
plot(sg681_time,interp_adt_681,'k','LineWidth',3)
plot(sg681_time(514:535),interp_ssh_681_swot(514:535)-.03,'r','LineWidth',3)
plot(sg681_time(615:637),interp_ssh_681_swot(615:637),'r','LineWidth',3)
axis([-30 25 -1.15 0.5])
set(gca,'FontSize',16)
%xlabel('Time (days since 1 January 2024)','FontSize',18)
ylabel('Height (m)','FontSize',18)

pl2 = subplot(3,1,2)
pcolor(sg681_time(5:end),depthints,sg681_potden(5:end,:)'); axis ij; colorbar; shading flat
hold on
kk = find(~isnan(sg681_time));
contour(sg681_time(kk),depthints,sg681_potden_sm(kk,:)',[26.6:.1:27.7],'Color','w')
axis([-30 25 0 950])
ccc = cmocean('dens');
set(gca,'YTick',[0:200:800])
set(gca,'ColorMap',ccc); hold on
set(gca,'FontSize',16)
xlabel('Time (days since 1 January 2024)','FontSize',18)
ylabel('Depth (m)','FontSize',18)

pl3 = subplot(3,1,3)
pcolor(lon,lat,ssha(:,:,8)+mdt(:,:,8)); shading flat; axis xy
ccc = cmocean('balance'); 
set(gca,'ColorMap',ccc); hold on
caxis([-1.2 .5])
plot(sg681_lon(615:637),sg681_lat(615:637),'ko','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y')
axis([152 155 -55 -52])
colorbar
set(gca,'FontSize',13)
set(gca,'YTick',[-55:1:-52],'YTickLabel',['55';'54';'53';'52'])
xlabel('Longitude (°E)','FontSize',14)
ylabel('Latitude (°S)','FontSize',14)

set(pl1,'Position',[0.07 0.52 0.86 0.45])
set(pl2,'Position',[0.07 0.06 0.86 0.4])
set(pl3,'Position',[0.11    0.725    0.16    0.2325])

