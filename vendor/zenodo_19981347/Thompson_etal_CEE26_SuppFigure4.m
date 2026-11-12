% This is a script to generate the density field needed for xinvert code
%
% AFT, Feb, 2025

addpath /Users/andy/Work/Matlab/colmaps/
addpath /Users/andy/Work/Matlab/cmocean_v1.4/cmocean/
addpath /Users/andy/Work/Matlab/seawater_ver3_3.1/
addpath /Users/andy/Work/Matlab/gsw_matlab_v3_06_11/ 
addpath /Users/andy/Work/Matlab/gsw_matlab_v3_06_11/library/
addpath /Users/andy/Work/Matlab/gsw_matlab_v3_06_11/thermodynamics_from_t/


yn = input('Do you need to load the SG681 data (y==1, n==0)?  ');

if yn==1
    
     %%% Load SG681 data %%%
    
    % Make a matrix that holds lon, lat, time, depth, temperature, salinity
    % potential density
    
    sg681mat = NaN*ones(557,7,50000);
    cd /Users/andy/Work/SWOT/Benoit/GliderData/SG681/
    
    %for i=[29:77 90:99]
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
    %for i=[100:118 120:135 140:160 162:176 178:198 200:246 248:393 395:577]
    %for i=[100:198 200:393 395:577]    
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

    
    sg681_theta = NaN*ones(786,490);
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
    
    cd ..
    
end


yn = input('Do you need to load the SG659 data (y==1, n==0)?  ');

if yn==1
    
     %%% Load SG659 data %%%
    
    % Make a matrix that holds lon, lat, time, depth, temperature, salinity
    % potential density
    
    sg659mat = NaN*ones(557,7,50000);
    cd /Users/andy/Work/SWOT/Benoit/GliderData/SG659/
    
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
    %for i=[100:118 120:135 140:160 162:176 178:198 200:246 248:393 395:577]
    %for i=[100:198 200:393 395:577]    
    for i=[100:393]    
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
     
    sg659_lon = NaN*ones(1,786);
    sg659_lat = NaN*ones(1,786);
    sg659_time = NaN*ones(1,786);
    
    % Determine mean longitude, latitude and time if downcast and upcast
    for i=2:393
        
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

    
    sg659_temp = NaN*ones(786,490);
    sg659_salt = NaN*ones(786,490);
    sg659_potden = NaN*ones(786,490);
    
    % Grid the data at 2-m intervals to start
    for i=2:393
        
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
    
    cd ..
    
end

sg659_salt_intercal = sg659_salt - 0.27;
    
%%%%%%%%%% Map SG681 position to DUACS SSH %%%%%%%%%%%%

cd /Users/andy/Work/SWOT/Benoit/SSHData/
sshtime = ncread('cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.25deg_P1D_1731273896002.nc','time');
sshlon = ncread('cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.25deg_P1D_1731273896002.nc','longitude');
sshlat = ncread('cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.25deg_P1D_1731273896002.nc','latitude');
sshadt = ncread('cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.25deg_P1D_1731273896002.nc','adt');
sshsla = ncread('cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.25deg_P1D_1731273896002.nc','sla');

sshtime = (sshtime/60/60/24) + datenum(1970,1,1,0,0,0) - datenum(2024,1,1,0,0,0);

clear interp_adt_681 interp_sla_681

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

%%%%%%%%%% Map SG659 position to DUACS SSH %%%%%%%%%%%%

cd /Users/andy/Work/SWOT/Benoit/SSHData/
sshtime = ncread('cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.25deg_P1D_1731273896002.nc','time');
sshlon = ncread('cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.25deg_P1D_1731273896002.nc','longitude');
sshlat = ncread('cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.25deg_P1D_1731273896002.nc','latitude');
sshadt = ncread('cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.25deg_P1D_1731273896002.nc','adt');
sshsla = ncread('cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.25deg_P1D_1731273896002.nc','sla');

sshtime = (sshtime/60/60/24) + datenum(1970,1,1,0,0,0) - datenum(2024,1,1,0,0,0);

clear interp_adt_659 interp_sla_659

for i=1:length(sg659_time)
    
    if ~isnan(sg659_time(i))
    kk = find(sshtime>sg659_time(i),1);
    diff = sshtime(kk)-sg659_time(i);
    mergeadt = (1-diff)*sshadt(:,:,kk) + diff*sshadt(:,:,kk-1);
    mergesla = (1-diff)*sshsla(:,:,kk) + diff*sshsla(:,:,kk-1);
    interp_adt_659(i) = interp2(sshlat,sshlon,mergeadt,sg659_lat(i),sg659_lon(i));
    interp_sla_659(i) = interp2(sshlat,sshlon,mergesla,sg659_lat(i),sg659_lon(i));
    i
    else
        interp_adt_659(i) = NaN;
        interp_sla_659(i) = NaN;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a density/SSH relationship for SG681

yn = input('Do you need to define the density-SSH relationship (y==1, n==0)?  ');
if yn==1

% Find a stratification for each SSH bin

    sshmin = -1.10;
    sshmax = 0.20;
    %slamin = -0.8;
    %slamax = 0.8;
    sshbinedges = sshmin:.02:sshmax;
%    sshbinedges = slamin:.025:slamax;
    midsshbin = (sshbinedges(2:end)+sshbinedges(1:end-1))/2;
    
    for i=1:length(sshbinedges)-1

        if i==1
            kk = find(interp_adt_681<sshbinedges(i+1));
        elseif i==length(sshbinedges-1)
            kk = find(interp_adt_681>sshbinedges(i));
        else
            kk = find(interp_adt_681>sshbinedges(i) & interp_adt_681<sshbinedges(i+1));
        end
        
        for j=1:length(kk); sg681_potden(kk(j),:) = sort(sg681_potden(kk(j),:)); end
        
        nn_profs(i) = length(kk); ['No. profs = ' num2str(nn_profs(i))]
        mean_day(i) = nanmean(sg681_time(kk)); ['Mean day = ' num2str(mean_day(i))]
        std_time(i) = nanstd(sg681_time(kk)); ['Std dev time = ' num2str(std_time(i))]
        
        figure(44); hold off
        for j=1:length(kk)
            %plot(sigma0659(:,kk(j)),press659); hold on; axis ij
            scatter(squeeze(sg681_potden((kk(j)),:)),11:2:989,8,sg681_time((kk(j)))*ones(490,1),'filled'); hold on; axis ij
            %scatter(squeeze(sigma0659(:,subindy(kk(j)))),press659,8,time659(subindy(kk(j)))*ones(491,1),'filled'); hold on; axis ij
            caxis([-35 35])
            colorbar
        end
        hold on; plot(nanmean(sg681_potden(kk,:),1),11:2:989,'LineWidth',2,'Color','b')
        
        uppertime = mean_day(i)+std_time(i);
        lowertime = mean_day(i)-std_time(i);
        ll = find(sg681_time(kk) > lowertime & sg681_time(kk) < uppertime)

        nn_profs_kept(i) = length(ll); ['No. profs kept = ' num2str(nn_profs_kept(i))]
        mean_day_kept(i) = nanmean(sg681_time(kk(ll))); ['Mean day kept = ' num2str(mean_day_kept(i))]
        std_time_kept(i) = nanstd(sg681_time(kk(ll))); ['Std dev time kepth = ' num2str(std_time_kept(i))]

        if length(kk)==1
            sshbindens_681_v2(:,i) = (sg681_potden(kk,:));
        end
        if i==33
            sshbindens_681_v2(:,i) = nanmean(sg681_potden(kk(4:6),:),1);    
            nn_profs_kept(i) = 3; ['No. profs kept = ' num2str(3)]
            mean_day_kept(i) = nanmean(sg681_time(kk(4:6))); ['Mean day kept = ' num2str(mean_day_kept(i))]
            std_time_kept(i) = nanstd(sg681_time(kk(4:6))); ['Std dev time kepth = ' num2str(std_time_kept(i))]
        elseif i==40
            sshbindens_681_v2(:,i) = nanmean(sg681_potden(kk(7:12),:),1);    
            nn_profs_kept(i) = 6; ['No. profs kept = ' num2str(3)]
            mean_day_kept(i) = nanmean(sg681_time(kk(7:12))); ['Mean day kept = ' num2str(mean_day_kept(i))]
            std_time_kept(i) = nanstd(sg681_time(kk(7:12))); ['Std dev time kepth = ' num2str(std_time_kept(i))]
        else
            sshbindens_681_v2(:,i) = nanmean(sg681_potden(kk(ll),:),1);
        end
        hold on; plot(sshbindens_681_v2(:,i),11:2:989,'LineWidth',2,'Color','k')
        
        i
        pause
    end
    
   sshbindens_681_v2(:,8) = (sshbindens_681_v2(:,7)+sshbindens_681_v2(:,9))/2;

    
    % Smooth the density profiles vertically
    sshbindens_681_v2_sm = NaN*ones(size(sshbindens_681_v2));
    for i=1:length(midsshbin)
        %sshbindens_681_sm(:,i) = smooth(sshbindens_681_v2(:,i));
        sshbindens_681_v2_sm(:,i) = smoothdata(sshbindens_681_v2(:,i),1,'gaussian',5);

    end
    
    % Try sorting the density profiles (this is mostly to take care of
    % spurious deep points)
    figure(5)
    depthints = 11:2:989;
    for i=1:length(midsshbin)
    plot(sshbindens_681_v2_sm(:,i),depthints); hold on
    plot(sort(sshbindens_681_v2_sm(:,i)),depthints,'r'); hold off
    sshbindens_681_v2_sm_sort(:,i) = sort(sshbindens_681_v2_sm(:,i));
    i 
    %pause
    end
    
%     % Smooth the density profiles vertically
%     sshbindens_681_v2_sm_sort_horsm = NaN*ones(size(sshbindens_681_v2_sm_sort));
%     for i=1:length(depthints)
%         sshbindens_681_v2_sm_sort_horsm(i,:) = smooth(sshbindens_681_v2_sm_sort(i,:),10);
%     end

        
    
    figure(145)
    subplot(2,1,1)
    contourf(midsshbin,11:2:989,sshbindens_681_v2,50); shading flat; colorbar
    axis ij
    
    siggy = 26.6:.02:27.7;
    depthints = 11:2:989;
    depthsolve = NaN*ones(length(siggy),length(midsshbin));
    for i=1:length(siggy) 
        
%        depthsolve(i,:) = NaN*ones(1,length(midsshbin));
        for j=1:length(midsshbin)
            kk = find(~isnan(sshbindens_681_v2_sm(:,j)));
            if ~isempty(kk)
                depthsolve(i,j) = interp1(sshbindens_681_v2_sm(kk,j),depthints(kk),siggy(i));
            end
        end
        i
    end
    
    figure(83); close
    hh = figure(83)
    pcolor(midsshbin,siggy,depthsolve); shading flat; colorbar
    
    
    smoothed_depthsolve_09 = NaN*ones(length(siggy),length(midsshbin));
    smoothed_depthsolve_095 = NaN*ones(length(siggy),length(midsshbin));
    smoothed_depthsolve_098 = NaN*ones(length(siggy),length(midsshbin));
    smoothed_depthsolve_0999 = NaN*ones(length(siggy),length(midsshbin));

    for i=1:length(siggy)
        kk = find(~isnan(depthsolve(i,:)));
        if ~isempty(kk)
        curve = fit(midsshbin(kk)',depthsolve(i,kk)','smoothingspline','SmoothingParam',0.9); 
        smoothed_depthsolve_09(i,:) = curve(midsshbin);
        
        curve = fit(midsshbin(kk)',depthsolve(i,kk)','smoothingspline','SmoothingParam',0.95);         
        smoothed_depthsolve_095(i,:) = curve(midsshbin);   
        
        curve = fit(midsshbin(kk)',depthsolve(i,kk)','smoothingspline','SmoothingParam',0.98);           
        smoothed_depthsolve_098(i,:) = curve(midsshbin);
        
        curve = fit(midsshbin(kk)',depthsolve(i,kk)','smoothingspline','SmoothingParam',0.999);
        smoothed_depthsolve_0999(i,:) = curve(midsshbin);
        
        end
        i
    end
    
    
    for i=1:length(midsshbin)
    plot(siggy,smoothed_depthsolve_098(:,i)); hold on
    kk = find(smoothed_depthsolve_098(:,i) == max(smoothed_depthsolve_098(:,i)))
    if kk<56
        smoothed_depthsolve_098(kk+1:end,i)=NaN;
    end
    plot(siggy,sort(smoothed_depthsolve_098(:,i)),'r'); hold off
    smoothed_depthsolve_098_mod(:,i) = sort(smoothed_depthsolve_098(:,i));
    i
    pause
    end  
    
    for i=1:length(midsshbin)
        smoothed_depthsolve_098_mod(:,i) = smooth(smoothed_depthsolve_098_mod(:,i));
    end

    
    
    %smoothed_depthsolve = smoothed_depthsolve_09;
    %smoothed_depthsolve = smoothed_depthsolve_095;
    smoothed_depthsolve = smoothed_depthsolve_098_mod;
    %smoothed_depthsolve = smoothed_depthsolve_0999;
    
    smoothed_denssolve_098 = NaN*ones(length(midsshbin),length(depthints));
    smoothed_depthsolve(1:2,:) = NaN;
    
    for i=3:length(midsshbin)
        
        for j=1:length(depthints)
            
            smoothed_denssolve_098(i,j) = interp1(smoothed_depthsolve(3:end,i),siggy(3:end),depthints(j));
            
        end
        i
    end
    %kk = find(smoothed_denssolve_098==0); smoothed_denssolve_098(kk) = NaN;
    
    
    figure(84); close
    hh = figure(84)
    
    pl1 = subplot(2,2,1)
    pcolor(midsshbin,siggy,smoothed_depthsolve_09); shading flat; colorbar
    caxis([0 1000])
    
    pl2 = subplot(2,2,2)
    pcolor(midsshbin,siggy,smoothed_depthsolve_095); shading flat; colorbar
    caxis([0 1000])
    
    pl3 = subplot(2,2,3)
    pcolor(midsshbin,siggy,smoothed_depthsolve_098); shading flat; colorbar
    caxis([0 1000])
    
    pl4 = subplot(2,2,4)
    pcolor(midsshbin,siggy,smoothed_depthsolve_0999); shading flat; colorbar
    caxis([0 1000])
    
    
    figure(85)
    pcolor(midsshbin,depthints,smoothed_denssolve_098'); shading flat; colorbar
    caxis([min(siggy) max(siggy)]); axis ij
    
    
    %%% Final smoothing for the density profile %%%
    
    %sshbindens_681_final_curve_098 = NaN*ones(length(depthints),length(midsshbin));
    sshbindens_681_final_curve_098_v2 = NaN*ones(length(depthints),length(midsshbin));
    for i=1:length(depthints)
    figure(1); hold off
%        plot(midsshbin,sshbindens_681_v2_sm_sort_horsm(i,:),'kx'); hold on
        %plot(midsshbin,sshbindens_681_v2_sm_sort(i,:),'kx'); hold on
        plot(midsshbin,smoothed_denssolve_098(:,i),'kx'); hold on
%        curve = fit(midsshbin',sshbindens_681_v2_sm_sort_horsm(i,:)','smoothingspline','SmoothingParam',0.98)
        kk = find(~isnan(smoothed_denssolve_098(:,i)));
        %curve = fit(midsshbin(kk)',sshbindens_681_v2_sm_sort(i,kk)','smoothingspline','SmoothingParam',0.98)
        curve = fit(midsshbin(kk)',smoothed_denssolve_098(kk,i),'smoothingspline','SmoothingParam',0.98)
        plot(midsshbin,curve(midsshbin),'r')
        sshbindens_681_final_curve_098_v2(i,:) = curve(midsshbin);
    i
    end
    
    % Sort the new sshbins curves
    for i=1:length(midsshbin)
        sshbindens_681_final_curve_098_v2(:,i) = sort(sshbindens_681_final_curve_098_v2(:,i));
    end
    
    
    
    
    
    figure(90); close
    hh = figure(90); set(hh,'Position',[ 89         218        1219         648])
    
    pl1 = subplot(2,2,1)
    imagesc(midsshbin(3:126),depthints,sshbindens_681_sm_sort(:,3:126)); shading flat; colorbar
    caxis([26.6 27.7])
    axis([-1.1 0.18 0 1000])
    ccc = cmocean('dens');
    set(gca,'Colormap',ccc)
    set(gca,'FontSize',16)
    xlabel('SSH (m)','FontSize',18)
    ylabel('Depth (m)','FontSize',18)

    pl2 = subplot(2,2,2)
    pcolor(midsshbin,siggy,depthsolve); shading flat; cl1 = colorbar; axis ij
    caxis([0 1000])
    axis([-1.1 0.18 26.6 27.7])   
    set(cl1,'Direction','reverse')
    set(gca,'FontSize',16)
    xlabel('SSH (m)','FontSize',18)
    ylabel('Potential density (kg m^-^3)','FontSize',18)

    pl3 = subplot(2,2,4)
    pcolor(midsshbin,siggy,smoothed_depthsolve_098_mod); shading flat; cl2 = colorbar; axis ij
    caxis([0 1000])
    axis([-1.1 0.18 26.6 27.7])
    set(cl2,'Direction','reverse')
    set(gca,'FontSize',16)
    xlabel('SSH (m)','FontSize',18)
    ylabel('Potential density (kg m^-^3)','FontSize',18)
    
    pl4 = subplot(2,2,3)
    pcolor(midsshbin,depthints,smoothed_denssolve_098'); shading flat; colorbar
    caxis([min(siggy) max(siggy)]); axis ij
    caxis([26.6 27.7])
    axis([-1.1 0.18 0 1000])
    ccc = cmocean('dens');
    set(gca,'Colormap',ccc)
    set(gca,'FontSize',16)
    xlabel('SSH (m)','FontSize',18)
    ylabel('Depth (m)','FontSize',18)
    
    set(pl1,'Position',[0.06 0.58 0.38 0.4])
    set(pl2,'Position',[0.56 0.58 0.38 0.4])
    set(pl3,'Position',[0.56 0.08 0.38 0.4])
    set(pl4,'Position',[0.06 0.08 0.38 0.4])
    
end


figure(22); close
hh = figure(22); set(hh,'Position',[ 114   270   847   596])

pl1 = subplot(3,1,1)
plot(midsshbin,nn_profs,'b.-','LineWidth',2,'MarkerSize',18); hold on
plot(midsshbin,nn_profs_kept,'k.-','LineWidth',2,'MarkerSize',18)
set(gca,'FontSize',14)
xlabel('Sea surface height (SSH m)','FontSize',14)
ylabel('Number of profiles','FontSize',14)

pl2 = subplot(3,1,2)
plot(midsshbin,mean_day,'b.-','LineWidth',2,'MarkerSize',18); hold on
plot(midsshbin,mean_day_kept,'k.-','LineWidth',2,'MarkerSize',18)
set(gca,'FontSize',14)
xlabel('Sea surface height (SSH m)','FontSize',14)
ylabel('Mean day (days since 1 January 2024)','FontSize',14)

pl3 = subplot(3,1,3)
plot(midsshbin,std_time,'b.-','LineWidth',2,'MarkerSize',18); hold on
plot(midsshbin,std_time_kept,'k.-','LineWidth',2,'MarkerSize',18)
set(gca,'FontSize',14)
xlabel('Sea surface height (SSH m)','FontSize',14)
ylabel('Standard deviation (d)','FontSize',14)

set(pl1,'Position',[0.13 0.74 0.80 0.24])
set(pl2,'Position',[0.13 0.42 0.80 0.24])
set(pl3,'Position',[0.13 0.10 0.80 0.24])



%%% Reconstruct SG681 density %%%%
dsshbin = midsshbin(2)-midsshbin(1);
reconstruct_pd_681 = NaN*ones(length(depthints),length(sg681_time));
for i=5:786
    
    if ~isnan(interp_adt_681(i))
    
    target = interp_adt_681(i);
    kk = find(midsshbin>target,1);
    sep = (midsshbin(kk) - interp_adt_681(i))/dsshbin;
    
    if kk==1
        reconstruct_pd_681(:,i) = (1-sep)*sshbindens_681_final_curve_098(:,kk) + sep*sshbindens_681_final_curve_098(:,kk);
    else
        reconstruct_pd_681(:,i) = (1-sep)*sshbindens_681_final_curve_098(:,kk) + sep*sshbindens_681_final_curve_098(:,kk-1);
    end

    
    end
    i
    
end

%%% Reconstruct SG659 density %%%%
dsshbin = midsshbin(2)-midsshbin(1);
reconstruct_pd_659 = NaN*ones(length(depthints),length(sg659_time));
for i=5:786
    
    if ~isnan(interp_adt_659(i))
    
    target = interp_adt_659(i);
    kk = find(midsshbin>target,1);
    sep = (midsshbin(kk) - interp_adt_659(i))/dsshbin;
    
    if kk==1
        reconstruct_pd_659(:,i) = (1-sep)*sshbindens_681_final_curve_098(:,kk) + sep*sshbindens_681_final_curve_098(:,kk);
    else
        reconstruct_pd_659(:,i) = (1-sep)*sshbindens_681_final_curve_098(:,kk) + sep*sshbindens_681_final_curve_098(:,kk-1);
    end

    
    end
    i
    
end




figure(23); close
hh = figure(23); set(hh,'Position',[  57         175        1179         691])

pl1 = subplot(2,2,1)
plot(midsshbin,nn_profs_kept,'k.-','LineWidth',2,'MarkerSize',18); hold on
plot(-.48*ones(1,11),0:5:50,'k--','LineWidth',1)
axis([-1.1 0.2 0 50])
set(gca,'FontSize',14)
xlabel('Sea surface height (SSH m)','FontSize',16)
ylabel('Number of profiles','FontSize',16)


pl3 = subplot(2,2,3)
imagesc(midsshbin,depthints,sshbindens_681_final_curve_098); hold on
contour(midsshbin,depthints,sshbindens_681_final_curve_098,[26.6:.1:27.7],'Color',[.3 .3 .3],'LineWidth',.5)
colorbar
caxis([26.6 27.7])
set(gca,'FontSize',14)
xlabel('Sea surface height (SSH m)','FontSize',16)
ylabel('Depth (m)','FontSize',16)


pl2 = subplot(2,2,2)
dsshbin = midsshbin(2)-midsshbin(1);
reconstruct_pd_681 = NaN*ones(length(depthints),length(sg681_time));
for i=5:786
    
    if ~isnan(interp_adt_681(i))
    
    target = interp_adt_681(i);
    kk = find(midsshbin>target,1);
    sep = (midsshbin(kk) - interp_adt_681(i))/dsshbin;
    
    if kk==1
        reconstruct_pd_681(:,i) = (1-sep)*sshbindens_681_final_curve_098(:,kk) + sep*sshbindens_681_final_curve_098(:,kk);
    else
        reconstruct_pd_681(:,i) = (1-sep)*sshbindens_681_final_curve_098(:,kk) + sep*sshbindens_681_final_curve_098(:,kk-1);
    end

    
    end
    i
    
end

holder = reconstruct_pd_681(:,[5:396 399:593 595:end]);
imagesc(sg681_time([5:396 399:593 595:end]),depthints,holder - 27.2); hold on
%contour(sg681_time([5:396 399:593 595:end]),depthints,sg681_potden([5:396 399:593 595:end],:)',[26.6:.1:27.7],'Color',[0.5 0.5 0.5])
ccc = cmocean('dense');
ccc = cmocean('delta');
set(gca,'Colormap',ccc)
caxis([-.6 .6]); 
colorbar
axis([-30 25 0 950])
set(gca,'FontSize',14)
xlabel('Time (days since January 1, 2024)','FontSize',16)
ylabel('Depth (m)','FontSize',16)

pl4 = subplot(2,2,4)
imagesc(sg681_time([5:396 399:593 595:end]),depthints,reconstruct_pd_681(:,[5:396 399:593 595:end])-sg681_potden([5:396 399:593 595:end],:)')
ccc = cmocean('delta');
set(gca,'Colormap',ccc)
caxis([-0.3 0.3])
caxis([-.6 .6]); colorbar
axis([-30 25 0 950])
set(gca,'FontSize',14)
xlabel('Time (days since January 1, 2024)','FontSize',16)
ylabel('Depth (m)','FontSize',16)

set(pl1,'Position',[0.05 0.67 0.30 0.30])
set(pl2,'Position',[0.46 0.57 0.48 0.40])
set(pl3,'Position',[0.05 0.08 0.30 0.44])
set(pl4,'Position',[0.46 0.08 0.48 0.40])


figure(24); close
hh = figure(24); set(hh,'Position',[  82         386        1082         480])

pl1 = subplot(1,2,1)
imagesc(midsshbin,depthints,sshbindens_681_v2_sm); hold on
%contour(midsshbin,depthints,sshbindens_681_final_curve_098,[26.6:.1:27.7],'Color',[.3 .3 .3],'LineWidth',.5)
colorbar
caxis([26.6 27.7])
set(gca,'FontSize',14)
xlabel('Sea surface height (SSH m)','FontSize',16)
ylabel('Depth (m)','FontSize',16)


pl2 = subplot(1,2,2)
imagesc(midsshbin,depthints,sshbindens_681_final_curve_098); hold on
contour(midsshbin,depthints,sshbindens_681_final_curve_098,[26.6:.1:27.7],'Color',[.3 .3 .3],'LineWidth',.5)
colorbar
caxis([26.6 27.7])
set(gca,'FontSize',14)
xlabel('Sea surface height (SSH m)','FontSize',16)
ylabel('Depth (m)','FontSize',16)

set(pl1,'Position',[0.06 0.10 0.37 0.85])
set(pl2,'Position',[0.56 0.10 0.37 0.85])


figure(39); close
hh = figure(39); set(hh,'Position',[ 42          25        1039         841])

kk = ~isnan(interp_adt_681);
ll = ~isnan(sg659_theta(:,100));

depthy = 245;
%sg681_pden = sw_pden(sg681_theta,sg681_salt,ones(786,1)*depthints,0)
sg659_potden_intercal = sw_pden(sg659_salt_intercal,sg659_theta,ones(786,1)*depthints,0)-1000;

pl1 = subplot(4,2,1)
plot(sg681_time,interp_adt_681,'k','LineWidth',2)
set(gca,'FontSize',14)
ylabel('Sea surface height (m)','FontSize',14)
axis([-30 25 -1.2 .2])

pl3 = subplot(4,2,3)
%pcolor(sg681_time,depthints,sg681_potden' - nanmean(sg681_potden(:,depthy)) ); shading flat; colorbar; axis ij
pcolor(sg681_time(kk),depthints,sg681_potden(kk,:)' ); shading flat; colorbar; axis ij
ccc = cmocean('dense'); set(gca,'Colormap',ccc)
%caxis([-.6 .6])
caxis([26.4 28]); axis([-30 25 0 950])
set(gca,'FontSize',14)
ylabel('Depth (m)','FontSize',14)

pl5 = subplot(4,2,5)
%pcolor(sg681_time,depthints,reconstruct_pd_681 - nanmean(reconstruct_pd_681(depthy,:)) ); shading flat; colorbar; axis ij
pcolor(sg681_time(kk),depthints,reconstruct_pd_681(:,kk)  ); shading flat; colorbar; axis ij
ccc = cmocean('dense'); set(gca,'Colormap',ccc)
%caxis([-.6 .6])
caxis([26.4 28]); axis([-30 25 0 950])
set(gca,'FontSize',14)
ylabel('Depth (m)','FontSize',14)

pl7 = subplot(4,2,7)
pcolor(sg681_time(kk),depthints, (reconstruct_pd_681(:,kk) - nanmean(reconstruct_pd_681(depthy,kk)) ) - (sg681_potden(kk,:)' - nanmean(sg681_potden(kk,depthy))) ); shading flat; colorbar; axis ij
ccc = cmocean('delta'); set(gca,'Colormap',ccc)
caxis([-.65 .65]); axis([-30 25 0 950])
%caxis([26.5 28])
set(gca,'FontSize',14)
ylabel('Depth (m)','FontSize',14)
xlabel('Time (days)','FontSize',14)

pl2 = subplot(4,2,2)
plot(sg659_time,interp_adt_659,'k','LineWidth',2)
set(gca,'FontSize',14)
ylabel('Sea surface height (m)','FontSize',14)
axis([-30 25 -1.2 .2])

pl4 = subplot(4,2,4)
pcolor(sg659_time(ll),depthints,sg659_potden(ll,:)' ); shading flat; colorbar; axis ij
ccc = cmocean('dense'); set(gca,'Colormap',ccc)
%caxis([-.6 .6])
caxis([26.4 28]); axis([-30 25 0 950])
set(gca,'FontSize',14)
ylabel('Depth (m)','FontSize',14)

pl6 = subplot(4,2,6)
pcolor(sg659_time(ll),depthints,reconstruct_pd_659(:,ll)  ); shading flat; colorbar; axis ij
ccc = cmocean('dense'); set(gca,'Colormap',ccc)
%caxis([-.6 .6])
caxis([26.4 28]); axis([-30 25 0 950])
set(gca,'FontSize',14)
ylabel('Depth (m)','FontSize',14)

pl8 = subplot(4,2,8)
pcolor(sg659_time(ll),depthints, (reconstruct_pd_659(:,ll) - nanmean(reconstruct_pd_659(depthy,ll)) ) - (sg659_potden(ll,:)' - nanmean(sg659_potden(ll,depthy)))); shading flat; colorbar; axis ij
ccc = cmocean('delta'); set(gca,'Colormap',ccc)
caxis([-.65 .65]); axis([-30 25 0 950])
set(gca,'FontSize',14)
ylabel('Depth (m)','FontSize',14)
xlabel('Time (days)','FontSize',14)


set(pl1,'Position',[0.05 0.79 0.38 0.18])
set(pl2,'Position',[0.55 0.79 0.38 0.18])
set(pl3,'Position',[0.05 0.55 0.38 0.2])
set(pl4,'Position',[0.55 0.55 0.38 0.2])
set(pl5,'Position',[0.05 0.31 0.38 0.2])
set(pl6,'Position',[0.55 0.31 0.38 0.2])
set(pl7,'Position',[0.05 0.07 0.38 0.2])
set(pl8,'Position',[0.55 0.07 0.38 0.2])