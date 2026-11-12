%%% Third version of the vertical velocity script. %%%
%%% This generates the Q-vector for all the relevant dates. %%%

% Add required folders
addpath /Users/andy/Work/Matlab/colmaps/
addpath /Users/andy/Work/Matlab/cmocean_v1.4/cmocean/
addpath /Users/andy/Work/Matlab/seawater_ver3_3.1/
addpath /Users/andy/Work/Matlab/gsw_matlab_v3_06_11/ 
addpath /Users/andy/Work/Matlab/gsw_matlab_v3_06_11/library/
addpath /Users/andy/Work/Matlab/gsw_matlab_v3_06_11/thermodynamics_from_t/


%%%% The SSH-density relationship is dervied in MockDensity.m %%%%

%%%% Here are the list of SWOT passes.  %%%%

SWOT_swath_no = input('Please enter SWOT swath no.:  ');
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

% Choosing smoothing scale
gauss_smooth = 8; %%% This is what is currently used in the Nat. Comm. paper!
%gauss_smooth = 5;
%gauss_smooth = 10;
%gauss_smooth = 15;


%%%%%% Make the plot for SSH %%%%%%    

cd ../../GliderData
load sshbindens_681_final_curve_090.mat
%sshbindens_681_final_curve_098=sshbindens_681_final_curve_098';

sshmin = -1.10;
sshmax = 0.20;
sshbinedges = sshmin:.02:sshmax;
midsshbin = (sshbinedges(2:end)+sshbinedges(1:end-1))/2;

adt = ssha(:,:,choose_time) + mdt(:,:,choose_time);
%smoothed_denssolve = sshbindens_681_final_curve_098';
smoothed_denssolve = sshbindens_681_final_curve_098_v2';
%smoothed_denssolve = smoothed_denssolve_098;

siggy = 26.6:.02:27.7;
depthints = 11:2:989;

% Create the reconstructed density field
mockdensity_3 = NaN*ones(size(adt,1),size(adt,2),490);
dsshbin = midsshbin(2)-midsshbin(1);

    for i=1:size(lon,1)
        for j=1:size(lon,2)
            
            if ~isnan(adt(i,j))
        
                kk = find(midsshbin>adt(i,j),1);
                sep = (midsshbin(kk) - adt(i,j))/dsshbin;
                if kk>1
                    mockdensity_3(i,j,:) = (1-sep)*smoothed_denssolve(kk,:) + sep*smoothed_denssolve(kk-1,:);
                elseif kk==1
                    mockdensity_3(i,j,:) = smoothed_denssolve(1,:);
                elseif isempty(kk)
                    mockdensity_3(i,j,:) = smoothed_denssolve(end,:);
                end
                
            end
            
        end
    end
    

    figure(186); close
    hh = figure(186); set(hh,'Position',[ 59         100        1295 766  ])
    
    choose_sec = 15;
    
    pl1 = subplot(2,2,1)
    pcolor(lon,lat,adt); shading flat; colorbar; 
    hold on
    hax = haxby(64); 
    set(gca,'colormap',hax)
    set(gca,'YTick',[-57:-53],'YTickLabel',['57';'56';'55';'54';'53'])
    plot(lon(choose_sec,:),lat(choose_sec,:),'k-','LineWidth',2)
    %axis([151 156 -57.7 -52.5])
    caxis([-1.2 0.3])
    set(gca,'FontSize',16)
    xlabel('Longitude (°E)','FontSize',16)
    ylabel('Latitude (°S)','FontSize',16)
        
    pl2 = subplot(2,2,3)
    pcolor(lon,lat,squeeze(mockdensity_3(:,:,10))); shading flat; colorbar; hold on
    %axis([151 156 -57.7 -52.5])
    plot(lon(choose_sec,:),lat(choose_sec,:),'k-','LineWidth',2)
    ccc = cmocean('dens');
    set(gca,'YTick',[-57:-53],'YTickLabel',['57';'56';'55';'54';'53'])
    set(gca,'Colormap',ccc)
    caxis([26.5 max(siggy)]);
    set(gca,'FontSize',16)
    xlabel('Longitude (°E)','FontSize',16)
    ylabel('Latitude (°S)','FontSize',16)
    
    %pause
    
    % Get rid of NaN's near the surface and near the bottom
    holderdens = mockdensity_3;
    latrange = find(lat(1,:)>=-58 & lat(1,:)<-52);
    
    holderdens = holderdens(:,latrange,:);
    lon_buoy = lon(:,latrange);
    lat_buoy = lat(:,latrange);
    
    mm = size(lon_buoy,1);
    nn = size(lat_buoy,2);
    oo = length(depthints);
    
    orientation = input('Is the path oriented (SW-NE; 1) or (NW-SE; 2):  ')
    if orientation==2
        lon_buoy = lon_buoy(mm:-1:1,:,:); lon_buoy = lon_buoy(:,nn:-1:1,:);
        lat_buoy = lat_buoy(mm:-1:1,:,:); lat_buoy = lat_buoy(:,nn:-1:1,:);
        holderdens = holderdens(mm:-1:1,:,:); holderdens = holderdens(:,nn:-1:1,:);
    end
    
    
    % Get rid of NaN's near the surface
    for i=1:mm
        for j=1:nn
            kk = find(isnan(holderdens(i,j,1:30)));
            if ~isempty(kk)
                maxkk = max(kk);
            %    pf = polyfit((maxkk+1:30),holderdens(i,maxkk+1:30),1);
            %    holderdens(i,j,1:maxkk) = (1:maxkk)*pf(1)+pf(2);
                holderdens(i,j,1:maxkk) = holderdens(i,j,maxkk+1);
            end
        end
    end
    
    % Get rid of NaN's near the bottom
    for i=1:mm
        for j=1:nn
            kk = find(isnan(holderdens(i,j,401:490)))+400;
            if ~isempty(kk)
                minkk = min(kk);
                %curve = fit((401:minkk-1)',holderdens(i,401:minkk-1)','linearinterp');
                %holderdens(i,minkk:490) = curve(minkk:490);
                if ~isnan(holderdens(i,j,minkk-1))
                    pf = polyfit(401:minkk-1,squeeze(holderdens(i,j,401:minkk-1))',1);
                    holderdens(i,j,minkk:490) = (minkk:490)*pf(1)+pf(2);
                else
                    holderdens(i,j,minkk:490) = NaN;
                end
            end
        end
    end
    
    % Make sure the water column is stably stratified
    for i=1:mm
        for j=1:nn
            holderdens(i,j,:) = sort(holderdens(i,j,:));
        end
    end
    
    % Calculate the horizontal average of the density field and the
    % background N^2
    for i=1:length(depthints)
        rhomean(i) = squeeze(nanmean(nanmean(holderdens(:,:,i))));
    end
    rho0 = 1000+nanmean(rhomean);
    N2mean = NaN*ones(1,490);
    N2mean(2:end-1) = (rhomean(3:end)-rhomean(1:end-2))/4/rho0*9.81;
    
    % Calculate the buoyancy field
    buoyancy = NaN*ones(size(holderdens));
    for i=1:490
        %buoyancy_2(:,:,i) = (rhomean(i) - holderdens(:,:,i))*9.81/(rhomean(i)+1000);
        buoyancy(:,:,i) = (rho0 - (1000+holderdens(:,:,i)))/rho0*9.81;
    end
    
    % Calculate the horizontal and vertical buoyancy gradients
    M2_across = NaN*ones(mm,nn,oo);
    M2_along  = NaN*ones(mm,nn,oo);
    M2_across_sm = NaN*ones(size(M2_across));
    M2_along_sm  = NaN*ones(size(M2_along));
    N2 = NaN*ones(mm,nn,oo);
    N2_sm = NaN*ones(size(N2));
    
    
    % Calculate M^2, 
    'Calculate the M2 field'
    for j=1:nn
        M2_across(2:end-1,j,:) = (buoyancy(3:end,j,:)-buoyancy(1:end-2,j,:))/2/2000;
    end
    for i=1:mm
        M2_along(i,2:end-1,:) = (buoyancy(i,3:end,:)-buoyancy(i,1:end-2,:))/2/2000;
    end
        
    % Smooth the M^2 field
    'Smooth the M2 field'
    for i=1:oo
        M2_across_sm(:,:,i) = smoothdata(M2_across(:,:,i),1,'gaussian',gauss_smooth);
        M2_across_sm(:,:,i) = smoothdata(M2_across_sm(:,:,i),2,'gaussian',gauss_smooth);
        
        M2_along_sm(:,:,i) = smoothdata(M2_along(:,:,i),1,'gaussian',gauss_smooth);
        M2_along_sm(:,:,i) = smoothdata(M2_along_sm(:,:,i),2,'gaussian',gauss_smooth);
        %M2_across_sm(:,:,i) = conv2(M2_across(:,:,i),K,'same');
    end
    
    % Put the gap back in the middle that was filled in by interpolation
    kk = find(isnan(buoyancy(:,1,200)));
    M2_along_sm(kk,:,:) = NaN;
    M2_across_sm(kk,:,:) = NaN;
    

    % Calculate N^2
    'Calculate the N^2 field'
    N2(:,:,2:end-1) = -(buoyancy(:,:,3:end)-buoyancy(:,:,1:end-2))/2/2;
    N2mean_sm = smooth(N2mean,12);

    % Smooth the N^2 field
    %'Smooth the N2 field'
    %for i=1:mm
    %    for j=1:nn
    %        N2_sm(i,j,:) = smooth(N2(i,j,:),20);
    %    endu
    %    i
    %end

    
    %%%%% Plot a summary figure %%%%%
    
    figure(286); close
    hh = figure(286); set(hh,'Position',[ 59         100        1295 766  ])
    
    %choose_sec = input('Choose a longitude position for along-swath transects:  ');
    choose_sec = 15;
    
    pl1 = subplot(2,2,1)
    pcolor(lon,lat,adt); shading flat; colorbar; 
    hold on
    hax = haxby(64); 
    set(gca,'colormap',hax)
    set(gca,'YTick',[-57:-53],'YTickLabel',['57';'56';'55';'54';'53'])
    plot(lon_buoy(choose_sec,:),lat_buoy(choose_sec,:),'k-','LineWidth',2)
    axis([151 159 -58 -52])
    caxis([-1.2 0.3])
    set(gca,'FontSize',16)
    xlabel('Longitude (°E)','FontSize',16)
    ylabel('Latitude (°S)','FontSize',16)

    pl2 = subplot(2,2,3)
    pcolor(lon_buoy,lat_buoy,squeeze(buoyancy(:,:,250))); shading flat; colorbar; hold on
    axis([151 159 -58 -52])
    plot(lon_buoy(choose_sec,:),lat_buoy(choose_sec,:),'k-','LineWidth',2)
    ccc = cmocean('dens');
    set(gca,'YTick',[-57:-53],'YTickLabel',['57';'56';'55';'54';'53'])
    set(gca,'Colormap',ccc)
    %caxis([27 max(siggy)]);
    set(gca,'FontSize',16)
    xlabel('Longitude (°E)','FontSize',16)
    ylabel('Latitude (°S)','FontSize',16)
    
    pl3 = subplot(3,2,2)
    holderdens = squeeze(mockdensity_3(choose_sec,:,:));
    pcolor(lat_buoy(choose_sec,:),depthints,squeeze(buoyancy(choose_sec,:,:))'); shading flat; colorbar; axis ij; hold on
    ccc = cmocean('dens');
    set(gca,'Colormap',ccc)   
    %axis([-57.7 -52.7 0 950])
    set(gca,'FontSize',16)
    %xlabel('Latitude (°S)','FontSize',16)
    ylabel('Depth (m)','FontSize',16)
    
    latposition = lat(42,1:457);
         
    pl4 = subplot(2,2,4)
    pcolor(lat_buoy(choose_sec,:),depthints,squeeze(M2_along_sm(choose_sec,:,:))'); shading flat; colorbar; 
    axis ij
    ccc = cmocean('balance');
    set(gca,'colormap',ccc)
    %axis([-57.7 -52.5 0 1000])
    caxis([-1.5e-7 1.5e-7])
    set(gca,'FontSize',16)
    xlabel('Longitude (°E)','FontSize',16)
    ylabel('Latitude (°S)','FontSize',16)
    title('M2 along (b_y) (s^-^2)','FontSize',16) 
    
    set(pl1,'Position',[0.06 0.56 0.25 0.42])
    set(pl2,'Position',[0.06 0.06 0.25 0.42])
    set(pl3,'Position',[0.44 0.55 0.5 0.25])
    set(pl4,'Position',[0.44 0.08 0.5 0.25])
    %set(pl5,'Position',[0.44 0.06 0.5 0.25])
    %pause
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Calculate the vertical velocity shear using the cyclostrophic surface
% currents and geostrophic shear

ang = NaN*ones(size(lon_buoy));
cross_u = NaN*ones(size(lon)); cross_u_fit13 = NaN*ones(size(lon)); cross_uc_fit13 = NaN*ones(size(lon));
along_v = NaN*ones(size(lon)); along_v_fit13 = NaN*ones(size(lon)); along_vc_fit13 = NaN*ones(size(lon));
cross_uc_fit13_buoy = NaN*ones(size(lon_buoy)); along_vc_fit13_buoy = NaN*ones(size(lon_buoy));

% First calculate the along and across velocity fields
'Calculate along and across velocity fields'
holderu = ucgos_fit13(:,:,choose_time); 
holderv = vcgos_fit13(:,:,choose_time);
holderu_buoy = holderu(:,latrange); holderv_buoy = holderv(:,latrange);
    
for i=2:size(lon_buoy,1)-1
    for j=2:size(lon_buoy,2)-1
        
            [dist(i,j), ang(i,j)] = sw_dist([lat_buoy(i-1,j) lat_buoy(i+1,j)],[lon_buoy(i-1,j) lon_buoy(i+1,j)],'km');
        
%           cross_u(i,j) = ugeo(i,j,9)*cos(ang(i,j)/180*pi) + vgeo(i,j,9)*sin(ang(i,j)/180*pi);
%           along_v(i,j) = ugeo(i,j,9)*sin(ang(i,j)/180*pi) + vgeo(i,j,9)*cos(ang(i,j)/180*pi);
%        
%           cross_u_fit13(i,j) = ugeo_fit13(i,j,9)*cos(ang(i,j)/180*pi) + vgeo_fit13(i,j,9)*sin(ang(i,j)/180*pi);
%           along_v_fit13(i,j) = ugeo_fit13(i,j,9)*sin(ang(i,j)/180*pi) + vgeo_fit13(i,j,9)*cos(ang(i,j)/180*pi);
%         
            cross_uc_fit13_buoy(i,j) =  holderu_buoy(i,j)*cos(ang(i,j)/180*pi) + holderv_buoy(i,j)*sin(ang(i,j)/180*pi);
            along_vc_fit13_buoy(i,j) = -holderu_buoy(i,j)*sin(ang(i,j)/180*pi) + holderv_buoy(i,j)*cos(ang(i,j)/180*pi);
        
    end
end

% Calcuate the 2D Coriolis field
ff = sw_f(lat_buoy);
ffmat = NaN*ones(size(buoyancy)); for i=1:490; ffmat(:,:,i) = ff; end

% Use 2-km along and across distances
acrossdist = 2000;
alongdist = 2000;

% Calculate velocity shear structure
'Calcuate the velocity shear structure'
shear_along_sm = -cumsum(M2_across_sm,3)*2./ffmat;
shear_across_sm = cumsum(M2_along_sm,3)*2./ffmat;


% Find the total velocity
'Calculate the total velocity'
full_across_sm = NaN*ones(size(buoyancy)); full_along_sm = NaN*ones(size(buoyancy));
for i=1:490
    full_across_sm(:,:,i) = cross_uc_fit13_buoy + shear_across_sm(:,:,i);
    full_along_sm(:,:,i)  = along_vc_fit13_buoy + shear_along_sm(:,:,i);
end


% Calculate the velocity gradients
'Calculate the velocity gradients'
full_across_sm_x = NaN*ones(size(buoyancy)); full_across_sm_y = NaN*ones(size(buoyancy));
full_along_sm_x  = NaN*ones(size(buoyancy)); full_along_sm_y  = NaN*ones(size(buoyancy));

    full_across_sm_x(2:end-1,:,:) = (full_across_sm(3:end,:,:) - full_across_sm(1:end-2,:,:))/2/2000;
    full_across_sm_y(:,2:end-1,:) = (full_across_sm(:,3:end,:) - full_across_sm(:,1:end-2,:))/2/2000;
    full_along_sm_x(2:end-1,:,:)  = (full_along_sm(3:end,:,:)  - full_along_sm(1:end-2,:,:))/2/2000;
    full_along_sm_y(:,2:end-1,:)  = (full_along_sm(:,3:end,:)  - full_along_sm(:,1:end-2,:))/2/2000;


% Add some smoothing to the velocity gradients
'Smooth the velocity gradients'
full_across_sm_x_sm = NaN*ones(size(buoyancy)); full_across_sm_y_sm = NaN*ones(size(buoyancy));
full_along_sm_x_sm  = NaN*ones(size(buoyancy)); full_along_sm_y_sm  = NaN*ones(size(buoyancy));
for i=1:490
    full_along_sm_x_sm(:,:,i) = smoothdata(full_along_sm_x(:,:,i),1,'gaussian',gauss_smooth);
    full_along_sm_x_sm(:,:,i) = smoothdata(full_along_sm_x_sm(:,:,i),2,'gaussian',gauss_smooth);
    
    full_along_sm_y_sm(:,:,i) = smoothdata(full_along_sm_y(:,:,i),1,'gaussian',gauss_smooth);
    full_along_sm_y_sm(:,:,i) = smoothdata(full_along_sm_y_sm(:,:,i),2,'gaussian',gauss_smooth);
    
    full_across_sm_x_sm(:,:,i) = smoothdata(full_across_sm_x(:,:,i),1,'gaussian',gauss_smooth);
    full_across_sm_x_sm(:,:,i) = smoothdata(full_across_sm_x_sm(:,:,i),2,'gaussian',gauss_smooth);
    
    full_across_sm_y_sm(:,:,i) = smoothdata(full_across_sm_y(:,:,i),1,'gaussian',gauss_smooth);
    full_across_sm_y_sm(:,:,i) = smoothdata(full_across_sm_y_sm(:,:,i),2,'gaussian',gauss_smooth);
    
end

% Save some additional diagnostics --> vorticity, strain, divergence
'Calculate vorticity, strain and divergence fields'
vorticity = full_along_sm_x_sm - full_across_sm_y_sm;
strain_n = full_across_sm_x_sm - full_along_sm_y_sm;
strain_s = full_along_sm_x_sm + full_across_sm_y_sm;
strain = sqrt(strain_n.^2 + strain_s.^2);
divergence = full_across_sm_x_sm + full_along_sm_y_sm;


% Calculate the Q-vector terms
'Calculate the Q-vector terms'
%%%%%%%%%%%%%%%
Q1 = -full_across_sm_x_sm.*M2_across_sm - full_along_sm_x_sm.*M2_along_sm;
Q2 = -full_across_sm_y_sm.*M2_across_sm - full_along_sm_y_sm.*M2_along_sm;
%%%%%%%%%%%%%%%


% Calculate gradients in Q-vector terms
'Calculate the gradient of the Q-vector terms'
Q1_sm_x = NaN*ones(size(buoyancy)); Q2_sm_y = NaN*ones(size(buoyancy));

    Q1_sm_x(2:end-1,:,:) = (Q1(3:end,:,:)-Q1(1:end-2,:,:))/2/2000;
    Q2_sm_y(:,2:end-1,:) = (Q2(:,3:end,:)-Q2(:,1:end-2,:))/2/2000;

    
% Smooth the gradient-Q terms
'Do some smoothing on the Q gradient terms'
Q1_sm_x_sm = NaN*ones(size(buoyancy)); Q2_sm_y_sm = NaN*ones(size(buoyancy));
for i=1:490
    Q1_sm_x_sm(:,:,i) = smoothdata(Q1_sm_x(:,:,i),1,'gaussian',gauss_smooth);
    Q1_sm_x_sm(:,:,i) = smoothdata(Q1_sm_x_sm(:,:,i),2,'gaussian',gauss_smooth);
    
    Q2_sm_y_sm(:,:,i) = smoothdata(Q2_sm_y(:,:,i),1,'gaussian',gauss_smooth);
    Q2_sm_y_sm(:,:,i) = smoothdata(Q2_sm_y_sm(:,:,i),2,'gaussian',gauss_smooth);
    
end

% Define the Q-forcing field
'Define the Q-forcing field'
QF = 2*(Q1_sm_x + Q2_sm_y);
QFsmooth = 2*(Q1_sm_x_sm + Q2_sm_y_sm);
kk = find(isnan(buoyancy(:,200,200)));
QFsmooth(kk,:,:) = NaN;

% Tile the divergence of Q -- left and right swaths

'Indices:  '
lon_buoy(:,1)'
QFsmooth(:,1)'
if orientation==1
    LeftSwath_LeftEdge = 5;
    LeftSwath_RightEdge = 31;
    RightSwath_LeftEdge = 40;
    RightSwath_RightEdge = 66;
elseif orientation==2
    LeftSwath_LeftEdge = 5;
    LeftSwath_RightEdge = 30;
    RightSwath_LeftEdge = 40;
    RightSwath_RightEdge = 65;
end
    
%LeftSwath_LeftEdge = input('Enter the index of the left edge of the left swath:  ');
%LeftSwath_RightEdge = input('Enter the index of the right edge of the left swath:  ');
%RightSwath_LeftEdge = input('Enter the index of the left edge of the right swath:  ');
%RightSwath_RightEdge = input('Enter the index of the right edge of the right swath:  ');

tilesize_lon_left = (LeftSwath_RightEdge-LeftSwath_LeftEdge)*2;
tilesize_lon_right = (RightSwath_RightEdge-RightSwath_LeftEdge)*2;
tilesize_lat = 2*nn-2;


% Tile the forcing field
'Tile the forcing field'
QFsm_tile_left  = NaN*ones(tilesize_lon_left,tilesize_lat,oo);
QFsm_tile_right = NaN*ones(tilesize_lon_right,tilesize_lat,oo);

    QFsm_tile_left(1:(tilesize_lon_left)/2+1,1:(tilesize_lat)/2+1,:) = QFsmooth(LeftSwath_LeftEdge:LeftSwath_RightEdge,:,:); 
    QFsm_tile_left((tilesize_lon_left)/2+2:end,1:(tilesize_lat)/2+1,:) = QFsmooth(LeftSwath_RightEdge-1:-1:LeftSwath_LeftEdge+1,:,:);
    QFsm_tile_left(:,(tilesize_lat)/2+2:tilesize_lat,:) = QFsm_tile_left(:,tilesize_lat/2:-1:2,:);
    
    QFsm_tile_right(1:(tilesize_lon_right)/2+1,1:(tilesize_lat)/2+1,:) = QFsmooth(RightSwath_LeftEdge:RightSwath_RightEdge,:,:); 
    QFsm_tile_right((tilesize_lon_right)/2+2:end,1:(tilesize_lat)/2+1,:) = QFsmooth(RightSwath_RightEdge-1:-1:RightSwath_LeftEdge+1,:,:);
    QFsm_tile_right(:,(tilesize_lat)/2+2:tilesize_lat,:) = QFsm_tile_right(:,tilesize_lat/2:-1:2,:);
    
    % Clean up a few things!
    if SWOT_swath_no == '004' & date_index == 7;
        QFsm_tile_left(1,:,:) = (QFsm_tile_left(2,:,:)+QFsm_tile_left(end,:,:))/2;
    end
    if SWOT_swath_no == '004' & date_index == 8;
        QFsm_tile_left(1,:,:) = (QFsm_tile_left(2,:,:)+QFsm_tile_left(end,:,:))/2;
        QFsm_tile_right(1,:,:) = (QFsm_tile_right(2,:,:)+QFsm_tile_right(end,:,:))/2;
    end
    if SWOT_swath_no == '032' & date_index == 7;
        QFsm_tile_left(2,81,:) = (QFsm_tile_left(2,80,:)+QFsm_tile_left(2,82,:))/2;
        QFsm_tile_left(2,635,:) = (QFsm_tile_left(2,634,:)+QFsm_tile_left(2,636,:))/2;
        QFsm_tile_left(1,:,:) = (QFsm_tile_left(2,:,:)+QFsm_tile_left(end,:,:))/2;
    end
    if SWOT_swath_no == '032' & date_index == 8;
        QFsm_tile_left(1,:,:) = (QFsm_tile_left(2,:,:)+QFsm_tile_left(end,:,:))/2;
    end
    if SWOT_swath_no == '060' & date_index == 8;
        QFsm_tile_left(1,:,:) = (QFsm_tile_left(2,:,:)+QFsm_tile_left(end,:,:))/2;
    end
    if SWOT_swath_no == '099' & date_index == 8;
        QFsm_tile_left(1,:,:) = (QFsm_tile_left(2,:,:)+QFsm_tile_left(end,:,:))/2;
    end
    if SWOT_swath_no == '155' & date_index == 8;
        QFsm_tile_left(1,:,:) = (QFsm_tile_left(2,:,:)+QFsm_tile_left(end,:,:))/2;
        QFsm_tile_right(1,:,:) = (QFsm_tile_right(2,:,:)+QFsm_tile_right(end,:,:))/2;
    end
    if SWOT_swath_no == '155' & date_index == 9;
        QFsm_tile_left(1,:,:) = (QFsm_tile_left(2,:,:)+QFsm_tile_left(end,:,:))/2;
        QFsm_tile_right(1,:,:) = (QFsm_tile_right(2,:,:)+QFsm_tile_right(end,:,:))/2;
    end
     if SWOT_swath_no == '310' & date_index == 8;
        QFsm_tile_left(1,:,:) = (QFsm_tile_left(2,:,:)+QFsm_tile_left(end,:,:))/2;
        QFsm_tile_right(1,:,:) = (QFsm_tile_right(2,:,:)+QFsm_tile_right(end,:,:))/2;
     end
    if SWOT_swath_no == '377' & date_index == 8;
        QFsm_tile_right(1,:,:) = (QFsm_tile_right(2,:,:)+QFsm_tile_right(end,:,:))/2;
    end
    if SWOT_swath_no == '433' & date_index == 9;
        QFsm_tile_right(1,:,:) = (QFsm_tile_right(2,:,:)+QFsm_tile_right(end,:,:))/2;
    end
    if SWOT_swath_no == '560' & date_index == 8;
        QFsm_tile_left(1,:,:) = (QFsm_tile_left(2,:,:)+QFsm_tile_left(end,:,:))/2;
    end
    if SWOT_swath_no == '560' & date_index == 9;
        QFsm_tile_left(1,:,:) = (QFsm_tile_left(2,:,:)+QFsm_tile_left(end,:,:))/2;
    end
    
    figure(66); close
    hh = figure(66); set(hh,'Position',[63          83        1295         766])
    pl1 = subplot(2,1,1)
    maxie = max(abs(QFsm_tile_left(:)));
    imagesc(QFsm_tile_left(:,:,200)); caxis([-maxie maxie]); 
    
    pl2 = subplot(2,1,2)
    maxie = max(abs(QFsm_tile_right(:)));
    imagesc(QFsm_tile_right(:,:,200)); caxis([-maxie maxie]); 
    
    %pause
         

% for i=1:oo
%     for j=1:nn
%         kk = find(isnan(QFsm_tile_left(:,j,i)));
%         if ~isempty(kk) 
%             if kk==1
%                 QFsm_tile_left(kk,j,i) = (QFsm_tile_left(kk,j+1,i)+QFsm_tile_left(kk,j-1,i))/2;
%             else
%                 QFsm_tile_left(kk,j,i) = (QFsm_tile_left(kk-1,j,i)+QFsm_tile_left(kk+1,j,i))/2;
%             end
%         end
%         kk = find(isnan(QFsm_tile_right(:,j,i)));
%         if ~isempty(kk) 
%             if kk==1
%                 QFsm_tile_right(kk,j,i) = (QFsm_tile_right(kk,j+1,i)+QFsm_tile_right(kk,j-1,i))/2;
%             else
%                 QFsm_tile_right(kk,j,i) = (QFsm_tile_right(kk,j+1,i)+QFsm_tile_right(kk,j-1,i))/2;
%             end
%         end
%     end
% end

% Extend to deeper depths
'Extend to deeper depths'
pp = oo+oo/2;
QFsm_extend_left = NaN*ones(tilesize_lon_left,tilesize_lat,pp);
QFsm_extend_right = NaN*ones(tilesize_lon_right,tilesize_lat,pp);

QFsm_extend_left(:,:,1:oo) = QFsm_tile_left;
QFsm_extend_right(:,:,1:oo) = QFsm_tile_right;
for i=oo+1:pp; QFsm_extend_left(:,:,i) = QFsm_tile_left(:,:,oo); QFsm_extend_right(:,:,i) = QFsm_tile_right(:,:,oo); end

N2mean_extend = NaN*ones(1,pp); N2mean_extend(1:oo) = N2mean; 
kk = find(depthints>900 & depthints<950); N2mean_extend(oo+1:end) = nanmean(N2mean(kk));
depth_extend = NaN*ones(1,pp);
depth_extend(1:oo) = depthints; depth_extend(oo+1:pp) = 991:2:1479;

% Truncate the vertical direction for computational efficiency
QFsm_trunc_left = NaN*ones(tilesize_lon_left,tilesize_lat,pp/3);
QFsm_trunc_right = NaN*ones(tilesize_lon_right,tilesize_lat,pp/3);
for i=1:pp/3
    QFsm_trunc_left(:,:,i) = nanmean(QFsm_extend_left(:,:,(i*3)-2:(i*3)),3);
    QFsm_trunc_right(:,:,i) = nanmean(QFsm_extend_right(:,:,(i*3)-2:(i*3)),3);
end
depth_trunc = depth_extend(2:3:end-1);
for i=1:pp/3; N2mean_trunc(i) = nanmean(N2mean_extend( (i*3)-2:i*3 )); end


% Write to nc file
'Write to nc file'
nc_filename = ['FOCUS_data_trunc_swath_' SWOT_swath_no '_' num2str(date_index) '_GS8_SSHDens2.nc'];
nccreate(nc_filename,'QF_trunc_left','Dimensions',{'crosstile',tilesize_lon_left,'alongtile',tilesize_lat,'depth_trunc_size',pp/3});
nccreate(nc_filename,'QF_trunc_right','Dimensions',{'crosstile',tilesize_lon_right,'alongtile',tilesize_lat,'depth_trunc_size',pp/3});
nccreate(nc_filename,'N2mean_trunc','Dimensions',{'depth_trunc_size',245});
nccreate(nc_filename,'vorticity','Dimensions',{'crossfull',size(lon_buoy,1),'alongfull',size(lon_buoy,2),'depthfull',length(depthints)});
nccreate(nc_filename,'strain','Dimensions',{'crossfull',size(lon_buoy,1),'alongfull',size(lon_buoy,2),'depthfull',length(depthints)});
nccreate(nc_filename,'divergence','Dimensions',{'crossfull',size(lon_buoy,1),'alongfull',size(lon_buoy,2),'depthfull',length(depthints)});
nccreate(nc_filename,'buoyancy','Dimensions',{'crossfull',size(lon_buoy,1),'alongfull',size(lon_buoy,2),'depthfull',length(depthints)});
nccreate(nc_filename,'w_vel','Dimensions',{'crossfull',size(lon_buoy,1),'alongfull',size(lon_buoy,2),'depth_trunc_size',pp/3});
nccreate(nc_filename,'lon_buoy','Dimensions',{'crossfull',size(lon_buoy,1),'alongfull',size(lon_buoy,2)});
nccreate(nc_filename,'lat_buoy','Dimensions',{'crossfull',size(lat_buoy,1),'alongfull',size(lat_buoy,2)});
nccreate(nc_filename,'depth_trunc','Dimensions',{'depth_trunc_size',245});
nccreate(nc_filename,'LeftSwath_LeftIndex','Dimensions',{'swathindex',1});
nccreate(nc_filename,'LeftSwath_RightIndex','Dimensions',{'swathindex',1});
nccreate(nc_filename,'RightSwath_LeftIndex','Dimensions',{'swathindex',1});
nccreate(nc_filename,'RightSwath_RightIndex','Dimensions',{'swathindex',1});
nccreate(nc_filename,'YearDay','Dimensions',{'yearday',1});
nccreate(nc_filename,'orientation','Dimensions',{'swathindex',1})

ncwrite(nc_filename,'QF_trunc_left',QFsm_trunc_left);
ncwrite(nc_filename,'QF_trunc_right',QFsm_trunc_right);
ncwrite(nc_filename,'N2mean_trunc',N2mean_trunc);
ncwrite(nc_filename,'vorticity',vorticity);
ncwrite(nc_filename,'strain',strain);
ncwrite(nc_filename,'divergence',divergence);
ncwrite(nc_filename,'buoyancy',buoyancy);
ncwrite(nc_filename,'lon_buoy',lon_buoy);
ncwrite(nc_filename,'lat_buoy',lat_buoy);
ncwrite(nc_filename,'depth_trunc',depth_trunc);
ncwrite(nc_filename,'LeftSwath_LeftIndex',LeftSwath_LeftEdge);
ncwrite(nc_filename,'LeftSwath_RightIndex',LeftSwath_RightEdge);
ncwrite(nc_filename,'RightSwath_LeftIndex',RightSwath_LeftEdge);
ncwrite(nc_filename,'RightSwath_RightIndex',RightSwath_RightEdge);
ncwrite(nc_filename,'YearDay',time(1,date_index));
ncwrite(nc_filename,'orientation',orientation);

%%%%%%%%%%%%%%%
'Add the vertical velocities to the .nc file'
 
wleft = ncread(['w_' num2str(SWOT_swath_no) '_' num2str(date_index) '_GS8_DensSSH2.nc'],'w_left');
wright = ncread(['w_' num2str(SWOT_swath_no) '_' num2str(date_index) '_GS8_DensSSH2.nc'],'w_right');

wfull = NaN*ones(size(buoyancy,1),size(buoyancy,2),length(depth_trunc));
if orientation==1
    wfull(5:31,:,:) = wleft(1:27,1:nn,:);
    wfull(40:66,:,:) = wright(1:27,1:nn,:);
elseif orientation==2
    wfull(5:30,:,:) = wleft(1:26,1:nn,:);
    wfull(40:65,:,:) = wright(1:26,1:nn,:);
end

ncwrite(nc_filename,'w_vel',wfull);


% w_right = ncread('merge_w_032_8.nc','w_right');
% 
% kk = find(depth_trunc>950,1);
% w_swath = NaN*ones(size(lon_buoy,1),size(lat_buoy,2),length(depth_trunc(1:kk-1)));
% 
% for i=1:kk-1
%     w_swath(3:34,1:370,i) = w_left(1:32,1:370,i);
%     w_swath(37:68,1:370,i) = w_right(1:32,1:370,i);
%     i
% end

cd /Users/andy/Work/Manuscripts/Andy_FOCUS/
