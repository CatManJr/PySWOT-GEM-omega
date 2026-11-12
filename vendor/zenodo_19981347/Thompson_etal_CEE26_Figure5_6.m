% This script generates the figure showing statistical properties of
% vertical velocity
%
% AFT, Last updated: May 2, 2026
cd /Users/andy/Work/SWOT/Benoit/GliderData/

yn = input('Do you need to load the data (y==1, n==0)?  ')

if yn==1

% Select the SWOT swaths with clean vertical velocity distributions  
files = ['004_7';'004_8';'060_7';'060_8';'099_7';'099_8';'127_7';'155_8';'155_9';'183_9';'183_1';'282_8';'338_8';'377_8';'405_9';'433_9';'560_8';'560_9']

% Pameters related to the swath geometry
mm = 69;
nn = 358;
oo = 490;
qq = 245;
no_files = size(files,1);

% Set the Coriolis frequency
ff = sw_f(-55);


w_holder = NaN*ones(no_files*mm,nn,qq);
vort_holder = NaN*ones(no_files*mm,nn,oo);
strain_holder = NaN*ones(no_files*mm,nn,oo);

'Load data'
for i=1:no_files
    holderfile = files(i,:);
    if holderfile == '183_1'; holderfile = '183_10'; end
    
    w_holder(mm*(i-1)+1:mm*i,:,:)      = ncread(['FOCUS_data_trunc_swath_' holderfile '.nc'],'w_vel');
    vort_holder(mm*(i-1)+1:mm*i,:,:)   = ncread(['FOCUS_data_trunc_swath_' holderfile '.nc'],'vorticity');
    strain_holder(mm*(i-1)+1:mm*i,:,:) = ncread(['FOCUS_data_trunc_swath_' holderfile '.nc'],'strain');
    i

end

w_holder_950 = w_holder(:,:,1:157);
vort_holder_950 = vort_holder(:,:,2:3:470);
strain_holder_950 = strain_holder(:,:,2:3:470);


histedges_w = -600:4:600;
histmid_w = (histedges_w(1:end-1)+histedges_w(2:end))/2;
ll_w = length(histedges_w);

%histedges_vort = -1:.01:1;
histedges_vort = -1:.02:1;
histmid_vort = (histedges_vort(1:end-1)+histedges_vort(2:end))/2;
ll_vort = length(histedges_vort);

%histedges_strain = -1:.005:0;
histedges_strain = -1:.01:0;
histmid_strain = (histedges_strain(1:end-1)+histedges_strain(2:end))/2;
ll_strain = length(histedges_strain);

nn_wvel = NaN*ones(no_files,ll_w-1); mean_w = NaN*ones(1,no_files); skew_w = NaN*ones(1,no_files);
nn_vort = NaN*ones(no_files,ll_vort-1); mean_vort = NaN*ones(1,no_files); skew_vort = NaN*ones(1,no_files);
nn_strain = NaN*ones(no_files,ll_strain-1); mean_strain = NaN*ones(1,no_files); %skew_vort = NaN*ones(1,no_files);

for i=1:no_files
    holder_w = w_holder_950((i-1)*mm+1:i*mm,:,:);
    holder_vort = vort_holder_950((i-1)*mm+1:i*mm,:,:);
    holder_strain = strain_holder_950((i-1)*mm+1:i*mm,:,:);
    
    mean_w(i) = nanmean(holder_w(:)); mean_vort(i) = nanmean(holder_vort(:)); mean_strain(i) = nanmean(holder_strain(:));
    skew_w(i) = skewness(holder_w(:)); skew_vort(i) = skewness(holder_vort(:));
    nn_wvel(i,:) = hist(holder_w(:)*60*60*24,histmid_w);
    nn_vort(i,:) = hist(holder_vort(:)/ff,histmid_vort);
    nn_strain(i,:) = hist(holder_strain(:)/ff,histmid_strain);
    i
end
nn_wvel_tot = hist(w_holder_950(:)*60*60*24,histmid_w);
mean_w_tot = nanmean(w_holder_950(:)); skew_w_tot = skewness(w_holder_950(:));

nn_vort_tot = hist(vort_holder_950(:)/ff,histmid_vort);
%nn_vort_tot = histc(vort_holder_950(:)/ff,histedges_vort;
mean_vort_tot = nanmean(vort_holder_950(:)); skew_vort_tot = skewness(vort_holder_950(:));

nn_strain_tot = hist(strain_holder_950(:)/ff,histmid_strain);
mean_strain_tot = nanmean(strain_holder_950(:));

% Calculate a conditional mean of vertical velocity and vorticity
for i=1:length(histmid_vort);
    kk = find(vort_holder_950(:)/ff>histedges_vort(i) & vort_holder_950(:)/ff<histedges_vort(i+1));
    mean_w_cond_vort(i) = nanmean(w_holder_950(kk));
    std_w_cond_vort(i) = nanstd(w_holder_950(kk));
    i
end

% Calculate a conditional mean of vertical velocity on both strain and
% vorticty
mean_w_cond_vort_strain = NaN*ones(length(histmid_vort),length(histmid_strain));
median_w_cond_vort_strain = NaN*ones(length(histmid_vort),length(histmid_strain));
std_w_cond_vort_strain = NaN*ones(length(histmid_vort),length(histmid_strain));
nn_w_cond_vort_strain = NaN*ones(length(histmid_vort),length(histmid_strain));

for i=1:length(histedges_vort)-1;
    for j=1:length(histedges_strain)-1;
        kk = find(vort_holder_950(:)/ff>histedges_vort(i) & vort_holder_950(:)/ff<histedges_vort(i+1) & strain_holder_950(:)/ff>histedges_strain(j) & strain_holder_950(:)/ff<histedges_strain(j+1));
        mean_w_cond_vort_strain(i,j) = nanmean(w_holder_950(kk));
        median_w_cond_vort_strain(i,j) = nanmedian(w_holder_950(kk));
        std_w_cond_vort_strain(i,j) = nanstd(w_holder_950(kk));
        nn_w_cond_vort_strain(i,j) = length(kk);
    end
    i
end

% Calculate a conditional mean of vertical velocity on both strain and
% vorticty -- only between 100 and 400 m
mean_w_cond_vort_strain_100_400 = NaN*ones(length(histmid_vort),length(histmid_strain));
std_w_cond_vort_strain_100_400 = NaN*ones(length(histmid_vort),length(histmid_strain));
nn_w_cond_vort_strain_100_400 = NaN*ones(length(histmid_vort),length(histmid_strain));

for i=1:length(histedges_vort)-1;
    for j=1:length(histedges_strain)-1;
        holdervort = vort_holder_950(:,:,16:66);
        holderstrain = strain_holder_950(:,:,16:66);
        holderw = w_holder_950(:,:,16:66);
        kk = find(holdervort(:)/ff>histedges_vort(i) & holdervort(:)/ff<histedges_vort(i+1) & holderstrain(:)/ff>histedges_strain(j) & holderstrain(:)/ff<histedges_strain(j+1));
        mean_w_cond_vort_strain_100_400(i,j) = nanmean(holderw(kk));
        std_w_cond_vort_strain_100_400(i,j) = nanstd(holderw(kk));
        nn_w_cond_vort_strain_100_400(i,j) = length(kk);
    end
    i
end


% Calculate the standard deviation as a function of depth
for i=1:157
    
    holder_w = w_holder_950(:,:,i); holder_vort = vort_holder_950(:,:,i);
    w_std_depth(i) = nanstd(holder_w(:));
    vort_std_depth(i) = nanstd(holder_w(:));

end

end

%save w_data histmid_w histmid_vort histmid_strain nn_wvel nn_wvel_tot nn_vort nn_vort_tot w_std_depth nn_w_cond_vort_strain std_w_cond_vort_strain mean_w_cond_vort_strain skew_w skew_vort_tot skew_vort
load w_data

figure(88); close
hh = figure(88); set(hh,'Position',[ 60         118        1222         748])

pl1 = subplot(2,3,1)

for i=1:no_files
    summer = sum(nn_wvel(i,2:end-1))*4;
%    plot(histmid(2:end-1),nn_wvel(i,2:end-1)/summer,'Color',[.5 .5 .5],'LineWidth',.5); hold on
    plot(histmid_w(2:end-1),nn_wvel(i,2:end-1)/summer,'LineWidth',.5); hold on
    set(gca,'YScale','log')
    i
    %pause
end
summer = sum(nn_wvel_tot)*4;
plot(histmid_w(2:end-1),nn_wvel_tot(2:end-1)/summer,'-','LineWidth',2,'Color',[0 84 147]/256); hold on
%plot(histmid((ll-1)/2+1:end-1),nn(ll/2:-1:2),'r-','LineWidth',1)
set(gca,'FontSize',14);
xlabel('Vertical velocity (m day^-^1)','FontSize',16)
ylabel('Probability distribution function (day m^-^1)','FontSize',16)
t1 = text(-583.6289,    0.0586,'(a)','FontSize',14)

pl2 = subplot(2,3,2)

for i=1:no_files
    summer = sum(nn_vort(i,2:end-1))*.01;
%    plot(histmid(2:end-1),nn_wvel(i,2:end-1)/summer,'Color',[.5 .5 .5],'LineWidth',.5); hold on
    plot(histmid_vort(2:end-1),nn_vort(i,2:end-1)/summer,'LineWidth',.5); hold on
end
summer = sum(nn_vort_tot)*.01;
plot(histmid_vort(2:end-1),nn_vort_tot(2:end-1)/summer,'-','LineWidth',2,'Color',[148 23 81]/256); hold on
%plot(histmid((ll-1)/2+1:end-1),nn(ll/2:-1:2),'r-','LineWidth',1)
set(gca,'YScale','log')
set(gca,'FontSize',14);
xlabel('Vertical vorticity (\zeta/f)','FontSize',16)
ylabel('Probability distribution function','FontSize',16)
t2 = text(-0.9673,    5.5898,'(b)','FontSize',14)

pl3 = subplot(1,3,3)
plot(w_std_depth*60*60*24,depth_trunc(1:157),'ko','MarkerSize',12,'MarkerFaceColor',[0 84 147]/256); hold on; axis ij
%plot(-vort_std_depth/ff,depth_trunc(1:157),'ko','MarkerSize',12,'MarkerFaceColor',[148 23 81]/256); hold on
set(gca,'FontSize',14);
xlabel('Standard deviation (w, m day^-^1)','FontSize',16)
ylabel('Depth (m)','FontSize',16)
t3 = text(54.7239,   26.7380,'(e)','FontSize',14)

% pl4 = subplot(2,3,4)
% plot(histmid_vort,(mean_w_cond_vort)*60*60*24,'k.-','LineWidth',2); hold on
% plot(-1:.1:1,zeros(1,21),'--','Color',[.5 .5 .5],'LineWidth',.5)
% axis([-1 1 -50 100])
% set(gca,'FontSize',14);
% xlabel('Vertical vorticity (\zeta/f)','FontSize',16)
% ylabel('Mean vertical velocity (m day^-^1)','FontSize',16)

pl8 = subplot(4,3,8)
plot(1:18,skew_w,'ko','MarkerSize',12,'Color','k','MarkerFaceColor',[0 84 147]/256); hold on
plot(0:19,skew_vort_tot*ones(1,20),'k--','LineWidth',2)
set(gca,'XTick',[2:2:18])
axis([0 19 -4 4])
set(gca,'XTick',[2:2:18])
set(gca,'FontSize',14);
xlabel('SWOT pass','FontSize',16)
ylabel('w skewness','FontSize',16)
t4 = text(0.3261,    3.1973,'(c)','FontSize',14)

pl11 = subplot(4,3,11)
plot(1:18,skew_vort,'ks','MarkerSize',12,'Color','k','MarkerFaceColor',[148 23 81]/256); hold on
plot(0:19,skew_vort_tot*ones(1,20),'k--','LineWidth',2)
axis([0 19 -1 1])
set(gca,'XTick',[2:2:18])
set(gca,'FontSize',14);
xlabel('SWOT pass','FontSize',16)
ylabel('\zeta/f skewness','FontSize',16)
t5 = text(0.3261,    0.7993,'(d)','FontSize',14)


set(pl1,'Position',[0.05 0.52 0.3 0.45])
set(pl2,'Position',[0.43 0.52 0.3 0.45])
%set(pl4,'Position',[0.05 0.06 0.37 0.35])
set(pl8,'Position',[0.05 0.22 0.30 0.20])
set(pl11,'Position',[0.43 0.22 0.30 0.20])
set(pl3,'Position',[0.79 0.22 0.2 0.75])


figure(89); close
hh = figure(89); set(hh,'Position',[ 60         118        1222         748])

kk = find(nn_w_cond_vort_strain<=50);
nn_w_cond_vort_strain(kk) = NaN; mean_w_cond_vort_strain(kk) = NaN;
std_w_cond_vort_strain(kk) = NaN;

pl1 = subplot(1,3,1)
pcolor(histmid_vort,histmid_strain,mean_w_cond_vort_strain'*60*60*24); shading flat; axis ij; hold on
colorbar
ccc = cmocean('balance');
set(gca,'FontSize',16,'Colormap',ccc)
caxis([-120 120])
axis([-1 1 -.85 0])
xx = -1:.05:0; yy = -1:.05:0; plot(xx,yy,'k--','LineWidth',1)
xx = 0:.05:1; yy = 0:-.05:-1; plot(xx,yy,'k--','LineWidth',1)
xlabel('Vertical vorticity (\zeta/f)','FontSize',18)
ylabel('Strain (\alpha/f)','FontSize',18)
t1 = text(-0.9745,   -0.8253,'(a)','FontSize',16)
t1a = text( 1.2836,   -0.1490  ,'Conditional mean vertical velocity (m day^-^1)','FontSize',18,'Rotation',90)

pl2 = subplot(1,3,2)
pcolor(histmid_vort,histmid_strain,log10(nn_w_cond_vort_strain')); shading flat; axis ij; hold on
colorbar
ccc = cmocean('thermal');
set(gca,'FontSize',16,'Colormap',ccc)
axis([-1 1 -.85 0])
xx = -1:.05:0; yy = -1:.05:0; plot(xx,yy,'k--','LineWidth',1)
xx = 0:.05:1; yy = 0:-.05:-1; plot(xx,yy,'k--','LineWidth',1)
%xlabel('Vertical vorticity (\zeta/f)','FontSize',18)
ylabel('Strain (\alpha/f)','FontSize',18)
t2 = text(-0.9614,   -0.8030,'(b)','FontSize',16)
t2a = text(1.4164,   -0.2982,'Log_1_0 count','FontSize',18,'Rotation',90)


pl3 = subplot(1,3,3)
pcolor(histmid_vort,histmid_strain,std_w_cond_vort_strain'*60*60*24); shading flat; axis ij; hold on
colorbar
caxis([0 200])
axis([-1 1 -.85 0])
xx = -1:.05:0; yy = -1:.05:0; plot(xx,yy,'k--','LineWidth',1)
xx = 0:.05:1; yy = 0:-.05:-1; plot(xx,yy,'k--','LineWidth',1)
ccc = cmocean('tempo');
set(gca,'FontSize',16,'Colormap',ccc)
set(gca,'FontSize',16)
xlabel('Vertical vorticity (\zeta/f)','FontSize',18)
ylabel('Strain (\alpha/f)','FontSize',18)
t3 = text(-0.9614,   -0.8030,'(c)','FontSize',16)
t3a = text( 1.4164,   -0.1909 ,'Std dev w (m day^-^1)','FontSize',18,'Rotation',90)


set(pl1,'Position',[.05 .19 .45 .72])
set(pl2,'Position',[.65 .57 .28 .38])
set(pl3,'Position',[.65 .15 .28 .38])


cd /Users/andy/Work/Manuscripts/Andy_FOCUS/
