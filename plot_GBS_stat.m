% plot measurement statistics

%% old code, use plot_GBS_meas_stat.m

%%

% UT-GBS ZS measurements
load('/home/kristof/work/GBS/UT-GBS/2016/VCD/ozone/ozone_v2_86-90.mat')

day_sza_filt=floor(qdoas_raw(find(qdoas_raw(:,7) <=91),3));

% figure(222)
% [a,b]=unique(day_sza_filt);
% b=[diff(b)', size(day_sza_filt,1)-b(end)+1];
% a=a+0.25;
% 
% bar(a,b,0.5,'c'), hold on
% b=b./2
% a=a+0.5
% bar(a,b,0.5,'b'), hold on

height=710;
width=1158;

figure('Position', [100, 100, width, height]);
subplot(211)
histogram(day_sza_filt); hold on

day_sza_filt=floor(qdoas_raw(find(qdoas_raw(:,7) <=91 & qdoas_raw(:,7) >=86),3));
histogram(day_sza_filt)

legend('SZA<86','SZA=[86,91]')
ylabel('N.o. measurements')
title('Number of Zenith Sky measurements')

subplot(212)
% UT-GBS does MAX-DOAS for SZA<84 (<86 for PGBS)

day_sza_filt=floor(qdoas_raw(find(qdoas_raw(:,7) <=84),3));
histogram(day_sza_filt)
xlabel('Day of the year (UTC)')
ylabel('N.o. measurements')
title('Number of MAX-DOAS scans (5 elevations)')
xlim([50,95])

set(findall(gcf,'-property','FontSize'),'FontSize',15.6)
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman') 


% PEARL-GBS ZS measurements
load('/home/kristof/work/GBS/PEARL-GBS/2016/VCD/no2_oneDC/no2_v1_86-91.mat')
day_sza_filt=floor(qdoas_raw(find(qdoas_raw(:,7) <=91),3));

% figure(2)
figure('Position', [100, 100, width, height]);
subplot(211)
histogram(day_sza_filt), hold on

day_sza_filt=floor(qdoas_raw(find(qdoas_raw(:,7) <=91 & qdoas_raw(:,7) >=86),3));
histogram(day_sza_filt)

legend('SZA<86','SZA=[86,91]')
ylabel('N.o. measurements')
title('Number of Zenith Sky measurements')

load('/home/kristof/work/GBS/PEARL-GBS/2016/MAX-DOAS/dscd15.mat')

subplot(212)
h=histogram(dscd_S_15.day);
xlabel('Day of the year (UTC)')
ylabel('N.o. measurements')
title('Number of MAX-DOAS scans (8 elevations)')
xlim([60,95])

set(findall(gcf,'-property','FontSize'),'FontSize',15.6)
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman') 

clear
