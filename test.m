%% total avks
% aod=prctile(info_aer.col,0:10:100);
% for i=1:length(aod)-1
% 
%     ind=(info_aer.col>=aod(i) & info_aer.col < aod(i+1));
% 
%     plot(mean(avk_col(:,ind),2),alt), hold on
%     
% end    
plot(mean(avk_col,2), alt), hold on
tmp=max(mean(avk_col,2));
e=2.718281828459046;
plot([tmp/e,tmp/e],[0,4],'k-')
grid on
grid minor

%% compare Ramina's QDOAS output to mine
% ind1=(data.O3RMS<0.003 & data.SZA<92 & data.SZA>86);
% ind2=(data_mine.O3RMS<0.003 & data_mine.SZA<92 & data_mine.SZA>86);
% plot(data.Fractionalday(ind1),data.O3SlColo3(ind1),'ko'), hold on
% plot(data_mine.Fractionalday(ind2),data_mine.O3SlColo3(ind2),'rx'), hold on
% 
% [~, inda,indb]=intersect(data.Fractionalday,data_mine.Fractionalday);
% d1=data(inda,:);
% d2=data_mine(indb,:);
% ind=(d1.O3RMS<0.003 & d1.SZA<92 & d1.SZA>86 & d2.O3RMS<0.003 & d2.SZA<92 & d2.SZA>86);
% % ind=(d1.O3RMS<0.003 & d2.O3RMS<0.003);
% 
% figure
% % plot(d1.Fractionalday(ind),d1.O3SlColo3(ind),'ko'), hold on
% % plot(d2.Fractionalday(ind),d2.O3SlColo3(ind),'rx'), hold on
% plot(d1.Fractionalday(ind),((d1.O3SlColo3(ind)-d2.O3SlColo3(ind))./d1.O3SlColo3(ind))*100,'kx')

%% break up 2019-2020 berwer data (from Xiaoyi)
% for bw_num=unique(o3data_ews.serial)'
%     ind=o3data_ews.serial==bw_num;
%     tmp=table();
%     tmp.DateTime=o3data_ews.GMTDate(ind);
%     tmp.ColumnO3=o3data_ews.o3(ind);
%     tmp.StdDevO3=o3data_ews.err_o3(ind);
%     tmp.Airmass=o3data_ews.mu(ind);
%     tmp.ZA=o3data_ews.za(ind);
%     
%     brewer=tmp;
%     tmp=unique(tmp.DateTime.Year);
%     if length(tmp)==1
%         save(['brewer' num2str(bw_num) '_' num2str(tmp) '.mat'],'brewer')
%     else
%         save(['brewer' num2str(bw_num) '_' num2str(tmp(1)) '-' num2str(tmp(2)) '.mat'],'brewer')
%     end
% end
% ind=(brewer_69.Datetime.Year==2020 & brewer_69.Airmass<5 & brewer_69.StdDevO3<2.5);
% plot(brewer_69.Datetime(ind),brewer_69.ColumnO3(ind),'x'), hold on
% ind=(brewer_192.Datetime.Year==2020 & brewer_192.Airmass<5 & brewer_192.StdDevO3<2.5);
% plot(brewer_192.Datetime(ind),brewer_192.ColumnO3(ind),'x'), hold on
% ind=(brewer_223.Datetime.Year==2020 & brewer_223.Airmass<5 & brewer_223.StdDevO3<2.5);
% plot(brewer_223.Datetime(ind),brewer_223.ColumnO3(ind),'x'), hold on
%     
% legend('69','192','223')

%% bruker tg errors
% valid_tg={'O3','NO2','HCl','HNO3','ClONO2','HF'};
% for i=valid_tg
%     load(['/home/kristof/work/bruker/PEARL_ozone_depletion/bruker_' lower(i{1}) '.mat'])
%     err_sys=bruker.tot_col_err_sys./bruker.tot_col;
%     err_rand=bruker.tot_col_err_rand./bruker.tot_col;
%     err_tot=mean(sqrt(err_sys.^2 + err_rand.^2)*100);
%     disp([i{1} ': ' num2str(round(err_tot,1)) '%'])
% end
%% plot pandora RMS
% rms_lim=0.003;
% for ea=[1,2,3,5,8,10,15,20,30,40,50]
%     
%     ind=data.Elevviewingangle==ea;
%     figure
%     dscatter(data.SZA(ind),data.NO2_VisRMS(ind))
%     title(['elev: ' num2str(ea) ' deg'])
%     ylim([0,0.008])
%     grid on
%     
%     tmp=sum(data.NO2_VisRMS(ind)>rms_lim)/sum(ind);
%     disp(['Filtered ' num2str(tmp*100) '% at ea = ' num2str(ea)])
%     
% end

%% compare pandora dSCDs
% [fd,inda,indb]=intersect(data.Fractionalday,new.Fractionalday);
% 
% comp1=data.NO2_VisSlColno2(inda);
% comp2=new.NO2_VisSlColno2(indb);
% 
% rms1=data.NO2_VisRMS(inda);
% rms2=new.NO2_VisRMS(indb);
% elevs=data.Elevviewingangle(inda);
% 
% rel_diff=((comp2./comp1)-1)*100;
% 
% mm=[];
% sig=[];
% eas=[1,2,3,5,8,10,15,20,30,40,50];
% for i=eas
%    
%     ind=(abs(rel_diff)<100 & (rms1<0.003 & rms2<0.003) & elevs==i);
%     mm=[mm,mean(rel_diff(ind))]; 
%     sig=[sig,std(rel_diff(ind))];
% 
% end
% 
% ind=(abs(rel_diff)<100 & (rms1<0.003 & rms2<0.003));
% disp([mean(rel_diff(ind)), std(rel_diff(ind))])
% 
% plot(eas,mm,'kx-'), hold on
% plot(eas,mm+sig,'kx--')
% plot(eas,mm-sig,'kx--')
% grid on
% xlabel('Elevation angle')
% ylabel('Relative difference (%)')
% legend('mean','std','location','northwest')
% % plot(fd(ind),rel_diff(ind),'kx')

%% BEE stats
% % load('/home/kristof/work/BEEs/BEE_dataset_all.mat')
% % bee_dataset(bee_dataset.times.Year==2015,:)=[];
% 
% % % percentage of each wdir
% % disp((sum(bee_dataset.N_SE_rest==1)/10150)*100)
% % disp((sum(bee_dataset.N_SE_rest==2)/10150)*100)
% % disp((sum(bee_dataset.N_SE_rest==3)/10150)*100)
% % disp((sum(isnan(bee_dataset.N_SE_rest))/10150)*100)
% 
% % percentage of each wdir, for BrO above median
% tmp=prctile(bee_dataset.bro_col,75);
% disp((sum(bee_dataset.N_SE_rest==1 & bee_dataset.bro_col>tmp)/10150)*100)
% disp((sum(bee_dataset.N_SE_rest==2 & bee_dataset.bro_col>tmp)/10150)*100)
% disp((sum(bee_dataset.N_SE_rest==3 & bee_dataset.bro_col>tmp)/10150)*100)
% disp((sum(isnan(bee_dataset.N_SE_rest) & bee_dataset.bro_col>tmp)/10150)*100)

%% GBS yearly plots for campaign meeting
% % load('/home/kristof/work/GBS/VCD_results/NDACC_RD/UT-GBS_O3_VCD_2019all.mat')
% % try 
% %     gscatter(mjd2k_to_date(reanalysis.mjd2k),reanalysis.mean_vcd,reanalysis.ampm,'br','..')
% % end
% 
% try 
%     load('/home/kristof/work/GBS/VCD_results/NDACC_RD/UT-GBS_NO2_VCD_2019all.mat')    
%     gscatter(mjd2k_to_date(reanalysis.mjd2k),reanalysis.mean_vcd,reanalysis.ampm,'br','..')
% end
% try
%     hold on
%     load('/home/kristof/work/GBS/VCD_results/NDACC_RD/PEARL-GBS_NO2_UV_VCD_2019all.mat')   
%     gscatter(mjd2k_to_date(reanalysis.mjd2k),reanalysis.mean_vcd,reanalysis.ampm,'br','oo')
% end
% 
% grid on
% xlim([datenum(2019,2,15),datenum(2019,10,15)])
% xlabel('Date, 2019 (UTC)')
% % ylabel('Ozone VCD (molec/cm^2)')
% % legend('UT-GBS am','UT-GBS pm')
% ylabel('NO_2 VCD (molec/cm^2)')
% legend('UT-GBS am','UT-GBS pm','PEARL-GBS am','PEARL-GBS pm')
% 
% set(findall(gcf,'-property','FontSize'),'FontSize',15)
% set(findall(gcf,'-property','FontName'),'FontName','Arial') 
% set(gcf, 'Position', [100, 100, 1100, 650]);

%% add time info to sea ice contact files, or modify contact time

% area={'FYSI','MYSI','water','land'};
% 
% for i=1:length(area)
%     for j=1:5
%     
%         fname=['FP_' area{i} '_contact_' num2str(j) 'day.mat'];
%         
%         load(fname);
%         
% %         tmp=table();
% % 
% %         tmp.run_times=run_times';
% %         tmp.run_start=run_start';
% %         tmp.run_end=run_end';
% %         tmp.contact=FP_SI_contact;
% % 
% %         FP_SI_contact=tmp;
%         
%         FP_SI_contact.contact=FP_SI_contact.contact/12534^2;
%         
%         save(fname,'FP_SI_contact')
%         
%     end
% end

%% test box plot
% x = [1,2,3,4,5,1,2,3,4,6];
% group = [1,1,2,2,2,3,3,3,4,4];
% positions = [1 1.25 2 2.25];
% boxplot(x,group, 'positions', positions);
% 
% set(gca,'xtick',[mean(positions(1:2)) mean(positions(3:4)) ])
% set(gca,'xticklabel',{'Direct care','Housekeeping'})
% 
% color = ['c', 'y', 'c', 'y'];
% h = findobj(gca,'Tag','Box');
% for j=1:length(h)
%    patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
% end
% 
% c = get(gca, 'Children');
% 
% hleg1 = legend(c(1:2), 'Feature1', 'Feature2' );

%% pandora ozone and NO2 timeseries
% load('/home/kristof/work/documents/conferences/Pandora_workshop_2019/data/all_2019.mat')
% 
% figure
% ut_o3(ut_o3.mjd2k<7020,:)=[];
% ut_no2(ut_no2.mjd2k<7020,:)=[];
% 
% plot(ut_o3.DateTime,ut_o3.mean_vcd,'ks','markersize',8), hold on
% 
% ind=brewer_ds.Serial_NO==111;
% plot(brewer_ds.DateTime(ind),brewer_ds.ColumnO3(ind),'b.','markersize',12)
% 
% plot(Pandora_o3.DateTime,Pandora_o3.O3_VCD,'r.','markersize',12)
% 
% xlabel('Date, 2019 [UTC]')
% ylabel('O_3 VCD [DU]')
% 
% legend('UT-GBS am/pm','Brewer #111','Pandora #144')
% 
% figure
% plot(ut_no2.DateTime,ut_no2.mean_vcd,'ks','markersize',8), hold on
% 
% ind=find(Pandora_no2.RMS<=6e-4 & Pandora_no2.INT<=200);
% 
% Pandora_no2=Pandora_no2(ind,:);
% 
% plot(Pandora_no2.DateTime,Pandora_no2.NO2_VCD,'c.','markersize',6)
% 
% [ft,years]=fracdate(Pandora_no2.DateTime);
% Pandora_no2.mjd2k=ft_to_mjd2k(ft,years);
% 
% Pandora_no2.tot_col=Pandora_no2.NO2_VCD;
% Pandora_no2.tot_col_err=Pandora_no2.NO2_VCD_err;
% Pandora_no2.sza=Pandora_no2.SZA;
% 
% [ vcd1, vcd2, ~, ~, times]=find_coincidences_time( ut_no2, Pandora_no2, 3);
% 
% plot(mjd2k_to_date(times),vcd2,'r.','markersize',12)
% 
% 
% xlabel('Date, 2019 [UTC]')
% ylabel('NO_2 VCD [DU]')
% 
% legend('UT-GBS am/pm','Pandora #144 all','Pandora #144 paired',...
%        'location','southeast')

%% SI contactfor FLEXPART back trajectories
% % for age=[1,2,0,20]
% for age=[1,2]
%     age
%     for bt_len=[1,2]
%         bt_len
%         get_SI_contact(bt_len,age)
%     end
% end

%% check if FLEXPART to EASE grid mapping works
% load('/home/kristof/work/BEEs/sea_ice_data/EASE_grid_SI_age.mat')
% load('/home/kristof/berg/FLEXPART_10.02/grid_data/fine_FLEXPART_mesh.mat')
% 
% % surf(age(:,:,1)','EdgeColor','None', 'facecolor', 'flat');
% % view(2)
% % xlim([0,722])
% % ylim([0,722])
% 
% age=age(:,:,1)';
% lat_age=lat_age';
% lon_age=lon_age';
% lon_age(lon_age<0)=lon_age(lon_age<0)+360;
% 
% % ic=339;
% % ir=196;
% % 
% % lin_ind=ir+((ic-1)*722);
% % 
% % [row,col]=find(FP_fine_mask==lin_ind);
% % disp([mean(lat_FP_fine(row)), mean(lon_FP_fine(col))]);
% % disp([lat_age(ir,ic),lon_age(ir,ic)])
% % disp([lat_age(lin_ind),lon_age(lin_ind)])
% tic
% fysi=find(age==1);
% 
% % disp([mean(lat_age(fysi)), mean(lon_age(fysi))]);
% x=ismember(FP_fine_mask,fysi);
% 
% % disp([mean(lat_FP_mesh(x)),mean(lon_FP_mesh(x))]);
% toc

%% 0.5 deg in longitude as a function of latitude
% 
% R_e = 6378.1;
% circ=2*R_e*pi;
% 
% lats=0:90;
% 
% dists=(2*pi*(R_e*sind(90-lats)))/720;
% 
% areaint([80,80,79.5,79.5],[-86.5,-86,-86,-86.5])*(4*pi*R_e^2)
%
%% plot SI age
%
% load('EASE_grid_SI_age.mat')
% 
% age(age==0)=NaN;
% 
% age(age>1 & age<10)=15;
% 
% for i=1:208
%     
%     if date_age.Month(i)<3
%         continue
%     elseif date_age.Month(i)>5
%         continue
%     end
%     
%     surf(fliplr(age(:,:,i)),'EdgeColor','None', 'facecolor', 'flat');
%     view(2)
%     xlim([0,722])
%     ylim([0,722])
%     title(datestr(date_age(i)))
%     try
%         tmp=1;
%         while tmp % loop so key presses (return 1) are not accepted
%             tmp=waitforbuttonpress;
%         end
%     catch
%         % returns error if figure is closed, exit when that happens
%         return
%     end
% end

%% CINDI-2 profiling version comparison
% 
% species='tracegas';
% 
% figure
% load(['/home/kristof/work/profile_retrievals/profile_results/CINDI-2/matlab_files',...
%      '/medianDSCD_uv_v2_all/',species,'/medianDSCD_uv_v2_all_all.mat'])
% 
% uv_v2=prof_nd(1,:);
%  
% plot(times,prof_nd(1,:),'ko'), hold on 
%  
% load(['/home/kristof/work/profile_retrievals/profile_results/CINDI-2/matlab_files',...
%      '/medianDSCD_uv_v3_all/',species,'/medianDSCD_uv_v3_all_all.mat'])
% 
% uv_v3=prof_nd(1,:);
%  
% plot(times,prof_nd(1,:),'rx'), hold on 
% 
% figure
% load(['/home/kristof/work/profile_retrievals/profile_results/CINDI-2/matlab_files',...
%      '/medianDSCD_vis_v2_all/',species,'/medianDSCD_vis_v2_all_all.mat'])
% 
% vis_v2=prof_nd(1,:);
% 
% plot(times,prof_nd(1,:),'ko'), hold on 
%  
% load(['/home/kristof/work/profile_retrievals/profile_results/CINDI-2/matlab_files',...
%      '/medianDSCD_vis_v3_all/',species,'/medianDSCD_vis_v3_all_all.mat'])
% 
% vis_v3=prof_nd(1,:);
% 
% plot(times,prof_nd(1,:),'rx'), hold on 
% 
% % ind=(times.Day<16 | times.Day>18);
% ind=(times.Day<30);
% 
% mean((uv_v3(ind)./uv_v2(ind))-1)*100
% mean((vis_v3(ind)./vis_v2(ind))-1)*100

%% OClO plots
% 
% % ind=(oclo_2017.SZA>86 & oclo_2017.SZA<92);
% % plot(oclo_2017.Fractionalday(ind),oclo_2017.NO2SlColoclo(ind),'ro'), hold on
% % 
% % ind=(oclo_2018.SZA>86 & oclo_2018.SZA<92);
% % plot(oclo_2018.Fractionalday(ind),oclo_2018.NO2SlColoclo(ind),'bx'), hold on
% % 
% % ind=(oclo_2019.SZA>86 & oclo_2019.SZA<92);
% % plot(oclo_2019.Fractionalday(ind),oclo_2019.NO2SlColoclo(ind),'g+'), hold on
% 
% figure
% ind=(oclo_2017.SZA>86 & oclo_2017.SZA<92);
% plot(oclo_2017.Fractionalday(ind),oclo_2017.NO2RMS(ind),'ro'), hold on
% 
% ind=(oclo_2018.SZA>86 & oclo_2018.SZA<92);
% plot(oclo_2018.Fractionalday(ind),oclo_2018.NO2RMS(ind),'bx'), hold on
% 
% ind=(oclo_2019.SZA>86 & oclo_2019.SZA<92);
% plot(oclo_2019.Fractionalday(ind),oclo_2019.NO2RMS(ind),'g+'), hold on

%% globe plots
% load coast
% figure
% hold on
% ax=axesm('MapProjection','globe');
% axis off
% gridm on
% framem on
% 
% % mlabel on
% plabel on;
% setm(gca,'MLabelParallel',50)
% % setm(gca,'PLabelMeridian',180)
% % setm(gca,'MLabelParallel',0)
% 
% for i=1:36
%     geoshow(ax, [0,90,0,0], [10*i,10*i,10*i-10,10*i],'DisplayType','polygon',...
%             'facecolor','w', 'edgecolor','w')
% end
% 
% geoshow(ax, lat, long,'DisplayType', 'polygon', 'FaceColor', [0.7,0.7,0.7])
% 
% view(4,80)
% zoom(2.7)

%% plot flexpart release boxes
% load coast
% ax = worldmap([79,82], [-92.4,-80.4]);
% geoshow(ax, lat, long,'DisplayType', 'polygon', 'FaceColor', [0.7,0.7,0.7]), hold on
% 
% % Erik, for columns
% plotm([79.55,80.55],[-86.92,-86.92],'k-')
% plotm([79.55,80.55],[-85.92,-85.92],'k-')
% plotm([80.55,80.55],[-85.92,-86.92],'k-')
% plotm([79.55,79.55],[-85.92,-86.92],'k-')
% 
% % Kristof, for MAX-DOAS BrO
% plotm([79.85,80.65],[-87.92,-87.92],'b-')
% plotm([79.85,80.65],[-84.92,-84.92],'b-')
% plotm([80.65,80.65],[-84.92,-87.92],'b-')
% plotm([79.85,79.85],[-84.92,-87.92],'b-')
% 
% plotm(80.053,-86.416,'ro')

%% test reading in flexpart data
% cc={'r','g','b','y','c'};
% for i=0:4
%     plotm(trajectories(:,18+i*5),trajectories(:,17+i*5),cc{i+1},'linewidth',2)
% end

% for i=2:2:72
%     tracer_plot=tracer(:,:,i)';
%     tracer_plot(tracer_plot==0)=NaN;
% 
%     clf
% %     ax = worldmap([70,85], [-126.4,-46.4]);
%     ax = worldmap([70,90], [-146.4,-26.4]);
%     geoshow(ax, lat, long,'DisplayType', 'polygon', 'FaceColor', [0.7,0.7,0.7]);
% 
%     surfm(latitude,longitude,tracer_plot,'facecolor', 'interp','facealpha',0.8);
%     colormap(flipud(hot))
%     cb = colorbar();
%     cb.Ruler.Scale = 'log';
%     cb.Ruler.MinorTick = 'on';
%     
%     % pause loop until figure is clicked on
%     try
%         tmp=1;
%         while tmp % loop so key presses (return 1) are not accepted
%             tmp=waitforbuttonpress;
%         end
%     catch
%         % returns error if figure is closed, exit when that happens
%         return
%     end
%     
% end

%% compare 2011 BrO profiles retrieved using different a priori
% cd('/home/kristof/work/profile_retrievals/profile_results')
% dirnames={'_+1_elev_corr','_double_scaleh','_bolz'};
% plot_symb={'ks','rx','bo','g+','cd'};
% 
% load('/home/kristof/work/profile_retrievals/profile_results/xiaoyi_2011_BrO_cols.mat')
% 
% for i=1:length(dirnames)
%     
%     load(['eureka_2011' dirnames{i} '/tracegas/profiles_2011_filt.mat'])
% %     info(times.Day>2,:)=[];
% %     times(times.Day>2,:)=[];
%     if i==1
%         info_orig=info;
%     else
% %         mean((info.col./info_orig.col)-1)*100
% %         mean((info.col-info_orig.col))
%     end
%     
% %     figure
% %     prof(end,1)=30*1e-6;
% %     surf(ft,alt,prof*1e6,'EdgeColor','None', 'facecolor', 'interp'); hold on
% %     view(2)
% %     colorbar
% %     title(dirnames{i})
% 
%     plot(times,info.DOFS,plot_symb{i},'markersize',8), hold on
% %     plot(times,prof(1,:),plot_symb{i}), hold on
% %     try 
% %         mean(dscd.O4meas(dscd.O4meas~=0)-dscd.O4retr(dscd.O4meas~=0))
% %     catch
% %         mean(dscd.meas(dscd.meas~=0)-dscd.retr(dscd.meas~=0))
% %     end
% 
% end
% 
% % plot(xiaoyi_all.DateTime,xiaoyi_all.BrO_VCDmoleccm2,plot_symb{i+1},'markersize',8), hold on
% 
% ll={'my baseline','double scale h','box apriori','from Xiaoyi'};
% legend(ll,'location','best')
% 
% xlabel('Date, 2011 (UTC)')
% ylabel('BrO VCD (molec/cm^2)')
% 
% grid on
% 
% set(gcf, 'Position', [100, 100, 1000, 600])

%% compare 20110404 BrO profiles retrieved using different a priori
% cd('/home/kristof/work/profile_retrievals/profile_results/eureka_2011_tests/tracegas')
% fnames={'_orig','_2km_sh_250m_res','_3km_sh_250m_res','_15ppt','_bolz'};
% plot_symb={'ks','rx','bo','g+','cd'};
% 
% for i=1:length(fnames)
%     
%     load(['20110404' fnames{i} '.mat'])
%     
%     if i==1
%         info_orig=info;
%     else
% %         mean((info.col./info_orig.col)-1)*100
%     end
%     
%     figure
%     prof(end,1)=30*1e-6;
%     surf(ft,alt,prof*1e6,'EdgeColor','None', 'facecolor', 'interp'); hold on
%     view(2)
%     colorbar
%     title(fnames{i})
% 
% %     plot(times,info.col,plot_symb{i}), hold on
%     % plot(times,prof(1,:),plot_symb{i}), hold on
%     % mean(dscd.O4meas(dscd.O4meas~=0)-dscd.O4retr(dscd.O4meas~=0))
% 
% end

% ll={'orig','2km sh','3km sh','15ppt','bolz'};
% legend(ll)

%% compare G and K from my retrieval to one of Xiaoyi's old versions (wholerun 110404)
% new=gain201104041604(:,2:end)*wf201104041604(:,2:end)';
% 
% figure
% for i=1:20
%     plot(new(:,i),[0:0.2:3.8]), hold on
% end
% 
% old=gain201104041610(:,2:end)*wf201104041610(:,2:end)';
% 
% figure
% for i=1:20
%     plot(old(:,i),[0:0.2:3.8]), hold on
% end
% 
% tmp=gain201104041610(:,2:end)*wf201104041604(:,2:end)';
% 
% figure
% for i=1:20
%     plot(tmp(:,i),[0:0.2:3.8]), hold on
% end

%% plot Xiaoyi's 2011 BrO
% load('/home/kristof/work/profile_retrievals/profile_results/xiaoyi_2011_BrO_cols.mat')
% plot(BrOretrieval20110401.Time,BrOretrieval20110401.BrO_VCDmoleccm2,'rx'), hold on
% plot(BrOretrieval20110402.Time,BrOretrieval20110402.BrO_VCDmoleccm2,'rx')
% plot(BrOretrieval20110403.Time,BrOretrieval20110403.BrO_VCDmoleccm2,'rx')
% plot(BrOretrieval20110404.Time,BrOretrieval20110404.BrO_VCDmoleccm2,'rx')
% plot(BrOretrieval20110405.Time,BrOretrieval20110405.BrO_VCDmoleccm2,'rx')


%% write BrO profiles as text file

% % load('profiles_2018_filt.mat');
% % load('tracegas_profiles_2013_normal_aer_only.mat')
% load('profiles_2011_filt.mat');
% 
% prof=prof';
% prof_err=prof_err';
% prof_nd=prof_nd';
% prof_nd_err=prof_nd_err';
% 
% for i=1:length(alt)
%     
%     prof_head{i}=['prof_' num2str(alt(i)*1000) 'm'];
%     prof_err_head{i}=['prof_err_' num2str(alt(i)*1000) 'm'];
%     
% end
% 
% out=array2table([ft,info.DOFS,info.col,info.col_err,prof*1e6,prof_err*1e6],...
%                 'VariableNames',[{'Fractional_time','DOFS','column','column_err'},...
%                 prof_head, prof_err_head]);
%             
% writetable(out,'Eureka_BrO_profiles_ppt_2011','delimiter',',');
% 
% add manually to file:
% # Fractional time is 0 at Jan 1, 00:00
% # Columns are 0-4 km, in molec/cm^2
% # profiles and errors are parts per trillion by volume

%% compare OSIRIS v3 to v6
% 
% figure()
% hold on
% for i=2002:2013
%     v3=sum(osiris_v3.year==i);
%     v6=sum(osiris_v6.year==i);
%     plot(i,v3,'rx')
%     plot(i,v6,'bx')
%     disp([num2str(v3/v6)])
% end
% 
% legend('v3','v6')
% ylim([0,1050])
% xlabel('Year')
% ylabel('N.o. part. cols near Eureka')
% grid on

%% compare DMP values for ACE/DOAS at select altitudes
% day_range=[40.25,80.25];
% ind=find(gbs_o3.fractional_time<day_range(1)-1 | gbs_o3.fractional_time>day_range(2));
% gbs_o3(ind,:)=[];
% 
% figure
% subplot(211)
% plot([0,7000],[1.2,1.2]*1e-4,'k-')
% hold on
% plot([0,7000],[1.6,1.6]*1e-4,'k-')
% 
% plot(ace_fts_o3.mjd2k,ace_fts_o3.spv(:,2),'ko')
% plot(gbs_o3.mjd2k,gbs_o3.spv(:,2),'rx')
% 
% ylim([0.5,2.5]*1e-4)
% 
% subplot(212)
% plot([0,7000],[1.2,1.2]*1e-4,'k-')
% hold on
% plot([0,7000],[1.6,1.6]*1e-4,'k-')
% 
% plot(ace_fts_o3.mjd2k,ace_fts_o3.spv(:,5),'ko')
% plot(gbs_o3.mjd2k,gbs_o3.spv(:,5),'rx')
% 
% ylim([0.5,2.5]*1e-4)

%% plot GBS DMPs by year

% year=2016;
% 
% load(['DOAS_O3_VIS_DMP_table_' num2str(year) '.mat'])
% figure(20)
% clf
% fractional_time(fractional_time>80)=[];
% m_spv=zeros(60,1);
% 
% plot([1.2,1.2]*1e-4,[.8,59.5],'k-','linewidth',1.5), hold on
% plot([1.6,1.6]*1e-4,[.8,59.5],'k-','linewidth',1.5)
% 
% for i=1:length(fractional_time)
%     m_spv=m_spv+dmp_all{i}.spv;
%     plot(dmp_all{i}.spv,dmp_all{i}.alt,'color',[0 0.7 0.7])
% end
% m_spv=m_spv./length(fractional_time); 
% 
% plot(m_spv,fliplr([.8,1.5:59.5]),'k-')
% ylim([13,52])

%% plot ACE DMPs by year
% figure(21)
% clf
% ind=find(dmpstruct.year==year & dmpstruct.fractional_time<80 & ...
%          dmpstruct.lat(30,:)'>75 & dmpstruct.lon(30,:)'>-120 & dmpstruct.lon(30,:)'<-60);
% 
% plot([1.2,1.2]*1e-4,[.8,59.5],'k-','linewidth',1.5), hold on
% plot([1.6,1.6]*1e-4,[.8,59.5],'k-','linewidth',1.5)
% 
% for i=ind
%     
%     plot(dmpstruct.spv(:,i),dmpstruct.altitude_km(:,i),'color',[0 0.7 0.7])
%     
% end
% 
% plot(mean(dmpstruct.spv(:,ind),2),dmpstruct.altitude_km(:,1),'k-')
% ylim([13,52])

%% match ACE DMP bad values to GLC bad values
% 
% [~,ind_dmp,ind_glc]=intersect(dmpstruct.mjd2k, glcstruct.mjd2k);
% 
% dmp_lat=dmpstruct.lat(:,ind_dmp);
% dmp_lon=dmpstruct.lon(:,ind_dmp);
% 
% glc_lat=glcstruct.lat(:,ind_glc);
% glc_lon=glcstruct.lon(:,ind_glc);
% 
% bad_lat=find(abs(mean(glc_lat))>90 | sum(isnan(glc_lat))==150);
% bad_lon=find(abs(mean(glc_lon))>180 | sum(isnan(glc_lon))==150);
% 
% bad_lat2=find(abs(mean(dmp_lat))>90 | sum(isnan(dmp_lat))==75);
% bad_lon2=find(abs(mean(dmp_lon))>180 | sum(isnan(dmp_lon))==75);
% 
% nan_t=find(sum(isnan(dmpstruct.T))==75);
%
%% match all ACE data to DMPs to plot coords of NaN profiles
% nan_t=find(sum(isnan(dmpstruct.T))==75);
% nan_mjd2k=dmpstruct.mjd2k(nan_t);
%
% date=tanstruct.date_mjd-date2mjd(2000,1,1,0,0,0);
% 
% for i=1:length(nan_mjd2k)
%     [tmp,ind]=sort(abs(date-nan_mjd2k(i)));
%     if tmp(1)<1e-4 
%         ind=ind(1);
%         ace_lat(i)=(tanstruct.lat_tangent(ind));
%         ace_lon(i)=(tanstruct.lon_tangent(ind));
%     end
% end
%
% figure()
% load coast
% ax = worldmap([-90,90], [-180,180]);
% geoshow(ax, lat, long,'DisplayType', 'polygon', 'FaceColor', [0.7,0.7,0.7])
% 
% plotm(ace_lat,ace_lon, 'k.')

%% merge UT-GBS VCD batches for VCD reprocessing (for Gaia)
% vt_new=table;
% rcd_new=struct;
% rcd_new.R2=[];
% rcd_new.fd_min=[];
% rcd_new.fd_max=[];
% rcd_new.sza_max=[];
% rcd_new.sza_min=[];
% rcd_new.nbr=[];
% 
% dscd_new=struct;
% dscd_new.fd=[];
% dscd_new.tot_tint=[];
% dscd_new.amf=[];
% 
% for i=11:12
%     
%     load(['UT-GBS_NO2_VCD_2018_' num2str(i) '.mat'])
%     
%     vt_new=[vt_new; VCD_table];
%     rcd_new.R2=[rcd_new.R2; rcd_S.R2];
%     rcd_new.fd_min=[rcd_new.fd_min; rcd_S.fd_min];
%     rcd_new.fd_max=[rcd_new.fd_max; rcd_S.fd_max];
%     rcd_new.sza_min=[rcd_new.sza_min; rcd_S.sza_min];
%     rcd_new.sza_max=[rcd_new.sza_max; rcd_S.sza_max];
%     rcd_new.nbr=[rcd_new.nbr; rcd_S.nbr];
%     
%     dscd_new.fd=[dscd_new.fd; dscd_S.fd];
%     dscd_new.tot_tint=[dscd_new.tot_tint; dscd_S.tot_tint];
%     dscd_new.amf=[dscd_new.amf; dscd_S.amf];
%     
% end
% 
% 
% VCD_table=vt_new;
% rcd_S=rcd_new;
% dscd_S=dscd_new;
% 
% clearvars -except rcd_S VCD_table dscd_S
% 
% save UT-GBS_NO2_VCD_2018_99.mat

%% plot bruker no2 profiles
% ind=find(bruker.sza>=60);
% 
% out=[];
% figure()
% for i=ind'
%     
%     pks = findpeaks(prof(i,:));
%     
%     if length(pks)<=2, continue, end
%     
%     out=[out;[i,bruker.sza(i)]];
%     
%     plot(prof(i,:), alt_km, 'k-'), hold on
%     plot(apriori(i,:), alt_km, 'r--')
% 
%     title(['i = ' num2str(i) ', SZA = ' num2str(bruker.sza(i))])
%     xlim([0,2.5]*1e9)
%     ylim([0,60])
%     legend('profile','a priori')
%     
%     waitforbuttonpress;
% %     pause(0.1)
%     cla;
% 
% end


%% check VCD data filter for each year -- SAOZ data
% % load data manually
% figure(99)
% clf
% [ind_goodvcd,VCD_table_new]=filter_VCD_output( 4, VCD_table, rcd_S, 3, false);
% 
% plot(VCD_table.fd,VCD_table.mean_vcd,'ko'), hold on
% plot(VCD_table_new.fd(ind_goodvcd),VCD_table_new.mean_vcd(ind_goodvcd),'rx',...
%      'linewidth',2,'markersize',10), hold on

% check VCD data filter for each year
% load data manually
% figure(98)
% clf
% [ind_goodvcd,VCD_table_new]=filter_VCD_output( 1, VCD_table, rcd_S, 1, true, true);
% if isempty(ind_goodvcd), return, end
% 
% plot(VCD_table.fd,VCD_table.mean_vcd,'ko'), hold on
% plot(VCD_table_new.fd(ind_goodvcd),VCD_table_new.mean_vcd(ind_goodvcd),'rx',...
%      'linewidth',2,'markersize',10), hold on
%  
% sys=(VCD_table_new.sigma_mean_vcd(ind_goodvcd)./VCD_table_new.mean_vcd(ind_goodvcd))*100;
% rand=(VCD_table_new.std_vcd(ind_goodvcd)./VCD_table_new.mean_vcd(ind_goodvcd))*100;

% % figure(2)
% % subplot(311)
% % plot(VCD_table_new.fd(ind_goodvcd),sys, 'kx')
% % subplot(312)
% % plot(VCD_table_new.fd(ind_goodvcd),rand, 'kx')
% % subplot(313)
% % plot(VCD_table_new.fd(ind_goodvcd),sqrt(sys.^2 + rand.^2), 'kx')
% % 
% % median(sys)
% % median(rand)

%% correct ozone VIS dscds for 2006-2010 PGBS data
% data.orig_O3SlColo3=data.O3SlColo3;
% data.O3SlColo3=data.O3SlColo3+(data.O3SlColX*4.4e18);

%% find files that use function fn in current directory
% fn = 'mjd2k_ft'; % function you want track down
% pp = pwd(); % project path
% d = dir([pp filesep '*.m']); % list M-files in folder. Not recursive
% for k = 1:numel(d);
%    file = d(k).name;
%    h = getcallinfo(file);
%    for n = 1:numel(h) % distinguish calls from subfunctions etc
%       name = h(n).functionPrefix;
%       lines = h(n).calls.fcnCalls.lines(strcmp(h(n).calls.fcnCalls.names, fn));
%       if lines
%          disp(['Function ' fn ' is called by function ' name ' at lines ' num2str(lines)])
%       end
%    end
% end



%% rename table headers
% dscd_uv=NO2uv;
% dscd_uv.Properties.VariableNames{7} = 'NO2_dSCD';
% dscd_uv.Properties.VariableNames{8} = 'NO2_dSCD_Err';
% dscd_uv=[dscd_uv,O4uv(:,7:8)];
% dscd_uv.Properties.VariableNames{9} = 'O4_dSCD';
% dscd_uv.Properties.VariableNames{10} = 'O4_dSCD_Err';
% dscd_uv.NO2_dSCD=dscd_uv.NO2_dSCD*1e15;
% dscd_uv.NO2_dSCD_Err=dscd_uv.NO2_dSCD_Err*1e15;
% dscd_uv.O4_dSCD=dscd_uv.O4_dSCD*1e42;
% dscd_uv.O4_dSCD_Err=dscd_uv.O4_dSCD_Err*1e42;


%% stuff
% for i=15:30:365
%     x=albedo_lookup(80.05, -86.4, 2016, i);
%     y=albedo_lookup(-80.05, -86.4, 2016, i);
%     disp([x,y])
% end



% 
% plot(qdoas_31(:,col_31.fd),(qdoas_31(:,col_31.o3_dscd)-...
%      qdoas_2109(:,col_2109.o3_dscd))./qdoas_31(:,col_31.o3_dscd), 'k.'), hold on
% 

% % % figure(99)
% % % for day=7:21
% % %     fd=x-59+1;
% % % 
% % %     ind=find(fd>day & fd<day+1);
% % % 
% % %     for i=ind
% % %         plot(T_arr(1:400,i)-273,y(1:400)), hold on
% % %     end
% % %     plot([-50,-15],[610,610],'k-')
% % %     
% % %     xlim([-50,-15])
% % %     ylim([0,2000])
% % %     title(['March ' num2str(day)])
% % % 
% % %     k = waitforbuttonpress;
% % %     clf
% % % end
% % % 
% % % figure(1)
% % % subplot(212)


% fd=x-58;
% 
% 
% surf(x,y(1:200)./1000,T_arr(1:200,:)-273,'edgecolor','none')
% 
% view(2)
% colormap('jet')
% colorbar
% 
% xlim([5,25])
% ylim([0,2])

% % % convert times to fractional time
% % [~,ft_opc]=fracdate(time,'dd-eee-yyyy hh:mm:ss');
% % 
% % sm=sum(Dp_data(:,3:5),2);
% % 
% % % bad data, 2017
% % % sm(7554:7567)=[];
% % % ft_opc(7554:7564)=[];
% % 
% % ind=find(sm<1100 & sm>0);
% % 
% % [~,ft_opc]=fracdate(time,'dd-eee-yyyy hh:mm:ss');
% % plot(ft_opc(ind)-59,sm(ind),'linewidth',1.2)
% % 
% % % xlim([7,21.5])
% % xlim([19,23.5])
% % ylim([0,950])
% % ylabel('D_p>1\mum (cm^-^3)')
% % grid on
% % 
