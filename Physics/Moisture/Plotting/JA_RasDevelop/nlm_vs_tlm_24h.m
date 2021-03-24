close all
clear
clc

load coast
coast_lat = lat; clear lat
coast_lon = long; clear long

cd /discover/nobackup/drholdaw/ExperimentData/Journal_Articles/moist_dev_mwr/

lon = ncread('fvpert.eta.dry.20120318_00z.nc4','lon');
lat = ncread('fvpert.eta.dry.20120318_00z.nc4','lat');

U_dry  = ncread('fvpert.eta.dry.20120318_00z.nc4','u');
V_dry  = ncread('fvpert.eta.dry.20120318_00z.nc4','v');
TV_dry = ncread('fvpert.eta.dry.20120318_00z.nc4','tv');
DP_dry = ncread('fvpert.eta.dry.20120318_00z.nc4','delp');
QV_dry = ncread('fvpert.eta.dry.20120318_00z.nc4','sphu');
QL_dry = ncread('fvpert.eta.dry.20120318_00z.nc4','qltot');
QI_dry = ncread('fvpert.eta.dry.20120318_00z.nc4','qitot');
O3_dry = ncread('fvpert.eta.dry.20120318_00z.nc4','ozone');

U_moi  = ncread('fvpert.eta.moi.20120318_00z.nc4','u');
V_moi  = ncread('fvpert.eta.moi.20120318_00z.nc4','v');
TV_moi = ncread('fvpert.eta.moi.20120318_00z.nc4','tv');
DP_moi = ncread('fvpert.eta.moi.20120318_00z.nc4','delp');
QV_moi = ncread('fvpert.eta.moi.20120318_00z.nc4','sphu');
QL_moi = ncread('fvpert.eta.moi.20120318_00z.nc4','qltot');
QI_moi = ncread('fvpert.eta.moi.20120318_00z.nc4','qitot');
O3_moi = ncread('fvpert.eta.moi.20120318_00z.nc4','ozone');

U_nlm_free  = ncread('iau590.prog.eta.free.20120318_00z.nc4','u');
V_nlm_free  = ncread('iau590.prog.eta.free.20120318_00z.nc4','v');   
TV_nlm_free = ncread('iau590.prog.eta.free.20120318_00z.nc4','tv');
DP_nlm_free = ncread('iau590.prog.eta.free.20120318_00z.nc4','delp');
QV_nlm_free = ncread('iau590.prog.eta.free.20120318_00z.nc4','sphu');
QL_nlm_free = ncread('iau590.prog.eta.free.20120318_00z.nc4','qltot');
QI_nlm_free = ncread('iau590.prog.eta.free.20120318_00z.nc4','qitot');
O3_nlm_free = ncread('iau590.prog.eta.free.20120318_00z.nc4','ozone');

U_nlm_del_0_2  = ncread('iau590.prog.eta.del1p0.20120318_00z.nc4','u');
V_nlm_del_0_2  = ncread('iau590.prog.eta.del1p0.20120318_00z.nc4','v');   
TV_nlm_del_0_2 = ncread('iau590.prog.eta.del1p0.20120318_00z.nc4','tv');
DP_nlm_del_0_2 = ncread('iau590.prog.eta.del1p0.20120318_00z.nc4','delp');
QV_nlm_del_0_2 = ncread('iau590.prog.eta.del1p0.20120318_00z.nc4','sphu');
QL_nlm_del_0_2 = ncread('iau590.prog.eta.del1p0.20120318_00z.nc4','qltot');
QI_nlm_del_0_2 = ncread('iau590.prog.eta.del1p0.20120318_00z.nc4','qitot');
O3_nlm_del_0_2 = ncread('iau590.prog.eta.del1p0.20120318_00z.nc4','ozone');

U_nlm = U_nlm_del_0_2 - U_nlm_free;
V_nlm = V_nlm_del_0_2 - V_nlm_free;
TV_nlm = TV_nlm_del_0_2 - TV_nlm_free;
DP_nlm = DP_nlm_del_0_2 - DP_nlm_free;
QV_nlm = QV_nlm_del_0_2 - QV_nlm_free;
QL_nlm = QL_nlm_del_0_2 - QL_nlm_free;
QI_nlm = QI_nlm_del_0_2 - QI_nlm_free;
O3_nlm = O3_nlm_del_0_2 - O3_nlm_free;

cd /home/drholdaw/Lin_Moist_Physics/Moist_Dev_MWR/

%Make center of colormap white
cmap = colormap;
cmap(32,:) = [1 1 1];
cmap(33,:) = [1 1 1];
close

fontsize = 11;

line_wid_cont = 1.0;
line_wid_det = 0.6;

grey = 0.75;

plotlevel = 50;


plot_fields = 0;
plot_diffrs = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   PLOT THE ACTUAL FIELDS   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_fields == 1
    
    TVdiffmax = max(max([TV_nlm(:,:,plotlevel); TV_dry(:,:,plotlevel); TV_moi(:,:,plotlevel)]));
    TVdiffmin = min(min([TV_nlm(:,:,plotlevel); TV_dry(:,:,plotlevel); TV_moi(:,:,plotlevel)]));
    QVdiffmax = max(max([QV_nlm(:,:,plotlevel); QV_dry(:,:,plotlevel); QV_moi(:,:,plotlevel)]));
    QVdiffmin = min(min([QV_nlm(:,:,plotlevel); QV_dry(:,:,plotlevel); QV_moi(:,:,plotlevel)]));

    TVdiffabsmax = max(abs([TVdiffmax TVdiffmin]));
    QVdiffabsmax = max(abs([QVdiffmax QVdiffmin]));
    
    cint_tv = 2*TVdiffabsmax/10;
    cint_qv = 2*QVdiffabsmax/15;
    
    figure('Position',[ 164 30 1115 889])
    
    subplot(3,2,1)
    [C,h] = contour(lon,lat,TV_nlm(:,:,plotlevel)','LineWidth',line_wid_cont);
    caxis([-TVdiffabsmax TVdiffabsmax])
    
    colormap(cmap)
    pos1 = get(gca,'Position');
    colbar = colorbar;
    set(gca,'Position',[pos1(1)-0.05 pos1(2) pos1(3) pos1(4)])
    title('T_v Nonlinear Difference','FontSize',fontsize,'FontName','TimesNewRoman')
    title(colbar,'(K)','FontSize',fontsize,'FontName','TimesNewRoman')
    hold on; 
    plot(coast_lon,coast_lat,'Color',[grey grey grey])
    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);
    set(gca,'XTickLabel',{'180W', '120W', '60W', '0', '60E', '120E', '180E'}...
       ,'YtickLabel',{'90S' '60S' '30S' '0' '30N' '60N' '90N'});
    set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

    set(h,'LevelStep',cint_tv);
    
    subplot(3,2,3)
    [C,h] = contour(lon,lat,TV_dry(:,:,plotlevel)','LineWidth',line_wid_cont);
    caxis([-TVdiffabsmax TVdiffabsmax])
    
    colormap(cmap)
    pos2 = get(gca,'Position');
    colbar = colorbar;
    set(gca,'Position',[pos2(1)-0.05 pos2(2) pos2(3) pos2(4)])
    title('T_v TLM with Dry Physics','FontSize',fontsize,'FontName','TimesNewRoman')
    title(colbar,'(K)','FontSize',fontsize,'FontName','TimesNewRoman')
    hold on; 
    plot(coast_lon,coast_lat,'Color',[grey grey grey])
    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);
    set(gca,'XTickLabel',{'180W', '120W', '60W', '0', '60E', '120E', '180E'}...
       ,'YtickLabel',{'90S' '60S' '30S' '0' '30N' '60N' '90N'});
    set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

    set(h,'LevelStep',cint_tv);
    
    subplot(3,2,5)
    [C,h] = contour(lon,lat,TV_moi(:,:,plotlevel)','LineWidth',line_wid_cont);
    caxis([-TVdiffabsmax TVdiffabsmax])
    
    colormap(cmap)
    pos2 = get(gca,'Position');
    colbar = colorbar;
    set(gca,'Position',[pos2(1)-0.05 pos2(2) pos2(3) pos2(4)])
    title('T_v NLM with Moist Physics','FontSize',fontsize,'FontName','TimesNewRoman')
    title(colbar,'(K)','FontSize',fontsize,'FontName','TimesNewRoman')
    hold on; 
    plot(coast_lon,coast_lat,'Color',[grey grey grey])
    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);
    set(gca,'XTickLabel',{'180W', '120W', '60W', '0', '60E', '120E', '180E'}...
       ,'YtickLabel',{'90S' '60S' '30S' '0' '30N' '60N' '90N'});
    set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

    set(h,'LevelStep',cint_tv);
 
 
    subplot(3,2,2)
    [C,h] = contour(lon,lat,QV_nlm(:,:,plotlevel)','LineWidth',line_wid_cont);
    caxis([-QVdiffabsmax QVdiffabsmax])
    
    colormap(cmap)
    pos3 = get(gca,'Position');
    colbar = colorbar;
    set(gca,'Position',[pos3(1) pos3(2) pos3(3) pos3(4)])
    title('q Nonlinear Difference','FontSize',fontsize,'FontName','TimesNewRoman')
    title(colbar,'(kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
    hold on; 
    plot(coast_lon,coast_lat,'Color',[grey grey grey])
    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);
    set(gca,'XTickLabel',{'180W', '120W', '60W', '0', '60E', '120E', '180E'}...
       ,'YtickLabel',{'90S' '60S' '30S' '0' '30N' '60N' '90N'});
    set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

    set(h,'LevelStep',cint_qv);

    subplot(3,2,4)
    [C,h] = contour(lon,lat,QV_dry(:,:,plotlevel)','LineWidth',line_wid_cont);
    caxis([-QVdiffabsmax QVdiffabsmax])
    
    colormap(cmap)
    pos4 = get(gca,'Position');
    colbar = colorbar;
    set(gca,'Position',[pos4(1) pos4(2) pos4(3) pos4(4)])
    title('q TLM with Dry Physics','FontSize',fontsize,'FontName','TimesNewRoman')
    title(colbar,'(kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
    hold on; 
    plot(coast_lon,coast_lat,'Color',[grey grey grey])
    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);
    set(gca,'XTickLabel',{'180W', '120W', '60W', '0', '60E', '120E', '180E'}...
       ,'YtickLabel',{'90S' '60S' '30S' '0' '30N' '60N' '90N'});
    set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

    set(h,'LevelStep',cint_qv);    
    
    subplot(3,2,6)
    [C,h] = contour(lon,lat,QV_moi(:,:,plotlevel)','LineWidth',line_wid_cont);
    caxis([-QVdiffabsmax QVdiffabsmax])
    
    colormap(cmap)
    pos4 = get(gca,'Position');
    colbar = colorbar;
    set(gca,'Position',[pos4(1) pos4(2) pos4(3) pos4(4)])
    title('q TLM with Moist Physics','FontSize',fontsize,'FontName','TimesNewRoman')
    title(colbar,'(kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
    hold on; 
    plot(coast_lon,coast_lat,'Color',[grey grey grey])
    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);
    set(gca,'XTickLabel',{'180W', '120W', '60W', '0', '60E', '120E', '180E'}...
       ,'YtickLabel',{'90S' '60S' '30S' '0' '30N' '60N' '90N'});
    set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

    set(h,'LevelStep',cint_qv);


end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    PLOT THE DIFFERENCES    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_diffrs == 1
    
    dry_tlm_TV_diff = TV_nlm(:,:,plotlevel)'-TV_dry(:,:,plotlevel)';
    moi_tlm_TV_diff = TV_nlm(:,:,plotlevel)'-TV_moi(:,:,plotlevel)';
    dry_tlm_QV_diff = QV_nlm(:,:,plotlevel)'-QV_dry(:,:,plotlevel)';
    moi_tlm_QV_diff = QV_nlm(:,:,plotlevel)'-QV_moi(:,:,plotlevel)';
    
%     dry_tlm_QV_diff_zone_av = zeros(181,72);
%     moi_tlm_QV_diff_zone_av = zeros(181,72);
%     nlm_tlm_QV_diff_zone_av = zeros(181,72);
%     A = zeros(181,72); B = zeros(181,72); C = zeros(181,72);
%     for i = 1:288
%         A(:,:) = QV_nlm(i,:,:);
%         B(:,:) = QV_dry(i,:,:);
%         C(:,:) = QV_moi(i,:,:);
%         nlm_tlm_QV_diff_zone_av = nlm_tlm_QV_diff_zone_av + (A);
%         dry_tlm_QV_diff_zone_av = dry_tlm_QV_diff_zone_av + (A-B);
%         moi_tlm_QV_diff_zone_av = moi_tlm_QV_diff_zone_av + (A-C);
%     end
%     nlm_tlm_QV_diff_zone_av = nlm_tlm_QV_diff_zone_av/i;
%     dry_tlm_QV_diff_zone_av = dry_tlm_QV_diff_zone_av/i;
%     moi_tlm_QV_diff_zone_av = moi_tlm_QV_diff_zone_av/i;
%     
% %     QVdiffmaxzone = max(max([nlm_tlm_QV_diff_zone_av; dry_tlm_QV_diff_zone_av; moi_tlm_QV_diff_zone_av]));
% %     QVdiffminzone = min(min([nlm_tlm_QV_diff_zone_av; dry_tlm_QV_diff_zone_av; moi_tlm_QV_diff_zone_av]));
%     QVdiffmaxzone = max(max([dry_tlm_QV_diff_zone_av; moi_tlm_QV_diff_zone_av]));
%     QVdiffminzone = min(min([dry_tlm_QV_diff_zone_av; moi_tlm_QV_diff_zone_av]));
%     QVdiffabsmaxzone = max(abs([QVdiffmaxzone QVdiffminzone]));
% 
%     figure
%     [C,h] = contour((nlm_tlm_QV_diff_zone_av)');
%     colormap(cmap)
%     colorbar
%     set(gca,'YDir','reverse')
%     caxis([-QVdiffabsmaxzone QVdiffabsmaxzone])
%     ylim([45 72])
%     xlim([35 145])
%     
%     figure
%     [C,h] = contour((dry_tlm_QV_diff_zone_av)');
%     colormap(cmap)
%     colorbar
%     set(gca,'YDir','reverse')
%     caxis([-QVdiffabsmaxzone QVdiffabsmaxzone])
%     ylim([45 72])
%     xlim([35 145])
%     
%     figure
%     [C,h] = contour((moi_tlm_QV_diff_zone_av)');
%     colormap(cmap)
%     colorbar
%     set(gca,'YDir','reverse')
%     caxis([-QVdiffabsmaxzone QVdiffabsmaxzone])
%     ylim([45 72])
%     xlim([35 145])
%     
%     asd
    
    TVdiffmax = max(max([dry_tlm_TV_diff; moi_tlm_TV_diff]));
    TVdiffmin = min(min([dry_tlm_TV_diff; moi_tlm_TV_diff]));
    QVdiffmax = max(max([dry_tlm_QV_diff; moi_tlm_QV_diff]));
    QVdiffmin = min(min([dry_tlm_QV_diff; moi_tlm_QV_diff]));

    TVdiffabsmax = max(abs([TVdiffmax TVdiffmin]));
    QVdiffabsmax = max(abs([QVdiffmax QVdiffmin]));
    
    cint_tv = 2*TVdiffabsmax/10;
    cint_qv = 2*QVdiffabsmax/25;
    
    figure('Position',[  172         593        1058         256])
%     subplot(2,2,1)
%     [C,h] = contour(lon,lat,dry_tlm_TV_diff,'LineWidth',line_wid_cont);
%     caxis([-TVdiffabsmax TVdiffabsmax])
%     
%     colormap(cmap)
%     pos1 = get(gca,'Position');
%     colbar = colorbar;
%     set(gca,'Position',[pos1(1)-0.05 pos1(2) pos1(3) pos1(4)])
%     title('T_v Difference (Dry Physics)','FontSize',fontsize,'FontName','TimesNewRoman')
%     title(colbar,'(K)','FontSize',fontsize,'FontName','TimesNewRoman')
%     hold on; 
%     plot(coast_lon,coast_lat,'Color',[grey grey grey])
%     axis equal; box on; xlim([-180 180]); ylim([-90 90])
%     set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);
%     set(gca,'XTickLabel',{'180W', '120W', '60W', '0', '60E', '120E', '180E'}...
%        ,'YtickLabel',{'90S' '60S' '30S' '0' '30N' '60N' '90N'});
%     set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% 
%     set(h,'LevelStep',cint_tv);
%     
%     subplot(2,2,3)
%     [C,h] = contour(lon,lat,moi_tlm_TV_diff,'LineWidth',line_wid_cont);
%     caxis([-TVdiffabsmax TVdiffabsmax])
%     
%     colormap(cmap)
%     pos2 = get(gca,'Position');
%     colbar = colorbar;
%     set(gca,'Position',[pos2(1)-0.05 pos2(2) pos2(3) pos2(4)])
%     title('T_v Difference (Moist Physics)','FontSize',fontsize,'FontName','TimesNewRoman')
%     title(colbar,'(K)','FontSize',fontsize,'FontName','TimesNewRoman')
%     hold on; 
%     plot(coast_lon,coast_lat,'Color',[grey grey grey])
%     axis equal; box on; xlim([-180 180]); ylim([-90 90])
%     set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);
%     set(gca,'XTickLabel',{'180W', '120W', '60W', '0', '60E', '120E', '180E'}...
%        ,'YtickLabel',{'90S' '60S' '30S' '0' '30N' '60N' '90N'});
%     set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% 
%     set(h,'LevelStep',cint_tv);
    
 
    subplot(1,2,1)
    [C,h] = contour(lon,lat,dry_tlm_QV_diff,'LineWidth',line_wid_cont);
    caxis([-QVdiffabsmax QVdiffabsmax])
    
    colormap(cmap)
    pos3 = get(gca,'Position');
    colbar = colorbar;
    set(gca,'Position',[pos3(1)-0.05 pos3(2) pos3(3) pos3(4)])
    title('(a) q Difference (Dry Physics)','FontSize',fontsize,'FontName','TimesNewRoman')
    title(colbar,'(kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
    hold on; 
    plot(coast_lon,coast_lat,'Color',[grey grey grey])
    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);
    set(gca,'XTickLabel',{'180W', '120W', '60W', '0', '60E', '120E', '180E'}...
       ,'YtickLabel',{'90S' '60S' '30S' '0' '30N' '60N' '90N'});
    set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

    set(h,'LevelStep',cint_qv);

    subplot(1,2,2)
    [C,h] = contour(lon,lat,moi_tlm_QV_diff,'LineWidth',line_wid_cont);
    caxis([-QVdiffabsmax QVdiffabsmax])
    
    colormap(cmap)
    pos4 = get(gca,'Position');
    colbar = colorbar;
    set(gca,'Position',[pos4(1) pos4(2) pos4(3) pos4(4)])
    title('(b) q Difference (Moist Physics)','FontSize',fontsize,'FontName','TimesNewRoman')
    title(colbar,'(kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
    hold on; 
    plot(coast_lon,coast_lat,'Color',[grey grey grey])
    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);
    set(gca,'XTickLabel',{'180W', '120W', '60W', '0', '60E', '120E', '180E'}...
       ,'YtickLabel',{'90S' '60S' '30S' '0' '30N' '60N' '90N'});
    set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

    set(h,'LevelStep',cint_qv);
    


end

saveas(gcf,'tlm_vs_nlm_24h.eps', 'psc2')

