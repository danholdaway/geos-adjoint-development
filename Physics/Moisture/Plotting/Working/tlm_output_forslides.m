close all
clear
clc

load coast
coast_lat = lat; clear lat
coast_lon = long; clear long

% cd /discover/nobackup/drholdaw/tmp.8823/sens.20120318.000000/
cd /discover/nobackup/drholdaw/tmp.31275/sens.20120318.000000/

lon = ncread('fvpert.eta.nc4','lon');
lat = ncread('fvpert.eta.nc4','lat');

U_dry  = ncread('fvpert.eta.dry.20120317_06z.nc4','u');
V_dry  = ncread('fvpert.eta.dry.20120317_06z.nc4','v');
TV_dry = ncread('fvpert.eta.dry.20120317_06z.nc4','tv');
DP_dry = ncread('fvpert.eta.dry.20120317_06z.nc4','delp');
QV_dry = ncread('fvpert.eta.dry.20120317_06z.nc4','sphu');
QL_dry = ncread('fvpert.eta.dry.20120317_06z.nc4','qltot');
QI_dry = ncread('fvpert.eta.dry.20120317_06z.nc4','qitot');
O3_dry = ncread('fvpert.eta.dry.20120317_06z.nc4','ozone');

% cd /discover/nobackup/drholdaw/tmp.8823/sens.20120318.000000/
cd /discover/nobackup/drholdaw/tmp.31275/sens.20120318.000000/

U_moi  = ncread('fvpert.eta.moi.1.20120317_06z.nc4','u');
V_moi  = ncread('fvpert.eta.moi.1.20120317_06z.nc4','v');
TV_moi = ncread('fvpert.eta.moi.1.20120317_06z.nc4','tv');
DP_moi = ncread('fvpert.eta.moi.1.20120317_06z.nc4','delp');
QV_moi = ncread('fvpert.eta.moi.1.20120317_06z.nc4','sphu');
QL_moi = ncread('fvpert.eta.moi.1.20120317_06z.nc4','qltot');
QI_moi = ncread('fvpert.eta.moi.1.20120317_06z.nc4','qitot');
O3_moi = ncread('fvpert.eta.moi.1.20120317_06z.nc4','ozone');

cd /discover/nobackup/drholdaw/tmp.8823/nlm_runs/free/

U_nlm_free  = ncread('iau590.prog.eta.20120317_06z.nc4','u');
V_nlm_free  = ncread('iau590.prog.eta.20120317_06z.nc4','v');   
TV_nlm_free = ncread('iau590.prog.eta.20120317_06z.nc4','tv');
DP_nlm_free = ncread('iau590.prog.eta.20120317_06z.nc4','delp');
QV_nlm_free = ncread('iau590.prog.eta.20120317_06z.nc4','sphu');
QL_nlm_free = ncread('iau590.prog.eta.20120317_06z.nc4','qltot');
QI_nlm_free = ncread('iau590.prog.eta.20120317_06z.nc4','qitot');
O3_nlm_free = ncread('iau590.prog.eta.20120317_06z.nc4','ozone');

cd /discover/nobackup/drholdaw/tmp.8823/nlm_runs/del_1.0/

U_nlm_del_0_2  = ncread('iau590.prog.eta.20120317_06z.nc4','u');
V_nlm_del_0_2  = ncread('iau590.prog.eta.20120317_06z.nc4','v');   
TV_nlm_del_0_2 = ncread('iau590.prog.eta.20120317_06z.nc4','tv');
DP_nlm_del_0_2 = ncread('iau590.prog.eta.20120317_06z.nc4','delp');
QV_nlm_del_0_2 = ncread('iau590.prog.eta.20120317_06z.nc4','sphu');
QL_nlm_del_0_2 = ncread('iau590.prog.eta.20120317_06z.nc4','qltot');
QI_nlm_del_0_2 = ncread('iau590.prog.eta.20120317_06z.nc4','qitot');
O3_nlm_del_0_2 = ncread('iau590.prog.eta.20120317_06z.nc4','ozone');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/

U_nlm = U_nlm_del_0_2 - U_nlm_free;
V_nlm = V_nlm_del_0_2 - V_nlm_free;
TV_nlm = TV_nlm_del_0_2 - TV_nlm_free;
DP_nlm = DP_nlm_del_0_2 - DP_nlm_free;
QV_nlm = QV_nlm_del_0_2 - QV_nlm_free;
QL_nlm = QL_nlm_del_0_2 - QL_nlm_free;
QI_nlm = QI_nlm_del_0_2 - QI_nlm_free;
O3_nlm = O3_nlm_del_0_2 - O3_nlm_free;


U_moi_error  = U_moi - U_nlm;
V_moi_error  = V_moi - V_nlm;
TV_moi_error = TV_moi - TV_nlm;
DP_moi_error = DP_moi - DP_nlm;
QV_moi_error = QV_moi - QV_nlm;

U_dry_error  = U_dry - U_nlm;
V_dry_error  = V_dry - V_nlm;
TV_dry_error = TV_dry - TV_nlm;
DP_dry_error = DP_dry - DP_nlm;
QV_dry_error = QV_dry - QV_nlm;

plotlevel = 50;
line_wid = 1;
fontsize = 10;

plot_fields = 0;
plot_diffrs = 1;

Umax = max(max([U_dry(:,:,plotlevel); U_moi(:,:,plotlevel); U_nlm(:,:,plotlevel)]));
Umin = min(min([U_dry(:,:,plotlevel); U_moi(:,:,plotlevel); U_nlm(:,:,plotlevel)]));
Vmax = max(max([V_dry(:,:,plotlevel); V_moi(:,:,plotlevel); V_nlm(:,:,plotlevel)]));
Vmin = min(min([V_dry(:,:,plotlevel); V_moi(:,:,plotlevel); V_nlm(:,:,plotlevel)]));
TVmax = max(max([TV_dry(:,:,plotlevel); TV_moi(:,:,plotlevel); TV_nlm(:,:,plotlevel)]));
TVmin = min(min([TV_dry(:,:,plotlevel); TV_moi(:,:,plotlevel); TV_nlm(:,:,plotlevel)]));
DPmax = max(max([DP_dry(:,:,plotlevel); DP_moi(:,:,plotlevel); DP_nlm(:,:,plotlevel)]));
DPmin = min(min([DP_dry(:,:,plotlevel); DP_moi(:,:,plotlevel); DP_nlm(:,:,plotlevel)]));
QVmax = max(max([QV_dry(:,:,plotlevel); QV_moi(:,:,plotlevel); QV_nlm(:,:,plotlevel)]));
QVmin = min(min([QV_dry(:,:,plotlevel); QV_moi(:,:,plotlevel); QV_nlm(:,:,plotlevel)]));
QLmax = max(max([QL_dry(:,:,plotlevel); QL_moi(:,:,plotlevel); QL_nlm(:,:,plotlevel)]));
QLmin = min(min([QL_dry(:,:,plotlevel); QL_moi(:,:,plotlevel); QL_nlm(:,:,plotlevel)]));
QImax = max(max([QI_dry(:,:,plotlevel); QI_moi(:,:,plotlevel); QI_nlm(:,:,plotlevel)]));
QImin = min(min([QI_dry(:,:,plotlevel); QI_moi(:,:,plotlevel); QI_nlm(:,:,plotlevel)]));
O3max = max(max([O3_dry(:,:,plotlevel); O3_moi(:,:,plotlevel); O3_nlm(:,:,plotlevel)]));
O3min = min(min([O3_dry(:,:,plotlevel); O3_moi(:,:,plotlevel); O3_nlm(:,:,plotlevel)]));

%Make center of colormap white
cmap = colormap;
cmap(32,:) = [1 1 1];
cmap(33,:) = [1 1 1];
cmap(31,:) = [1 1 1];
cmap(34,:) = [1 1 1];
cmap(30,:) = [1 1 1];
cmap(35,:) = [1 1 1];
cmap(29,:) = [1 1 1];
cmap(36,:) = [1 1 1];
% cmap(28,:) = [1 1 1];
% cmap(37,:) = [1 1 1];
% cmap(27,:) = [1 1 1];
% cmap(38,:) = [1 1 1];
close

grey = 0.4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   PLOT THE ACTUAL FIELDS   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_fields == 1

    figure
    subplot(3,1,1)
    [C,h] = contourf(lon,lat,U_nlm(:,:,plotlevel)','LineStyle','none');
    caxis([Umin Umax])
    cint = get(h,'LevelStep');
    colorbar
    hold on; load topo.mat;
    topo_grenwich(1:180,181:360) = topo(1:180,1:180);
    topo_grenwich(1:180,1:180) = topo(1:180,181:360);
    contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);

    subplot(3,1,2)
    [C,h] = contourf(lon,lat,U_dry(:,:,plotlevel)','LineStyle','none','LevelStep',cint);
    caxis([Umin Umax])
    colorbar
    hold on; load topo.mat;
    topo_grenwich(1:180,181:360) = topo(1:180,1:180);
    topo_grenwich(1:180,1:180) = topo(1:180,181:360);
    contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);

    subplot(3,1,3)
    [C,h] = contourf(lon,lat,U_moi(:,:,plotlevel)','LineStyle','none','LevelStep',cint);
    caxis([Umin Umax])
    colorbar
    title('Q_{tot} \prime')
    hold on; load topo.mat;
    topo_grenwich(1:180,181:360) = topo(1:180,1:180);
    topo_grenwich(1:180,1:180) = topo(1:180,181:360);
    contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);

    figure
    subplot(3,1,1)
    [C,h] = contourf(lon,lat,V_nlm(:,:,plotlevel)','LineStyle','none');
    caxis([Vmin Vmax])
    cint = get(h,'LevelStep');
    colorbar
    hold on; load topo.mat;
    topo_grenwich(1:180,181:360) = topo(1:180,1:180);
    topo_grenwich(1:180,1:180) = topo(1:180,181:360);
    contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);

    subplot(3,1,2)
    [C,h] = contourf(lon,lat,V_dry(:,:,plotlevel)','LineStyle','none','LevelStep',cint);
    caxis([Vmin Vmax])
    colorbar
    hold on; load topo.mat;
    topo_grenwich(1:180,181:360) = topo(1:180,1:180);
    topo_grenwich(1:180,1:180) = topo(1:180,181:360);
    contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);

    subplot(3,1,3)
    [C,h] = contourf(lon,lat,V_moi(:,:,plotlevel)','LineStyle','none','LevelStep',cint);
    caxis([Vmin Vmax])
    colorbar
    title('Q_{tot} \prime')
    hold on; load topo.mat;
    topo_grenwich(1:180,181:360) = topo(1:180,1:180);
    topo_grenwich(1:180,1:180) = topo(1:180,181:360);
    contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);

    figure
    subplot(3,1,1)
    [C,h] = contourf(lon,lat,TV_nlm(:,:,plotlevel)','LineStyle','none');
    caxis([TVmin TVmax])
    cint = get(h,'LevelStep');
    colorbar
    hold on; load topo.mat;
    topo_grenwich(1:180,181:360) = topo(1:180,1:180);
    topo_grenwich(1:180,1:180) = topo(1:180,181:360);
    contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);

    subplot(3,1,2)
    [C,h] = contourf(lon,lat,TV_dry(:,:,plotlevel)','LineStyle','none','LevelStep',cint);
    caxis([TVmin TVmax])
    colorbar
    hold on; load topo.mat;
    topo_grenwich(1:180,181:360) = topo(1:180,1:180);
    topo_grenwich(1:180,1:180) = topo(1:180,181:360);
    contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);

    subplot(3,1,3)
    [C,h] = contourf(lon,lat,TV_moi(:,:,plotlevel)','LineStyle','none','LevelStep',cint);
    caxis([TVmin TVmax])
    colorbar
    title('Q_{tot} \prime')
    hold on; load topo.mat;
    topo_grenwich(1:180,181:360) = topo(1:180,1:180);
    topo_grenwich(1:180,1:180) = topo(1:180,181:360);
    contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);

    figure
    subplot(3,1,1)
    [C,h] = contourf(lon,lat,QV_nlm(:,:,plotlevel)','LineStyle','none');
    caxis([QVmin QVmax])
    cint = get(h,'LevelStep');
    colorbar
    title('Q_{tot} \prime')
    hold on; load topo.mat;
    topo_grenwich(1:180,181:360) = topo(1:180,1:180);
    topo_grenwich(1:180,1:180) = topo(1:180,181:360);
    contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);

    subplot(3,1,2)
    [C,h] = contourf(lon,lat,QV_dry(:,:,plotlevel)','LineStyle','none','LevelStep',cint);
    caxis([QVmin QVmax])
    colorbar
    title('Q_{tot} \prime')
    hold on; load topo.mat;
    topo_grenwich(1:180,181:360) = topo(1:180,1:180);
    topo_grenwich(1:180,1:180) = topo(1:180,181:360);
    contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);

    subplot(3,1,3)
    [C,h] = contourf(lon,lat,QV_moi(:,:,plotlevel)','LineStyle','none','LevelStep',cint);
    caxis([QVmin QVmax])
    colorbar
    title('Q_{tot} \prime')
    hold on; load topo.mat;
    topo_grenwich(1:180,181:360) = topo(1:180,1:180);
    topo_grenwich(1:180,1:180) = topo(1:180,181:360);
    contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);

    figure
    subplot(3,1,1)
    [C,h] = contourf(lon,lat,DP_nlm(:,:,plotlevel)','LineStyle','none');
    caxis([DPmin DPmax])
    cint = get(h,'LevelStep');
    colorbar
    hold on; load topo.mat;
    topo_grenwich(1:180,181:360) = topo(1:180,1:180);
    topo_grenwich(1:180,1:180) = topo(1:180,181:360);
    contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);

    subplot(3,1,2)
    [C,h] = contourf(lon,lat,DP_dry(:,:,plotlevel)','LineStyle','none','LevelStep',cint);
    caxis([DPmin DPmax])
    colorbar
    hold on; load topo.mat;
    topo_grenwich(1:180,181:360) = topo(1:180,1:180);
    topo_grenwich(1:180,1:180) = topo(1:180,181:360);
    contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);

    subplot(3,1,3)
    [C,h] = contourf(lon,lat,DP_moi(:,:,plotlevel)','LineStyle','none','LevelStep',cint);
    caxis([DPmin DPmax])
    colorbar
    title('Q_{tot} \prime')
    hold on; load topo.mat;
    topo_grenwich(1:180,181:360) = topo(1:180,1:180);
    topo_grenwich(1:180,1:180) = topo(1:180,181:360);
    contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);

end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    PLOT THE DIFFERENCES    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_diffrs == 1

    dry_tlm_U_diff = U_nlm(:,:,plotlevel)'-U_dry(:,:,plotlevel)';
    moi_tlm_U_diff = U_nlm(:,:,plotlevel)'-U_moi(:,:,plotlevel)';
    dry_tlm_V_diff = V_nlm(:,:,plotlevel)'-V_dry(:,:,plotlevel)';
    moi_tlm_V_diff = V_nlm(:,:,plotlevel)'-V_moi(:,:,plotlevel)';
    dry_tlm_TV_diff = TV_nlm(:,:,plotlevel)'-TV_dry(:,:,plotlevel)';
    moi_tlm_TV_diff = TV_nlm(:,:,plotlevel)'-TV_moi(:,:,plotlevel)';
    dry_tlm_QV_diff = QV_nlm(:,:,plotlevel)'-QV_dry(:,:,plotlevel)';
    moi_tlm_QV_diff = QV_nlm(:,:,plotlevel)'-QV_moi(:,:,plotlevel)';
    dry_tlm_DP_diff = DP_nlm(:,:,plotlevel)'-DP_dry(:,:,plotlevel)';
    moi_tlm_DP_diff = DP_nlm(:,:,plotlevel)'-DP_moi(:,:,plotlevel)';

    Udiffmax = max(max([dry_tlm_U_diff; moi_tlm_U_diff]));
    Udiffmin = min(min([dry_tlm_U_diff; moi_tlm_U_diff]));
    Vdiffmax = max(max([dry_tlm_V_diff; moi_tlm_V_diff]));
    Vdiffmin = min(min([dry_tlm_V_diff; moi_tlm_V_diff]));
    TVdiffmax = max(max([dry_tlm_TV_diff; moi_tlm_TV_diff]));
    TVdiffmin = min(min([dry_tlm_TV_diff; moi_tlm_TV_diff]));
    QVdiffmax = max(max([dry_tlm_QV_diff; moi_tlm_QV_diff]));
    QVdiffmin = min(min([dry_tlm_QV_diff; moi_tlm_QV_diff]));
    DPdiffmax = max(max([dry_tlm_DP_diff; moi_tlm_DP_diff]));
    DPdiffmin = min(min([dry_tlm_DP_diff; moi_tlm_DP_diff]));

    
    QVdiffmaxabs = max(max([QVdiffmax; QVdiffmin]));
    cint_q = 2*QVdiffmaxabs/20;
    
    
%     figure
%     subplot(2,1,1)
%     [C,h] = contourf(lon,lat,dry_tlm_U_diff,'LineStyle','none');
%     cint = get(h,'LevelStep');
%     caxis([Udiffmin Udiffmax]);
%     colorbar
%     title('Dry TLM vs Nonlinear Model','FontSize',fontsize,'FontName','TimesNewRoman')
%     hold on; load topo.mat;
%     topo_grenwich(1:180,181:360) = topo(1:180,1:180);
%     topo_grenwich(1:180,1:180) = topo(1:180,181:360);
%     contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
%     axis equal; box on; xlim([-180 180]); ylim([-90 90])
%     set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);
%     set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% 
%     subplot(2,1,2)
%     [C,h] = contourf(lon,lat,moi_tlm_U_diff,'LineStyle','none','LevelStep',cint);
%     % cint = get(h,'LevelStep')
%     caxis([Udiffmin Udiffmax])
%     colorbar
%     title('Moist TLM vs Nonlinear Model','FontSize',fontsize,'FontName','TimesNewRoman')
%     hold on; load topo.mat;
%     topo_grenwich(1:180,181:360) = topo(1:180,1:180);
%     topo_grenwich(1:180,1:180) = topo(1:180,181:360);
%     contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
%     axis equal; box on; xlim([-180 180]); ylim([-90 90])
%     set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);
%     set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
%     
%     figure
%     subplot(2,1,1)
%     [C,h] = contourf(lon,lat,dry_tlm_V_diff,'LineStyle','none');
%     cint = get(h,'LevelStep');
%     caxis([Vdiffmin Vdiffmax]);
%     colorbar
%     title('Dry TLM vs Nonlinear Model','FontSize',fontsize,'FontName','TimesNewRoman')
%     hold on; load topo.mat;
%     topo_grenwich(1:180,181:360) = topo(1:180,1:180);
%     topo_grenwich(1:180,1:180) = topo(1:180,181:360);
%     contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
%     axis equal; box on; xlim([-180 180]); ylim([-90 90])
%     set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);
%     set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% 
%     subplot(2,1,2)
%     [C,h] = contourf(lon,lat,moi_tlm_V_diff,'LineStyle','none','LevelStep',cint);
%     % cint = get(h,'LevelStep')
%     caxis([Vdiffmin Vdiffmax])
%     colorbar
%     title('Moist TLM vs Nonlinear Model','FontSize',fontsize,'FontName','TimesNewRoman')
%     hold on; load topo.mat;
%     topo_grenwich(1:180,181:360) = topo(1:180,1:180);
%     topo_grenwich(1:180,1:180) = topo(1:180,181:360);
%     contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
%     axis equal; box on; xlim([-180 180]); ylim([-90 90])
%     set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);
%     set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% 
%     figure
%     subplot(2,1,1)
%     [C,h] = contourf(lon,lat,dry_tlm_TV_diff,'LineStyle','none');
%     cint = get(h,'LevelStep');
%     caxis([TVdiffmin TVdiffmax]);
%     colorbar
%     title('Dry TLM vs Nonlinear Model','FontSize',fontsize,'FontName','TimesNewRoman')
%     hold on; load topo.mat;
%     topo_grenwich(1:180,181:360) = topo(1:180,1:180);
%     topo_grenwich(1:180,1:180) = topo(1:180,181:360);
%     contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
%     axis equal; box on; xlim([-180 180]); ylim([-90 90])
%     set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);
%     set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% 
%     subplot(2,1,2)
%     [C,h] = contourf(lon,lat,moi_tlm_TV_diff,'LineStyle','none','LevelStep',cint);
%     % cint = get(h,'LevelStep')
%     caxis([TVdiffmin TVdiffmax])
%     colorbar
%     title('Moist TLM vs Nonlinear Model','FontSize',fontsize,'FontName','TimesNewRoman')
%     hold on; load topo.mat;
%     topo_grenwich(1:180,181:360) = topo(1:180,1:180);
%     topo_grenwich(1:180,1:180) = topo(1:180,181:360);
%     contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
%     axis equal; box on; xlim([-180 180]); ylim([-90 90])
%     set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);
%     set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
    
    figure
    subplot(2,1,1)
    [C,h] = contourf(lon,lat,dry_tlm_QV_diff,'LineStyle','none');
    colormap(cmap)
    set(h,'LevelStep',cint_q);
    caxis([-QVdiffmaxabs QVdiffmaxabs]);
    colorbar
    hold on; 
    plot(coast_lon,coast_lat,'Color',[grey grey grey])
    title('Dry TLM vs Nonlinear Model','FontSize',fontsize,'FontName','TimesNewRoman')
    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);
    set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

    subplot(2,1,2)
    [C,h] = contourf(lon,lat,moi_tlm_QV_diff,'LineStyle','none');
    colormap(cmap)
    set(h,'LevelStep',cint_q);
    caxis([-QVdiffmaxabs QVdiffmaxabs]);
    colorbar
    hold on; 
    plot(coast_lon,coast_lat,'Color',[grey grey grey])
    title('Moist TLM vs Nonlinear Model','FontSize',fontsize,'FontName','TimesNewRoman')
    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);
    set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
    
%     figure
%     subplot(2,1,1)
%     [C,h] = contourf(lon,lat,dry_tlm_DP_diff,'LineStyle','none');
%     cint = get(h,'LevelStep');
%     caxis([DPdiffmin DPdiffmax]);
%     colorbar
%     title('Dry TLM vs Nonlinear Model','FontSize',fontsize,'FontName','TimesNewRoman')
%     hold on; load topo.mat;
%     topo_grenwich(1:180,181:360) = topo(1:180,1:180);
%     topo_grenwich(1:180,1:180) = topo(1:180,181:360);
%     contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
%     axis equal; box on; xlim([-180 180]); ylim([-90 90])
%     set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);
%     set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% 
%     subplot(2,1,2)
%     [C,h] = contourf(lon,lat,moi_tlm_DP_diff,'LineStyle','none','LevelStep',cint);
%     % cint = get(h,'LevelStep')
%     caxis([DPdiffmin DPdiffmax])
%     colorbar
%     title('Moist TLM vs Nonlinear Model','FontSize',fontsize,'FontName','TimesNewRoman')
%     hold on; load topo.mat;
%     topo_grenwich(1:180,181:360) = topo(1:180,1:180);
%     topo_grenwich(1:180,1:180) = topo(1:180,181:360);
%     contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
%     axis equal; box on; xlim([-180 180]); ylim([-90 90])
%     set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);
%     set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')


end


