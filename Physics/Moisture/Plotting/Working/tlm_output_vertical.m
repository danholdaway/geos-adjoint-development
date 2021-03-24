close all
clear
clc

cd /discover/nobackup/drholdaw/tmp.8823/sens.20120318.000000/

lon = ncread('fvpert.eta.nc4','lon');
lat = ncread('fvpert.eta.nc4','lat');
lev = ncread('fvpert.eta.nc4','lev');

U_dry  = ncread('fvpert.eta.dry.20120317_12z.nc4','u');
V_dry  = ncread('fvpert.eta.dry.20120317_12z.nc4','v');
TV_dry = ncread('fvpert.eta.dry.20120317_12z.nc4','tv');
DP_dry = ncread('fvpert.eta.dry.20120317_12z.nc4','delp');
QV_dry = ncread('fvpert.eta.dry.20120317_12z.nc4','sphu');
QL_dry = ncread('fvpert.eta.dry.20120317_12z.nc4','qltot');
QI_dry = ncread('fvpert.eta.dry.20120317_12z.nc4','qitot');
O3_dry = ncread('fvpert.eta.dry.20120317_12z.nc4','ozone');

U_moi  = ncread('fvpert.eta.moi.20120317_12z.nc4','u');
V_moi  = ncread('fvpert.eta.moi.20120317_12z.nc4','v');
TV_moi = ncread('fvpert.eta.moi.20120317_12z.nc4','tv');
DP_moi = ncread('fvpert.eta.moi.20120317_12z.nc4','delp');
QV_moi = ncread('fvpert.eta.moi.20120317_12z.nc4','sphu');
QL_moi = ncread('fvpert.eta.moi.20120317_12z.nc4','qltot');
QI_moi = ncread('fvpert.eta.moi.20120317_12z.nc4','qitot');
O3_moi = ncread('fvpert.eta.moi.20120317_12z.nc4','ozone');

cd /discover/nobackup/drholdaw/tmp.8823/nlm_runs/free/

U_nlm_free  = ncread('iau590.prog.eta.20120317_12z.nc4','u');
V_nlm_free  = ncread('iau590.prog.eta.20120317_12z.nc4','v');   
TV_nlm_free = ncread('iau590.prog.eta.20120317_12z.nc4','tv');
DP_nlm_free = ncread('iau590.prog.eta.20120317_12z.nc4','delp');
QV_nlm_free = ncread('iau590.prog.eta.20120317_12z.nc4','sphu');
QL_nlm_free = ncread('iau590.prog.eta.20120317_12z.nc4','qltot');
QI_nlm_free = ncread('iau590.prog.eta.20120317_12z.nc4','qitot');
O3_nlm_free = ncread('iau590.prog.eta.20120317_12z.nc4','ozone');

cd /discover/nobackup/drholdaw/tmp.8823/nlm_runs/del_0.2/

U_nlm_del_0_2  = ncread('iau590.prog.eta.20120317_12z.nc4','u');
V_nlm_del_0_2  = ncread('iau590.prog.eta.20120317_12z.nc4','v');   
TV_nlm_del_0_2 = ncread('iau590.prog.eta.20120317_12z.nc4','tv');
DP_nlm_del_0_2 = ncread('iau590.prog.eta.20120317_12z.nc4','delp');
QV_nlm_del_0_2 = ncread('iau590.prog.eta.20120317_12z.nc4','sphu');
QL_nlm_del_0_2 = ncread('iau590.prog.eta.20120317_12z.nc4','qltot');
QI_nlm_del_0_2 = ncread('iau590.prog.eta.20120317_12z.nc4','qitot');
O3_nlm_del_0_2 = ncread('iau590.prog.eta.20120317_12z.nc4','ozone');

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


lat = 80;
levs = 50:72;

Umax = max(max([U_dry(:,lat,levs); U_moi(:,lat,levs); U_nlm(:,lat,levs)]));
Umin = min(min([U_dry(:,lat,levs); U_moi(:,lat,levs); U_nlm(:,lat,levs)]));
Vmax = max(max([V_dry(:,lat,levs); V_moi(:,lat,levs); V_nlm(:,lat,levs)]));
Vmin = min(min([V_dry(:,lat,levs); V_moi(:,lat,levs); V_nlm(:,lat,levs)]));
TVmax = max(max([TV_dry(:,lat,levs); TV_moi(:,lat,levs); TV_nlm(:,lat,levs)]));
TVmin = min(min([TV_dry(:,lat,levs); TV_moi(:,lat,levs); TV_nlm(:,lat,levs)]));
DPmax = max(max([DP_dry(:,lat,levs); DP_moi(:,lat,levs); DP_nlm(:,lat,levs)]));
DPmin = min(min([DP_dry(:,lat,levs); DP_moi(:,lat,levs); DP_nlm(:,lat,levs)]));
QVmax = max(max([QV_dry(:,lat,levs); QV_moi(:,lat,levs); QV_nlm(:,lat,levs)]));
QVmin = min(min([QV_dry(:,lat,levs); QV_moi(:,lat,levs); QV_nlm(:,lat,levs)]));
QLmax = max(max([QL_dry(:,lat,levs); QL_moi(:,lat,levs); QL_nlm(:,lat,levs)]));
QLmin = min(min([QL_dry(:,lat,levs); QL_moi(:,lat,levs); QL_nlm(:,lat,levs)]));
QImax = max(max([QI_dry(:,lat,levs); QI_moi(:,lat,levs); QI_nlm(:,lat,levs)]));
QImin = min(min([QI_dry(:,lat,levs); QI_moi(:,lat,levs); QI_nlm(:,lat,levs)]));
O3max = max(max([O3_dry(:,lat,levs); O3_moi(:,lat,levs); O3_nlm(:,lat,levs)]));
O3min = min(min([O3_dry(:,lat,levs); O3_moi(:,lat,levs); O3_nlm(:,lat,levs)]));

TV_nlm_plot = zeros(length(lon),length(levs));
TV_nlm_plot(:,:) = TV_nlm(:,lat,levs);
TV_dry_plot = zeros(length(lon),length(levs));
TV_dry_plot(:,:) = TV_dry(:,lat,levs);
TV_moi_plot = zeros(length(lon),length(levs));
TV_moi_plot(:,:) = TV_moi(:,lat,levs);

QV_nlm_plot = zeros(length(lon),length(levs));
QV_nlm_plot(:,:) = QV_nlm(:,lat,levs);
QV_dry_plot = zeros(length(lon),length(levs));
QV_dry_plot(:,:) = QV_dry(:,lat,levs);
QV_moi_plot = zeros(length(lon),length(levs));
QV_moi_plot(:,:) = QV_moi(:,lat,levs);

figure
subplot(3,1,1)
[C,h] = contourf(lon,levs,TV_nlm_plot','LineStyle','none');
caxis([TVmin TVmax])
cint = get(h,'LevelStep');
set(gca,'YDir','reverse')
colorbar

subplot(3,1,2)
[C,h] = contourf(lon,levs,TV_dry_plot','LineStyle','none','LevelStep',cint);
caxis([TVmin TVmax])
set(gca,'YDir','reverse')
colorbar

subplot(3,1,3)
[C,h] = contourf(lon,levs,TV_moi_plot','LineStyle','none','LevelStep',cint);
caxis([TVmin TVmax])
set(gca,'YDir','reverse')
colorbar


figure
subplot(3,1,1)
[C,h] = contourf(lon,levs,QV_nlm_plot','LineStyle','none');
caxis([QVmin QVmax])
cint = get(h,'LevelStep');
set(gca,'YDir','reverse')
colorbar

subplot(3,1,2)
[C,h] = contourf(lon,levs,QV_dry_plot','LineStyle','none','LevelStep',cint);
caxis([QVmin QVmax])
set(gca,'YDir','reverse')
colorbar

subplot(3,1,3)
[C,h] = contourf(lon,levs,QV_moi_plot','LineStyle','none','LevelStep',cint);
caxis([QVmin QVmax])
set(gca,'YDir','reverse')
colorbar


%Plot the differences for the level.

dry_tlm_TV_diff = TV_nlm_plot-TV_dry_plot;
moi_tlm_TV_diff = TV_nlm_plot-TV_moi_plot;
dry_tlm_QV_diff = QV_nlm_plot-QV_dry_plot;
moi_tlm_QV_diff = QV_nlm_plot-QV_moi_plot;

TVdiffmax = max(max([dry_tlm_TV_diff; moi_tlm_TV_diff]));
TVdiffmin = min(min([dry_tlm_TV_diff; moi_tlm_TV_diff]));
QVdiffmax = max(max([dry_tlm_QV_diff; moi_tlm_QV_diff]));
QVdiffmin = min(min([dry_tlm_QV_diff; moi_tlm_QV_diff]));

figure
subplot(2,1,1)
[C,h] = contourf(lon,levs,dry_tlm_TV_diff','LineStyle','none');
cint = get(h,'LevelStep');
caxis([TVdiffmin TVdiffmax]);
set(gca,'YDir','reverse')
colorbar

subplot(2,1,2)
[C,h] = contourf(lon,levs,moi_tlm_TV_diff','LineStyle','none','LevelStep',cint);
caxis([TVdiffmin TVdiffmax])
set(gca,'YDir','reverse')
colorbar


figure
subplot(2,1,1)
[C,h] = contourf(lon,levs,dry_tlm_QV_diff','LineStyle','none');
cint = get(h,'LevelStep');
caxis([QVdiffmin QVdiffmax]);
set(gca,'YDir','reverse')
colorbar

subplot(2,1,2)
[C,h] = contourf(lon,levs,moi_tlm_QV_diff','LineStyle','none','LevelStep',cint);
caxis([QVdiffmin QVdiffmax])
set(gca,'YDir','reverse')
colorbar





% scrsz = get(0,'ScreenSize');
% figure('visible','on','Position',[0 scrsz(4) scrsz(3)/3 scrsz(4)/2])
% 
% subplot(2,4,1)
% contourf(lon,lat,U_dry(:,lat,levs)')
% % caxis([Umin Umax])
% colorbar
% title('U \prime')
% 
% subplot(2,4,2)
% contourf(lon,lat,V_dry(:,lat,levs)')
% % caxis([Vmin Vmax])
% colorbar
% title('V \prime')
% 
% subplot(2,4,3)
% contourf(lon,lat,TV_dry(:,lat,levs)')
% % caxis([TVmin TVmax])
% colorbar
% title('\theta \prime')
% 
% subplot(2,4,4)
% contourf(lon,lat,DP_dry(:,lat,levs)')
% % caxis([DPmin DPmax])
% colorbar
% title('P \prime')
% 
% subplot(2,4,5)
% contourf(lon,lat,QV_dry(:,lat,levs)')
% % caxis([QVmin QVmax])
% colorbar
% title('Q_{tot} \prime')
% 
% subplot(2,4,6)
% contourf(lon,lat,QL_dry(:,lat,levs)')
% % caxis([QLmin QLmax])
% colorbar
% title('Q_L \prime')
% 
% subplot(2,4,7)
% contourf(lon,lat,QI_dry(:,lat,levs)')
% % caxis([QImin QImax])
% colorbar
% title('Q_I \prime')
% 
% subplot(2,4,8)
% contourf(lon,lat,O3_dry(:,lat,levs)')
% % caxis([O3min O3max])
% colorbar
% title('O^3 \prime')
% 
% 
% figure('visible','on','Position',[scrsz(3)/2 scrsz(4) scrsz(3)/3 scrsz(4)/2])
% subplot(2,4,1)
% contourf(lon,lat,U_moi(:,lat,levs)')
% % caxis([Umin Umax])
% colorbar
% title('U \prime')
% 
% subplot(2,4,2)
% contourf(lon,lat,V_moi(:,lat,levs)')
% % caxis([Vmin Vmax])
% colorbar
% title('V \prime')
% 
% subplot(2,4,3)
% contourf(lon,lat,TV_moi(:,lat,levs)')
% % caxis([TVmin TVmax])
% colorbar
% title('\theta \prime')
% 
% subplot(2,4,4)
% contourf(lon,lat,DP_moi(:,lat,levs)')
% % caxis([DPmin DPmax])
% colorbar
% title('P \prime')
% 
% subplot(2,4,5)
% contourf(lon,lat,QV_moi(:,lat,levs)')
% % caxis([QVmin QVmax])
% colorbar
% title('Q_{tot} \prime')
% 
% subplot(2,4,6)
% contourf(lon,lat,QL_moi(:,lat,levs)')
% % caxis([QLmin QLmax])
% colorbar
% title('Q_L \prime')
% 
% subplot(2,4,7)
% contourf(lon,lat,QI_moi(:,lat,levs)')
% % caxis([QImin QImax])
% colorbar
% title('Q_I \prime')
% 
% subplot(2,4,8)
% contourf(lon,lat,O3_moi(:,lat,levs)')
% % caxis([O3min O3max])
% colorbar
% title('O^3 \prime')

