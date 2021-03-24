close all
close all
clear
clc


J = ncread('JACOBIAN_DEEP_SERIES.nc4','J');
kcbl = ncread('JACOBIAN_DEEP_SERIES.nc4','KCBL');
top = ncread('JACOBIAN_DEEP_SERIES.nc4','TOP');
depth = ncread('JACOBIAN_DEEP_SERIES.nc4','DEPTH');

DUDt = ncread('DIAG.nc4','DUDt');
DVDt = ncread('DIAG.nc4','DVDt');
DTHDt = ncread('DIAG.nc4','DTHDt');
DQDt = ncread('DIAG.nc4','DQDt');
PRC3 = ncread('DIAG.nc4','PRC3');
UPFR = ncread('DIAG.nc4','UPFR');
ENTP = ncread('DIAG.nc4','ENTP');
MFCU = ncread('DIAG.nc4','MFCU');

cd /discover/nobackup/drholdaw/inputs_forRAS/high_vs_low/

U = ncread('high_vs_low_14491.geosgcm_rasinputs3d.20100702_1200z.nc4','U_preras');
V = ncread('high_vs_low_14491.geosgcm_rasinputs3d.20100702_1200z.nc4','V_preras');
TH = ncread('high_vs_low_14491.geosgcm_rasinputs3d.20100702_1200z.nc4','TH_preras');
Q = ncread('high_vs_low_14491.geosgcm_rasinputs3d.20100702_1200z.nc4','Q_preras');

cd /home/drholdaw/Lin_Moist_Physics/RAS_Jacobian_Paper/mfiles/

% pref_cont

pref = [  1.000000000000000E-002  2.000000000000000E-002  3.270000000000000E-002 ...
  4.758501000000000E-002  6.600001000000000E-002  8.934500999999999E-002 ...
  0.119703000000000       0.159495000000000       0.211348990000000      ...
  0.278526000000000       0.365041010000000       0.475806010000000      ...
  0.616779100000000       0.795134120000000        1.01944000000000      ...
   1.30050995000000        1.65078995000000        2.08496994000000      ...
   2.62020996000000        3.27643005000000        4.07657104000000      ...
   5.04680115000000        6.21680115000000        7.61984070000000      ...
   9.29294128000000        11.2768994100000        13.6433996600000      ...
   16.4570996100000        19.7916003400000        23.7304003900000      ...
   28.3678002900000        33.8100000000000        40.1754101600000      ...
   47.6439111300000        56.3879101600000        66.6034082000000      ...
   78.5123095700000        92.3657226600000        108.662998050000      ...
   127.837001950000        150.392998050000        176.930000000000      ...
   207.642583808750        242.824924028906        283.419286250625      ...
   329.185379803281        363.565729625000        397.836113664219      ...
   431.940720416953        465.952285957422        499.867555766953      ...
   533.743876319297        567.522190400938        601.300435201719      ...
   635.015439200000        657.472713397891        679.927450581406      ...
   702.378753241406        724.797384625547        747.216522878828      ...
   765.152362603281        778.596047543828        792.035999691406      ...
   805.472753127734        818.909516334062        832.347306449531      ...
   845.780018229453        859.208150490938        872.636185102422      ...
   886.063891034766        899.487238762969        912.909638117734      ...
   926.809140625000 ];

prefh = 0.5*(pref(1:end-1)+pref(2:end));

% prefh = prefh+75;

% cd /discover/nobackup/drholdaw/tmp.13503/sens.20111115.000000/
% U = ncread('fsens.eta.111114_06z_nom.nc4','U');
% V = ncread('fsens.eta.111114_06z_nom.nc4','V');
% TV = ncread('fsens.eta.111114_06z_nom.nc4','TV');
% QV = ncread('fsens.eta.111114_06z_nom.nc4','QV');
% cd /home/drholdaw/Lin_Moist_Physics/RAS_Jacobian_Paper/mfiles/

Ni = 116;
Nj = 61;
Nk = 72;

Nk_vec = 1:72;

vec1 = [top kcbl];
vec2 = ones(1,2);

U_plot = zeros(72,1);
V_plot = zeros(72,1);
TH_plot = zeros(72,1);
Q_plot = zeros(72,1);
DUDt_plot = zeros(72,1);
DVDt_plot = zeros(72,1);
DTHDt_plot = zeros(72,1);
DQDt_plot = zeros(72,1);
PRC3_plot = zeros(72,1);
UPFR_plot = zeros(72,1);
ENTP_plot = zeros(72,1);
MFCU_plot = zeros(72,1);

U_plot(:) = U(Ni,Nj,:);
V_plot(:) = V(Ni,Nj,:);
TH_plot(:) = TH(Ni,Nj,:);
Q_plot(:) = Q(Ni,Nj,:);
DUDt_plot(:) = DUDt(Ni,Nj,:);
DVDt_plot(:) = DVDt(Ni,Nj,:);
DTHDt_plot(:) = DTHDt(Ni,Nj,:);
DQDt_plot(:) = DQDt(Ni,Nj,:);
PRC3_plot(:) = PRC3(Ni,Nj,:);
UPFR_plot(:) = UPFR(Ni,Nj,:);
ENTP_plot(:) = ENTP(Ni,Nj,:);
MFCU_plot(:) = MFCU(Ni,Nj,:);

% x = zeros(4*Nk,1);
% x(0*Nk+1:1*Nk) = U(150,100,:);
% x(1*Nk+1:2*Nk) = V(150,100,:);
% x(2*Nk+1:3*Nk) = TV(150,100,:);
% x(3*Nk+1:4*Nk) = QV(150,100,:);
% J = jtimesx(J(1:4*Nk,1:4*Nk),x);

%%%%%%%%%%%%%%%%%%%
J_Auu = zeros(Nk,Nk);
J_Aut = zeros(Nk,Nk);
J_Auq = zeros(Nk,Nk);

J_Avv = zeros(Nk,Nk);
J_Avt = zeros(Nk,Nk);
J_Avq = zeros(Nk,Nk);

J_Ht  = zeros(Nk,Nk);
J_Hq  = zeros(Nk,Nk);

J_Mt  = zeros(Nk,Nk);
J_Mq  = zeros(Nk,Nk);

J_Auu(:,:) = J(1,0*Nk+1:1*Nk,0*Nk+1:1*Nk);
J_Aut(:,:) = J(1,0*Nk+1:1*Nk,2*Nk+1:3*Nk);
J_Auq(:,:) = J(1,0*Nk+1:1*Nk,3*Nk+1:4*Nk);

J_Avv(:,:) = J(1,1*Nk+1:2*Nk,1*Nk+1:2*Nk);
J_Avt(:,:) = J(1,1*Nk+1:2*Nk,2*Nk+1:3*Nk);
J_Avq(:,:) = J(1,1*Nk+1:2*Nk,3*Nk+1:4*Nk);

J_Ht(:,:)  = J(1,2*Nk+1:3*Nk,2*Nk+1:3*Nk);
J_Hq(:,:)  = J(1,2*Nk+1:3*Nk,3*Nk+1:4*Nk);

J_Mt(:,:)  = J(1,3*Nk+1:4*Nk,2*Nk+1:3*Nk);
J_Mq(:,:)  = J(1,3*Nk+1:4*Nk,3*Nk+1:4*Nk);
%%%%%%%%%%%%%%%%%%%

% percent_cut = 0.01;
% 
% %REMOVE SMALL ELEMENT OF J_Ht
% J_Ht = remove0(J_Ht,percent_cut);
% J_Hq = remove0(J_Hq,percent_cut);
% J_Mt = remove0(J_Mt,percent_cut);
% J_Mq = remove0(J_Mq,percent_cut);
% 
% J_Auu = remove0(J_Auu,percent_cut);
% J_Aut = remove0(J_Aut,percent_cut);
% J_Auq = remove0(J_Auq,percent_cut);
% 
% J_Avv = remove0(J_Avv,percent_cut);
% J_Avt = remove0(J_Avt,percent_cut);
% J_Avq = remove0(J_Avq,percent_cut);



line_wid = 1;
fontsize = 14;

topm = 7;


figure



plot(DTHDt_plot(top-topm:end)*10e-5,top-topm:Nk)
set(gca,'YDir','reverse')
ylim([top-topm Nk])
hold on
plot(J_Mt(62,top-topm:end),Nk_vec(top-topm:end),'r');
set(gca,'YDir','reverse')
ylim([top-topm Nk])

figure

plot((DQDt_plot(top-topm:end)),top-topm:Nk)
set(gca,'YDir','reverse')
ylim([top-topm Nk])


figure
plot(TH_plot(top-topm:end),top-topm:Nk)
set(gca,'YDir','reverse')
ylim([top-topm Nk])

figure
plot(Q_plot(top-topm:end),top-topm:Nk)
set(gca,'YDir','reverse')
ylim([top-topm Nk])

fontsize = 8;





figure

[C,h] = contour(Nk_vec(top-topm:end),Nk_vec(top-topm:end),J_Mt(top-topm:end,top-topm:end)');
climmax = max(abs(caxis)); caxis([-climmax climmax])
set(h,'LineWidth',line_wid);
set(gca,'YDir','reverse')
title('(c) \partial Q /\partial \theta (kgkg^{-1}s^{-1}K^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Output Index (Q)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Perturbation Index (\theta)','FontSize',fontsize,'FontName','TimesNewRoman')
colormap jet
colorbar
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

    levlist = get(h,'LevelList');
    indzero = find(levlist==0)
    levlist(indzero) = [];
    set(h,'LevelList',levlist)

    
figure
plot(J_Mt(62,top-topm:end),Nk_vec(top-topm:end));
set(gca,'YDir','reverse')
ylim([top-topm Nk])


