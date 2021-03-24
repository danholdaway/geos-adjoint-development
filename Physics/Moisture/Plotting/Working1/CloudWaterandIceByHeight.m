close all
clear
clc

cd /discover/nobackup/drholdaw/atmp.22292/prog/prog_free/

file = ['x0011dh_a.prog.eta.20130117_0300z.nc4'];

lon = ncread(file,'lon');
lat = ncread(file,'lat');
lev = ncread(file,'lev');

qi_free = ncread(file,'qitot');
ql_free = ncread(file,'qltot');

cd /home/drholdaw/Lin_Moist_Physics/Inputs/
pref

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

lm = length(lev);

qi_tot = zeros(lm,1);
ql_tot = zeros(lm,1);

qi_nh = zeros(lm,1);
ql_nh = zeros(lm,1);

qi_sh = zeros(lm,1);
ql_sh = zeros(lm,1);

for k = 1:lm
    
    qi_tot(k) = mean(mean(qi_free(:,:,k)));
    ql_tot(k) = mean(mean(ql_free(:,:,k)));
    
    qi_sh(k) = mean(mean(qi_free(:,1:180,k)));
    ql_sh(k) = mean(mean(ql_free(:,1:180,k)));
    
    qi_nh(k) = mean(mean(qi_free(:,182:361,k)));
    ql_nh(k) = mean(mean(ql_free(:,182:361,k)));
    
end

grey = 0.75;
fontsize = 11;
line_wid_cont = 1.0;
line_wid_det = 0.6;

figure
set(gcf,'position',[788   486   491   433])

subplot(1,2,1)
plot(ql_tot,p_ref(2:lm+1),'k')
hold on
plot(ql_nh,p_ref(2:lm+1),'k--')
plot(ql_sh,p_ref(2:lm+1),'Color',[grey grey grey])
set(gca,'Ydir','reverse')
ylim([50 1000])
legend('Global','N. Hemisphere','S. Hemisphere')

set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
title('(a)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('q_l (kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Pressure (hPa)','FontSize',fontsize,'FontName','TimesNewRoman')
box on

subplot(1,2,2)
plot(qi_tot,p_ref(2:lm+1),'k')
hold on
plot(qi_nh,p_ref(2:lm+1),'k--')
plot(qi_sh,p_ref(2:lm+1),'Color',[grey grey grey])
set(gca,'Ydir','reverse')
set(gca,'YAxisLocation','right')
ylim([50 1000])

set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
title('(b)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('q_i (kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Pressure (hPa)','FontSize',fontsize,'FontName','TimesNewRoman')
box on

pos = get(gca,'position')
set(gca,'position',[0.9*pos(1) pos(2) pos(3) pos(4)])


figure
set(gcf,'position',[788   486   491   433])

subplot(1,2,1)
plot(ql_tot,1:lm,'k')
hold on
plot(ql_nh,1:lm,'k--')
plot(ql_sh,1:lm,'Color',[grey grey grey])
set(gca,'Ydir','reverse')
ylim([30 72])
legend('Global','N. Hemisphere','S. Hemisphere')

set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
title('(a)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('q_l (kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Pressure (hPa)','FontSize',fontsize,'FontName','TimesNewRoman')
box on

subplot(1,2,2)
plot(qi_tot,1:lm,'k')
hold on
plot(qi_nh,1:lm,'k--')
plot(qi_sh,1:lm,'Color',[grey grey grey])
set(gca,'Ydir','reverse')
set(gca,'YAxisLocation','right')
ylim([30 72])

set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
title('(b)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('q_i (kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Pressure (hPa)','FontSize',fontsize,'FontName','TimesNewRoman')
box on

pos = get(gca,'position')
set(gca,'position',[0.9*pos(1) pos(2) pos(3) pos(4)])
