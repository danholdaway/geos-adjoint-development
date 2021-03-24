close all
clear
clc

load coast
coast_lat = lat; clear lat
coast_lon = long; clear long

cd /discover/nobackup/drholdaw/ExperimentData/Journal_Articles/moist_dev_mwr/storm_sensitivity/

lon = ncread('iau590a.prog.eta.20120324_00z.nc4','lon');
lat = ncread('iau590a.prog.eta.20120324_00z.nc4','lat');

SLP_1 = ncread('iau590a.tavg2d_met_x_unz.20120324_0130z.nc4','SLP');
SLP_1 = SLP_1*0.01;
SLP_2 = ncread('iau590a.tavg2d_met_x_unz.20120324_2230z.nc4','SLP');
SLP_2 = SLP_2*0.01;

PRECCON_1 = ncread('iau590a.tavg2d_met_x_unz.20120324_0130z.nc4','PRECCON');
PRECLSC_1 = ncread('iau590a.tavg2d_met_x_unz.20120324_0130z.nc4','PRECLSC');
PRECANV_1 = ncread('iau590a.tavg2d_met_x_unz.20120324_0130z.nc4','PRECANV');
PRECTOT_1 = ncread('iau590a.tavg2d_met_x_unz.20120324_0130z.nc4','PRECTOT');

PRECCON_2 = ncread('iau590a.tavg2d_met_x_unz.20120324_2230z.nc4','PRECCON');
PRECLSC_2 = ncread('iau590a.tavg2d_met_x_unz.20120324_2230z.nc4','PRECLSC');
PRECANV_2 = ncread('iau590a.tavg2d_met_x_unz.20120324_2230z.nc4','PRECANV');
PRECTOT_2 = ncread('iau590a.tavg2d_met_x_unz.20120324_2230z.nc4','PRECTOT');

U_dp_dn = ncread('fsens_dp_dn.20120323_15z+20120325_00z-20120324_00z.nc4','u');
V_dp_dn = ncread('fsens_dp_dn.20120323_15z+20120325_00z-20120324_00z.nc4','v');
T_dp_dn = ncread('fsens_dp_dn.20120323_15z+20120325_00z-20120324_00z.nc4','tv');
Q_dp_dn = ncread('fsens_dp_dn.20120323_15z+20120325_00z-20120324_00z.nc4','sphu');
P_dp_dn = ncread('fsens_dp_dn.20120323_15z+20120325_00z-20120324_00z.nc4','delp');

U_dp_mn = ncread('fsens_dp_mn.20120323_15z+20120325_00z-20120324_00z.nc4','u');
V_dp_mn = ncread('fsens_dp_mn.20120323_15z+20120325_00z-20120324_00z.nc4','v');
T_dp_mn = ncread('fsens_dp_mn.20120323_15z+20120325_00z-20120324_00z.nc4','tv');
Q_dp_mn = ncread('fsens_dp_mn.20120323_15z+20120325_00z-20120324_00z.nc4','sphu');
P_dp_mn = ncread('fsens_dp_mn.20120323_15z+20120325_00z-20120324_00z.nc4','delp');

U_mp_dn = ncread('fsens_mp_dn.20120323_15z+20120325_00z-20120324_00z.nc4','u');
V_mp_dn = ncread('fsens_mp_dn.20120323_15z+20120325_00z-20120324_00z.nc4','v');
T_mp_dn = ncread('fsens_mp_dn.20120323_15z+20120325_00z-20120324_00z.nc4','tv');
Q_mp_dn = ncread('fsens_mp_dn.20120323_15z+20120325_00z-20120324_00z.nc4','sphu');
P_mp_dn = ncread('fsens_mp_dn.20120323_15z+20120325_00z-20120324_00z.nc4','delp');

U_mp_mn = ncread('fsens_mp_mn.20120323_15z+20120325_00z-20120324_00z.nc4','u');
V_mp_mn = ncread('fsens_mp_mn.20120323_15z+20120325_00z-20120324_00z.nc4','v');
T_mp_mn = ncread('fsens_mp_mn.20120323_15z+20120325_00z-20120324_00z.nc4','tv');
Q_mp_mn = ncread('fsens_mp_mn.20120323_15z+20120325_00z-20120324_00z.nc4','sphu');
P_mp_mn = ncread('fsens_mp_mn.20120323_15z+20120325_00z-20120324_00z.nc4','delp');

QV = ncread('Jgradf_twe.eta.nc4','tv');

cd /home/drholdaw/Lin_Moist_Physics/Moist_Dev_MWR/

lati1 = 101;
lati2 = 171;
loni1 = 65;
loni2 = 145;

plotlevel = 50;

max_U_dp_dn = max(max(U_dp_dn(loni1:loni2,lati1:lati2,plotlevel)));
min_U_dp_dn = min(min(U_dp_dn(loni1:loni2,lati1:lati2,plotlevel)));
max_V_dp_dn = max(max(V_dp_dn(loni1:loni2,lati1:lati2,plotlevel)));
min_V_dp_dn = min(min(V_dp_dn(loni1:loni2,lati1:lati2,plotlevel)));
max_T_dp_dn = max(max(T_dp_dn(loni1:loni2,lati1:lati2,plotlevel)));
min_T_dp_dn = min(min(T_dp_dn(loni1:loni2,lati1:lati2,plotlevel)));
max_Q_dp_dn = max(max(Q_dp_dn(loni1:loni2,lati1:lati2,plotlevel)));
min_Q_dp_dn = min(min(Q_dp_dn(loni1:loni2,lati1:lati2,plotlevel)));
max_P_dp_dn = max(max(P_dp_dn(loni1:loni2,lati1:lati2,plotlevel)));
min_P_dp_dn = min(min(P_dp_dn(loni1:loni2,lati1:lati2,plotlevel)));

max_U_dp_mn = max(max(U_dp_mn(loni1:loni2,lati1:lati2,plotlevel)));
min_U_dp_mn = min(min(U_dp_mn(loni1:loni2,lati1:lati2,plotlevel)));
max_V_dp_mn = max(max(V_dp_mn(loni1:loni2,lati1:lati2,plotlevel)));
min_V_dp_mn = min(min(V_dp_mn(loni1:loni2,lati1:lati2,plotlevel)));
max_T_dp_mn = max(max(T_dp_mn(loni1:loni2,lati1:lati2,plotlevel)));
min_T_dp_mn = min(min(T_dp_mn(loni1:loni2,lati1:lati2,plotlevel)));
max_Q_dp_mn = max(max(Q_dp_mn(loni1:loni2,lati1:lati2,plotlevel)));
min_Q_dp_mn = min(min(Q_dp_mn(loni1:loni2,lati1:lati2,plotlevel)));
max_P_dp_mn = max(max(P_dp_mn(loni1:loni2,lati1:lati2,plotlevel)));
min_P_dp_mn = min(min(P_dp_mn(loni1:loni2,lati1:lati2,plotlevel)));

max_U_mp_dn = max(max(U_mp_dn(loni1:loni2,lati1:lati2,plotlevel)));
min_U_mp_dn = min(min(U_mp_dn(loni1:loni2,lati1:lati2,plotlevel)));
max_V_mp_dn = max(max(V_mp_dn(loni1:loni2,lati1:lati2,plotlevel)));
min_V_mp_dn = min(min(V_mp_dn(loni1:loni2,lati1:lati2,plotlevel)));
max_T_mp_dn = max(max(T_mp_dn(loni1:loni2,lati1:lati2,plotlevel)));
min_T_mp_dn = min(min(T_mp_dn(loni1:loni2,lati1:lati2,plotlevel)));
max_Q_mp_dn = max(max(Q_mp_dn(loni1:loni2,lati1:lati2,plotlevel)));
min_Q_mp_dn = min(min(Q_mp_dn(loni1:loni2,lati1:lati2,plotlevel)));
max_P_mp_dn = max(max(P_mp_dn(loni1:loni2,lati1:lati2,plotlevel)));
min_P_mp_dn = min(min(P_mp_dn(loni1:loni2,lati1:lati2,plotlevel)));

max_U_mp_mn = max(max(U_mp_mn(loni1:loni2,lati1:lati2,plotlevel)));
min_U_mp_mn = min(min(U_mp_mn(loni1:loni2,lati1:lati2,plotlevel)));
max_V_mp_mn = max(max(V_mp_mn(loni1:loni2,lati1:lati2,plotlevel)));
min_V_mp_mn = min(min(V_mp_mn(loni1:loni2,lati1:lati2,plotlevel)));
max_T_mp_mn = max(max(T_mp_mn(loni1:loni2,lati1:lati2,plotlevel)));
min_T_mp_mn = min(min(T_mp_mn(loni1:loni2,lati1:lati2,plotlevel)));
max_Q_mp_mn = max(max(Q_mp_mn(loni1:loni2,lati1:lati2,plotlevel)));
min_Q_mp_mn = min(min(Q_mp_mn(loni1:loni2,lati1:lati2,plotlevel)));
max_P_mp_mn = max(max(P_mp_mn(loni1:loni2,lati1:lati2,plotlevel)));
min_P_mp_mn = min(min(P_mp_mn(loni1:loni2,lati1:lati2,plotlevel)));

max_abs_u_dp_dn = max([max_U_dp_dn min_U_dp_dn]);
max_abs_v_dp_dn = max([max_V_dp_dn min_V_dp_dn]);
max_abs_t_dp_dn = max([max_T_dp_dn min_T_dp_dn]);
max_abs_q_dp_dn = max([max_Q_dp_dn min_Q_dp_dn]);
max_abs_p_dp_dn = max([max_P_dp_dn min_P_dp_dn]);

max_abs_u_dp_mn = max([max_U_dp_mn min_U_dp_mn]);
max_abs_v_dp_mn = max([max_V_dp_mn min_V_dp_mn]);
max_abs_t_dp_mn = max([max_T_dp_mn min_T_dp_mn]);
max_abs_q_dp_mn = max([max_Q_dp_mn min_Q_dp_mn]);
max_abs_p_dp_mn = max([max_P_dp_mn min_P_dp_mn]);

max_abs_u_mp_dn = max([max_U_mp_dn min_U_mp_dn]);
max_abs_v_mp_dn = max([max_V_mp_dn min_V_mp_dn]);
max_abs_t_mp_dn = max([max_T_mp_dn min_T_mp_dn]);
max_abs_q_mp_dn = max([max_Q_mp_dn min_Q_mp_dn]);
max_abs_p_mp_dn = max([max_P_mp_dn min_P_mp_dn]);

max_abs_u_mp_mn = max([max_U_mp_mn min_U_mp_mn]);
max_abs_v_mp_mn = max([max_V_mp_mn min_V_mp_mn]);
max_abs_t_mp_mn = max([max_T_mp_mn min_T_mp_mn]);
max_abs_q_mp_mn = max([max_Q_mp_mn min_Q_mp_mn]);
max_abs_p_mp_mn = max([max_P_mp_mn min_P_mp_mn]);

for i = loni1:loni2
    for j = lati1:lati2
        for k = plotlevel
        
            %dp_dn
            if abs(U_dp_dn(i,j,k)) < max_abs_u_dp_dn/10
                U_dp_dn(i,j,k) = 0.0;
            end
            if abs(V_dp_dn(i,j,k)) < max_abs_v_dp_dn/10
                V_dp_dn(i,j,k) = 0.0;
            end
            if abs(T_dp_dn(i,j,k)) < max_abs_t_dp_dn/10
                T_dp_dn(i,j,k) = 0.0;
            end
            if abs(Q_dp_dn(i,j,k)) < max_abs_q_dp_dn/10
                Q_dp_dn(i,j,k) = 0.0;
            end
            if abs(P_dp_dn(i,j,k)) < max_abs_p_dp_dn/10
                P_dp_dn(i,j,k) = 0.0;
            end
        
            %dp_mn
            if abs(U_dp_mn(i,j,k)) < max_abs_u_dp_mn/10
                U_dp_mn(i,j,k) = 0.0;
            end
            if abs(V_dp_mn(i,j,k)) < max_abs_v_dp_mn/10
                V_dp_mn(i,j,k) = 0.0;
            end
            if abs(T_dp_mn(i,j,k)) < max_abs_t_dp_mn/10
                T_dp_mn(i,j,k) = 0.0;
            end
            if abs(Q_dp_mn(i,j,k)) < max_abs_q_dp_mn/10
                Q_dp_mn(i,j,k) = 0.0;
            end
            if abs(P_dp_mn(i,j,k)) < max_abs_p_dp_mn/10
                P_dp_mn(i,j,k) = 0.0;
            end
        
            %mp_dn
            if abs(U_mp_dn(i,j,k)) < max_abs_u_mp_dn/10
                U_mp_dn(i,j,k) = 0.0;
            end
            if abs(V_mp_dn(i,j,k)) < max_abs_v_mp_dn/10
                V_mp_dn(i,j,k) = 0.0;
            end
            if abs(T_mp_dn(i,j,k)) < max_abs_t_mp_dn/10
                T_mp_dn(i,j,k) = 0.0;
            end
            if abs(Q_mp_dn(i,j,k)) < max_abs_q_mp_dn/10
                Q_mp_dn(i,j,k) = 0.0;
            end
            if abs(P_mp_dn(i,j,k)) < max_abs_p_mp_dn/10
                P_mp_dn(i,j,k) = 0.0;
            end
        
            %mp_mn        
            if abs(U_mp_mn(i,j,k)) < max_abs_u_mp_mn/10
                U_mp_mn(i,j,k) = 0.0;
            end
            if abs(V_mp_mn(i,j,k)) < max_abs_v_mp_mn/10
                V_mp_mn(i,j,k) = 0.0;
            end
            if abs(T_mp_mn(i,j,k)) < max_abs_t_mp_mn/10
                T_mp_mn(i,j,k) = 0.0;
            end
            if abs(Q_mp_mn(i,j,k)) < max_abs_q_mp_mn/10
                Q_mp_mn(i,j,k) = 0.0;
            end
            if abs(P_mp_mn(i,j,k)) < max_abs_p_mp_mn/10
                P_mp_mn(i,j,k) = 0.0;
            end
            
        end
    end
end

max_abs_u_dpdn_vs_mpdn = max([max_abs_u_mp_dn max_abs_u_dp_dn]);
max_abs_v_dpdn_vs_mpdn = max([max_abs_v_mp_dn max_abs_v_dp_dn]);
max_abs_t_dpdn_vs_mpdn = max([max_abs_t_mp_dn max_abs_t_dp_dn]);
max_abs_q_dpdn_vs_mpdn = max([max_abs_q_mp_dn max_abs_q_dp_dn]);
max_abs_p_dpdn_vs_mpdn = max([max_abs_p_mp_dn max_abs_p_dp_dn]);

max_abs_u_dpdn_vs_mpmn = max([max_abs_u_mp_mn max_abs_u_dp_dn]);
max_abs_v_dpdn_vs_mpmn = max([max_abs_v_mp_mn max_abs_v_dp_dn]);
max_abs_t_dpdn_vs_mpmn = max([max_abs_t_mp_mn max_abs_t_dp_dn]);
max_abs_q_dpdn_vs_mpmn = max([max_abs_q_mp_mn max_abs_q_dp_dn]);
max_abs_p_dpdn_vs_mpmn = max([max_abs_p_mp_mn max_abs_p_dp_dn]);

lon1 = -100;
lon2 = 0;

lat1 = 10;
lat2 = 80;

lati1 = 101;
lati2 = 171;
loni1 = 65;
loni2 = 145;

grey = 0.5;

line_wid = 1;
fontsize = 12;

boxhor = -42.5:-38.5;
boxver = 54:58;

% for plotlev = 40:72
%     
%     contourf(lon(loni1:loni2),lat(lati1:lati2),Q_mp_dn(loni1:loni2,lati1:lati2,plotlev)','LineStyle','none');
%     pause
%     close
%     
% end
% 
% asd

% % UNCOMMENT TO CHECK BOX IS ALLIGNED PROPERLY
% figure
% contour(lon,lat,QV(:,:,63)')
% hold on
% plot(boxhor,boxver(1)*ones(1,5),'k')
% plot(boxhor,boxver(end)*ones(1,5),'k')
% plot(boxhor(1)*ones(1,5),boxver,'k')
% plot(boxhor(end)*ones(1,5),boxver,'k')
% plot(coast_lon,coast_lat,'Color',[0.8 0.8 0.8])
% 
% colorbar
% xlim([lon1 lon2])
% ylim([lat1 lat2])
% 
% % UNCOMMENT TO CHECK PRESSURE AND CONVECTIVE RAIN AT END OF WINDOW
% figure
% [C,h] = contour(lon(loni1:loni2),lat(lati1:lati2),PRECCON_2(loni1:loni2,lati1:lati2)');
% caxis([0 max(max(PRECCON_2(loni1:loni2,lati1:lati2)))])
% colormap(winter)
% pos1 = get(gca,'Position');
% colbar = colorbar;
% hold on
% [C1,h1] = contour(lon(loni1:loni2),lat(lati1:lati2),SLP_2(loni1:loni2,lati1:lati2)','Color','k');
% set(h1,'ShowText','on')
% plot(boxhor,boxver(1)*ones(1,5),'k')
% plot(boxhor,boxver(end)*ones(1,5),'k')
% plot(boxhor(1)*ones(1,5),boxver,'k')
% plot(boxhor(end)*ones(1,5),boxver,'k')
% plot(coast_lon,coast_lat,'Color',[grey grey grey])
% asd


cmap1 = colormap(winter);
% cmap1(1,:) = [1 1 1];

cmap2 = colormap(autumn);
cmap2 = flipud(cmap2);
cmap2(32,:) = [1 1 1];
cmap2(33,:) = [1 1 1];

close

colbar1_ticks = zeros(4,10);


figure('position',[ 88 137 1191 782])

subplot(2,2,1)
pos_sp_1 = get(gca,'position');
subplot(2,2,2)
pos_sp_2 = get(gca,'position');
subplot(2,2,3)
pos_sp_3 = get(gca,'position');
subplot(2,2,4)
pos_sp_4 = get(gca,'position');

close

figure('position',[ 88 137 1191 782])


%%%%%% U PLOTS %%%%%%%
ax(1,1) = axes('position',[pos_sp_1(1) pos_sp_1(2) pos_sp_1(3) pos_sp_1(4) ]);
[C1,h1] = contourf(lon(loni1:loni2),lat(lati1:lati2),U_dp_dn(loni1:loni2,lati1:lati2,plotlevel)','LineStyle','none');
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
[C2,h2] = contour(lon(loni1:loni2),lat(lati1:lati2),SLP_1(loni1:loni2,lati1:lati2)','Color','k','ShowText','on');
colbar(1,1) = colorbar('location','WestOutside');

%Reduce contour Step
set(h1,'LevelStep',0.5*get(h1,'LevelStep'));

caxis([-max_abs_u_dpdn_vs_mpdn max_abs_u_dpdn_vs_mpdn])

poscbar1 = get(colbar(1,1),'Position');
set(colbar(1,1),'Position',[poscbar1(1)-0.07 poscbar1(2) poscbar1(3) poscbar1(4)])

colbar1_ticks(1,1:length(get(colbar(1,1),'Ytick'))) = get(colbar(1,1),'Ytick');
clim_old_1(:,1) = get(ax(1,1),'CLim');

%%%%%%%
ax(1,2) = axes('position',[pos_sp_1(1) pos_sp_1(2) pos_sp_1(3) pos_sp_1(4) ]);
[C3,h3] = contour(lon(loni1:loni2),lat(lati1:lati2),PRECCON_1(loni1:loni2,lati1:lati2)');
hold on
plot(boxhor,boxver(1)*ones(1,5),'k')
plot(boxhor,boxver(end)*ones(1,5),'k')
plot(boxhor(1)*ones(1,5),boxver,'k')
plot(boxhor(end)*ones(1,5),boxver,'k')
set(h3,'LevelStep',0.75*get(h3,'LevelStep'))
colormap(cmap1)

%DETERMINE AND RESET RANGES OF COLORMAP AND COLORBARS
colormap([cmap2;cmap1]);
cm_length = size(get(gcf,'colormap'),1);
range(1,1:2) = [1 size(cmap1,1)];
range(2,1:2) = [1 size(cmap2,1)]+range(1,2);
cmap_ind_min = range(:,1);
cmap_ind_max = range(:,2);
cdata_range(1,:) = get(ax(1,1),'CLim');
cdata_range(2,:) = get(ax(1,2),'CLim');
color_data_min = cdata_range(:,1);
color_data_max = cdata_range(:,2);

cmin=((cmap_ind_max-1).*color_data_min-(cmap_ind_min-1).*color_data_max)./(cmap_ind_max-cmap_ind_min);
cmax=((cmap_ind_max-1).*color_data_min-(cmap_ind_min-1).*color_data_max+cm_length*(color_data_max-color_data_min))./(cmap_ind_max-cmap_ind_min);

set(ax(1,1),'CLim',[cmin(1),cmax(1)]);
set(ax(1,2),'CLim',[cmin(2),cmax(2)]);

set(get(colbar(1,1),'children'),'CData',(cmap_ind_min(1):cmap_ind_max(1))');
% set(get(colbar(1,2),'children'),'CData',(cmap_ind_min(2):cmap_ind_max(2))');

title('(a) \partial J/\partial u\prime (Jkg^{-1}m^{-1}s)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca, 'color', 'none');
box on; 
set(ax(1,1),'XTick',[],'Ytick',[]);
set(ax(1,2),'YTick',[10 24 38 52 66 80],'Xtick',[-100 -80 -60 -40 -20 0]);
set(ax(1,2),'XTickLabel',{'100W', '80W', '60W', '40W', '20W', '0'},'YtickLabel',{'10N' '24N' '38N' '52N' '66N' '80N'});

set(ax(1,1),'FontSize',fontsize,'FontName','TimesNewRoman')
set(ax(1,2),'FontSize',fontsize,'FontName','TimesNewRoman')


%%%%%% V PLOTS %%%%%%%
ax(2,1) = axes('position',[pos_sp_2(1) pos_sp_2(2) pos_sp_2(3) pos_sp_2(4) ]);
[C1,h1] = contourf(lon(loni1:loni2),lat(lati1:lati2),V_dp_dn(loni1:loni2,lati1:lati2,plotlevel)','LineStyle','none');
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
[C2,h2] = contour(lon(loni1:loni2),lat(lati1:lati2),SLP_1(loni1:loni2,lati1:lati2)','Color','k','ShowText','on');
colbar(2,1) = colorbar('location','WestOutside');

set(h1,'LevelStep',0.5*get(h1,'LevelStep'));

caxis([-max_abs_v_dpdn_vs_mpdn max_abs_v_dpdn_vs_mpdn])

poscbar1 = get(colbar(2,1),'Position');
set(colbar(2,1),'Position',[poscbar1(1)-0.07 poscbar1(2) poscbar1(3) poscbar1(4)])

colbar1_ticks(2,1:length(get(colbar(2,1),'Ytick'))) = get(colbar(2,1),'Ytick');
clim_old_1(:,2) = get(ax(2,1),'CLim');

%%%%%%%
ax(2,2) = axes('position',[pos_sp_2(1) pos_sp_2(2) pos_sp_2(3) pos_sp_2(4) ]);
[C3,h3] = contour(lon(loni1:loni2),lat(lati1:lati2),PRECCON_1(loni1:loni2,lati1:lati2)');
hold on
plot(boxhor,boxver(1)*ones(1,5),'k')
plot(boxhor,boxver(end)*ones(1,5),'k')
plot(boxhor(1)*ones(1,5),boxver,'k')
plot(boxhor(end)*ones(1,5),boxver,'k')
set(h3,'LevelStep',0.75*get(h3,'LevelStep'))
colormap(cmap1)

%DETERMINE AND RESET RANGES OF COLORMAP AND COLORBARS
colormap([cmap2;cmap1]);
cm_length = size(get(gcf,'colormap'),1);
range(1,1:2) = [1 size(cmap1,1)];
range(2,1:2) = [1 size(cmap2,1)]+range(1,2);
cmap_ind_min = range(:,1);
cmap_ind_max = range(:,2);
cdata_range(1,:) = get(ax(2,1),'CLim');
cdata_range(2,:) = get(ax(2,2),'CLim');
color_data_min = cdata_range(:,1);
color_data_max = cdata_range(:,2);

cmin=((cmap_ind_max-1).*color_data_min-(cmap_ind_min-1).*color_data_max)./(cmap_ind_max-cmap_ind_min);
cmax=((cmap_ind_max-1).*color_data_min-(cmap_ind_min-1).*color_data_max+cm_length*(color_data_max-color_data_min))./(cmap_ind_max-cmap_ind_min);

set(ax(2,1),'CLim',[cmin(1),cmax(1)]);
set(ax(2,2),'CLim',[cmin(2),cmax(2)]);

set(get(colbar(2,1),'children'),'CData',(cmap_ind_min(1):cmap_ind_max(1))');

title('(b) \partial J/\partial v\prime (Jkg^{-1}m^{-1}s)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca, 'color', 'none');
box on; 
set(ax(2,1),'XTick',[],'Ytick',[]);
set(ax(2,2),'YTick',[10 24 38 52 66 80],'Xtick',[-100 -80 -60 -40 -20 0]);
set(ax(2,2),'XTickLabel',{'100W', '80W', '60W', '40W', '20W', '0'},'YtickLabel',{'10N' '24N' '38N' '52N' '66N' '80N'});

set(ax(2,1),'FontSize',fontsize,'FontName','TimesNewRoman')
set(ax(2,2),'FontSize',fontsize,'FontName','TimesNewRoman')


%%%%% T PLOTS %%%%%%%
ax(3,1) = axes('position',[pos_sp_3(1) pos_sp_3(2) pos_sp_3(3) pos_sp_3(4) ]);
[C1,h1] = contourf(lon(loni1:loni2),lat(lati1:lati2),T_dp_dn(loni1:loni2,lati1:lati2,plotlevel)','LineStyle','none');
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
[C2,h2] = contour(lon(loni1:loni2),lat(lati1:lati2),SLP_1(loni1:loni2,lati1:lati2)','Color','k','ShowText','on');
colbar(3,1) = colorbar('location','WestOutside');

set(h1,'LevelStep',0.5*get(h1,'LevelStep'));

caxis([-max_abs_t_dpdn_vs_mpdn max_abs_t_dpdn_vs_mpdn])

poscbar1 = get(colbar(3,1),'Position');
set(colbar(3,1),'Position',[poscbar1(1)-0.07 poscbar1(2) poscbar1(3) poscbar1(4)])

colbar1_ticks(3,1:length(get(colbar(3,1),'Ytick'))) = get(colbar(3,1),'Ytick');
clim_old_1(:,3) = get(ax(3,1),'CLim');

colbar1 = get(colbar(3,1),'Ytick');

%%%%%%%
ax(3,2) = axes('position',[pos_sp_3(1) pos_sp_3(2) pos_sp_3(3) pos_sp_3(4) ]);
[C3,h3] = contour(lon(loni1:loni2),lat(lati1:lati2),PRECCON_1(loni1:loni2,lati1:lati2)');
hold on
plot(boxhor,boxver(1)*ones(1,5),'k')
plot(boxhor,boxver(end)*ones(1,5),'k')
plot(boxhor(1)*ones(1,5),boxver,'k')
plot(boxhor(end)*ones(1,5),boxver,'k')
set(h3,'LevelStep',0.75*get(h3,'LevelStep'))
colormap(cmap1)

%DETERMINE AND RESET RANGES OF COLORMAP AND COLORBARS
colormap([cmap2;cmap1]);
cm_length = size(get(gcf,'colormap'),1);
range(1,1:2) = [1 size(cmap1,1)];
range(2,1:2) = [1 size(cmap2,1)]+range(1,2);
cmap_ind_min = range(:,1);
cmap_ind_max = range(:,2);
cdata_range(1,:) = get(ax(3,1),'CLim');
cdata_range(2,:) = get(ax(3,2),'CLim');
color_data_min = cdata_range(:,1);
color_data_max = cdata_range(:,2);

cmin=((cmap_ind_max-1).*color_data_min-(cmap_ind_min-1).*color_data_max)./(cmap_ind_max-cmap_ind_min);
cmax=((cmap_ind_max-1).*color_data_min-(cmap_ind_min-1).*color_data_max+cm_length*(color_data_max-color_data_min))./(cmap_ind_max-cmap_ind_min);

set(ax(3,1),'CLim',[cmin(1),cmax(1)]);
set(ax(3,2),'CLim',[cmin(2),cmax(2)]);

set(get(colbar(3,1),'children'),'CData',(cmap_ind_min(1):cmap_ind_max(1))');

title('(c) \partial J/\partial T_v\prime (Jkg^{-1}K^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca, 'color', 'none');
box on; 
set(ax(3,1),'XTick',[],'Ytick',[]);
set(ax(3,2),'YTick',[10 24 38 52 66 80],'Xtick',[-100 -80 -60 -40 -20 0]);
set(ax(3,2),'XTickLabel',{'100W', '80W', '60W', '40W', '20W', '0'},'YtickLabel',{'10N' '24N' '38N' '52N' '66N' '80N'});

set(ax(3,1),'FontSize',fontsize,'FontName','TimesNewRoman')
set(ax(3,2),'FontSize',fontsize,'FontName','TimesNewRoman')

%%%%%% Q PLOTS %%%%%%%
ax(4,1) = axes('position',[pos_sp_4(1) pos_sp_4(2) pos_sp_4(3) pos_sp_4(4) ]);
[C1,h1] = contourf(lon(loni1:loni2),lat(lati1:lati2),Q_dp_dn(loni1:loni2,lati1:lati2,plotlevel)','LineStyle','none');
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
[C2,h2] = contour(lon(loni1:loni2),lat(lati1:lati2),SLP_1(loni1:loni2,lati1:lati2)','Color','k','ShowText','on');
colbar(4,1) = colorbar('location','WestOutside');

set(h1,'LevelStep',0.5*get(h1,'LevelStep'));

caxis([-max_abs_q_dp_dn max_abs_q_dp_dn])

poscbar1 = get(colbar(4,1),'Position');
set(colbar(4,1),'Position',[poscbar1(1)-0.07 poscbar1(2) poscbar1(3) poscbar1(4)])

colbar1_ticks(4,1:length(get(colbar(4,1),'Ytick'))) = get(colbar(4,1),'Ytick');
clim_old_1(:,4) = get(ax(4,1),'CLim');


%%%%%%%
ax(4,2) = axes('position',[pos_sp_4(1) pos_sp_4(2) pos_sp_4(3) pos_sp_4(4) ]);
[C3,h3] = contour(lon(loni1:loni2),lat(lati1:lati2),PRECCON_1(loni1:loni2,lati1:lati2)');
hold on
plot(boxhor,boxver(1)*ones(1,5),'k')
plot(boxhor,boxver(end)*ones(1,5),'k')
plot(boxhor(1)*ones(1,5),boxver,'k')
plot(boxhor(end)*ones(1,5),boxver,'k')
set(h3,'LevelStep',0.75*get(h3,'LevelStep'))
colormap(cmap1)
colbar(4,2) = colorbar;

colbar2_ticks = get(colbar(4,2),'Ytick');
clim_old_2 = get(ax(4,2),'CLim');

%DETERMINE AND RESET RANGES OF COLORMAP AND COLORBARS
colormap([cmap2;cmap1]);
cm_length = size(get(gcf,'colormap'),1);
range(1,1:2) = [1 size(cmap1,1)];
range(2,1:2) = [1 size(cmap2,1)]+range(1,2);
cmap_ind_min = range(:,1);
cmap_ind_max = range(:,2);
cdata_range(1,:) = get(ax(4,1),'CLim');
cdata_range(2,:) = get(ax(4,2),'CLim');
color_data_min = cdata_range(:,1);
color_data_max = cdata_range(:,2);

cmin=((cmap_ind_max-1).*color_data_min-(cmap_ind_min-1).*color_data_max)./(cmap_ind_max-cmap_ind_min);
cmax=((cmap_ind_max-1).*color_data_min-(cmap_ind_min-1).*color_data_max+cm_length*(color_data_max-color_data_min))./(cmap_ind_max-cmap_ind_min);

set(ax(4,1),'CLim',[cmin(1),cmax(1)]);
set(ax(4,2),'CLim',[cmin(2),cmax(2)]);

set(get(colbar(1,1),'children'),'CData',(cmap_ind_min(1):cmap_ind_max(1))');
% set(get(colbar(1,2),'children'),'CData',(cmap_ind_min(2):cmap_ind_max(2))');

set(get(colbar(2,1),'children'),'CData',(cmap_ind_min(1):cmap_ind_max(1))');
% set(get(colbar(2,2),'children'),'CData',(cmap_ind_min(2):cmap_ind_max(2))');

set(get(colbar(3,1),'children'),'CData',(cmap_ind_min(1):cmap_ind_max(1))');
% set(get(colbar(3,2),'children'),'CData',(cmap_ind_min(2):cmap_ind_max(2))');

set(get(colbar(4,1),'children'),'CData',(cmap_ind_min(1):cmap_ind_max(1))');
set(get(colbar(4,2),'children'),'CData',(cmap_ind_min(2):cmap_ind_max(2))');

%RESET POSITIONS WITH NEW COLORBARS
pos1 = get(ax(4,1),'position');
set(ax(4,2),'position',pos1)

title('(d) \partial J/\partial q\prime (Jkg^{-1}kg^{-1}kg)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca, 'color', 'none');
box on; 
set(ax(4,1),'XTick',[],'Ytick',[]);
set(ax(4,2),'YTick',[10 24 38 52 66 80],'Xtick',[-100 -80 -60 -40 -20 0]);
set(ax(4,2),'XTickLabel',{'100W', '80W', '60W', '40W', '20W', '0'},'YtickLabel',{'10N' '24N' '38N' '52N' '66N' '80N'});

set(ax(4,1),'FontSize',fontsize,'FontName','TimesNewRoman')
set(ax(4,2),'FontSize',fontsize,'FontName','TimesNewRoman')

clim_new_1 = get(ax(1,1),'CLim');
clim_new_1_diff = clim_new_1(2)-clim_new_1(1);
clim_new_1_inc = (colbar1_ticks(1,1:length(get(colbar(1,1),'Ytick')))-clim_old_1(1,1))./(clim_old_1(2,1)-clim_old_1(1,1));
clim_new_1 = clim_new_1(1) + clim_new_1_inc*clim_new_1_diff;
set(colbar(1,1),'YTick',clim_new_1)

set(colbar(1,1),'YtickLabel',{'-4', '-3', '-2', '-1', '0', '1', '2', '3', '4'})
title(colbar(1,1),'\times 10^{-6}','FontSize',fontsize,'FontName','TimesNewRoman')

clim_new_1 = get(ax(2,1),'CLim');
clim_new_1_diff = clim_new_1(2)-clim_new_1(1);
clim_new_1_inc = (colbar1_ticks(2,1:length(get(colbar(2,1),'Ytick')))-clim_old_1(1,2))./(clim_old_1(2,2)-clim_old_1(1,2));
clim_new_1 = clim_new_1(1) + clim_new_1_inc*clim_new_1_diff;
set(colbar(2,1),'YTick',clim_new_1(1:end-1))

set(colbar(2,1),'YtickLabel',{'-2', '-1.5', '-1', '-0.5', '0', '0.5', '1.0', '1.5', '2.0'})
title(colbar(2,1),'\times 10^{-6}','FontSize',fontsize,'FontName','TimesNewRoman')

clim_new_1 = get(ax(3,1),'CLim');
clim_new_1_diff = clim_new_1(2)-clim_new_1(1);
clim_new_1_inc = (colbar1_ticks(3,1:length(colbar1))-clim_old_1(1,3))./(clim_old_1(2,3)-clim_old_1(1,3));
clim_new_1 = clim_new_1(1) + clim_new_1_inc*clim_new_1_diff;
set(colbar(3,1),'YTick',clim_new_1)

set(colbar(3,1),'YtickLabel',{'-1.5', '-1', '-0.5', '0', '0.5', '1.0', '1.5'})
title(colbar(3,1),'\times 10^{-5}','FontSize',fontsize,'FontName','TimesNewRoman')

clim_new_1 = get(ax(4,1),'CLim');
clim_new_1_diff = clim_new_1(2)-clim_new_1(1);
clim_new_1_inc = (colbar1_ticks(4,1:length(get(colbar(4,1),'Ytick')))-clim_old_1(1,4))./(clim_old_1(2,4)-clim_old_1(1,4));
clim_new_1 = clim_new_1(1) + clim_new_1_inc*clim_new_1_diff;
set(colbar(4,1),'YTick',clim_new_1)

set(colbar(4,1),'YtickLabel',{'-2', '-1.5', '-1', '-0.5', '0', '0.5', '1.0', '1.5', '2.0'})
title(colbar(4,1),'\times 10^{-5}','FontSize',fontsize,'FontName','TimesNewRoman')

clim_new_2 = get(ax(4,2),'CLim');
clim_new_2_diff = clim_new_2(2)-clim_new_2(1);
clim_new_2_inc = (colbar2_ticks-clim_old_2(1))./(clim_old_2(2)-clim_old_2(1));
clim_new_2 = clim_new_2(1) + clim_new_2_inc*clim_new_2_diff;
% set(colbar(1,2),'YTick',clim_new_2)
% set(colbar(2,2),'YTick',clim_new_2)
% set(colbar(3,2),'YTick',clim_new_2)
set(colbar(4,2),'YTick',clim_new_2)

poscbar2 = get(colbar(4,2),'Position');
set(colbar(4,2),'Position',[0.93 0.28 poscbar2(3) 0.45])
ylabel(colbar(4,2),'Convective Precipitation Rate (kgm^{-2}s^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')

tickstr = num2str(colbar2_ticks,1);
% set(colbar(1,2),'YtickLabel',{'2', '4', '6', '8', '10', '12'})
% title(colbar(1,2),'\times 10^{-4}','FontSize',fontsize,'FontName','TimesNewRoman')
% 
% set(colbar(2,2),'YtickLabel',{'2', '4', '6', '8', '10', '12'})
% title(colbar(2,2),'\times 10^{-4}','FontSize',fontsize,'FontName','TimesNewRoman')
% 
% 
% set(colbar(3,2),'YtickLabel',{'2', '4', '6', '8', '10', '12'})
% title(colbar(3,2),'\times 10^{-4}','FontSize',fontsize,'FontName','TimesNewRoman')


set(colbar(4,2),'YtickLabel',{'2', '4', '6', '8', '10', '12'})
title(colbar(4,2),'\times 10^{-4}','FontSize',fontsize,'FontName','TimesNewRoman')








saveas(gcf,'storm_sensitivity_dp_dn.eps', 'psc2')
















figure('position',[ 88 137 1191 782])

subplot(2,2,1)
pos_sp_1 = get(gca,'position');
subplot(2,2,2)
pos_sp_2 = get(gca,'position');
subplot(2,2,3)
pos_sp_3 = get(gca,'position');
subplot(2,2,4)
pos_sp_4 = get(gca,'position');

close

figure('position',[ 88 137 1191 782])


%%%%%% U PLOTS %%%%%%%
ax(1,1) = axes('position',[pos_sp_1(1) pos_sp_1(2) pos_sp_1(3) pos_sp_1(4) ]);
[C1,h1] = contourf(lon(loni1:loni2),lat(lati1:lati2),U_mp_dn(loni1:loni2,lati1:lati2,plotlevel)','LineStyle','none');
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
[C2,h2] = contour(lon(loni1:loni2),lat(lati1:lati2),SLP_1(loni1:loni2,lati1:lati2)','Color','k','ShowText','on');
colbar(1,1) = colorbar('location','WestOutside');

%Reduce contour Step
set(h1,'LevelStep',0.5*get(h1,'LevelStep'));

caxis([-max_abs_u_dpdn_vs_mpdn max_abs_u_dpdn_vs_mpdn])

poscbar1 = get(colbar(1,1),'Position');
set(colbar(1,1),'Position',[poscbar1(1)-0.07 poscbar1(2) poscbar1(3) poscbar1(4)])

colbar1_ticks(1,1:length(get(colbar(1,1),'Ytick'))) = get(colbar(1,1),'Ytick');
clim_old_1(:,1) = get(ax(1,1),'CLim');

%%%%%%%
ax(1,2) = axes('position',[pos_sp_1(1) pos_sp_1(2) pos_sp_1(3) pos_sp_1(4) ]);
[C3,h3] = contour(lon(loni1:loni2),lat(lati1:lati2),PRECCON_1(loni1:loni2,lati1:lati2)');
hold on
plot(boxhor,boxver(1)*ones(1,5),'k')
plot(boxhor,boxver(end)*ones(1,5),'k')
plot(boxhor(1)*ones(1,5),boxver,'k')
plot(boxhor(end)*ones(1,5),boxver,'k')
set(h3,'LevelStep',0.75*get(h3,'LevelStep'))
colormap(cmap1)

%DETERMINE AND RESET RANGES OF COLORMAP AND COLORBARS
colormap([cmap2;cmap1]);
cm_length = size(get(gcf,'colormap'),1);
range(1,1:2) = [1 size(cmap1,1)];
range(2,1:2) = [1 size(cmap2,1)]+range(1,2);
cmap_ind_min = range(:,1);
cmap_ind_max = range(:,2);
cdata_range(1,:) = get(ax(1,1),'CLim');
cdata_range(2,:) = get(ax(1,2),'CLim');
color_data_min = cdata_range(:,1);
color_data_max = cdata_range(:,2);

cmin=((cmap_ind_max-1).*color_data_min-(cmap_ind_min-1).*color_data_max)./(cmap_ind_max-cmap_ind_min);
cmax=((cmap_ind_max-1).*color_data_min-(cmap_ind_min-1).*color_data_max+cm_length*(color_data_max-color_data_min))./(cmap_ind_max-cmap_ind_min);

set(ax(1,1),'CLim',[cmin(1),cmax(1)]);
set(ax(1,2),'CLim',[cmin(2),cmax(2)]);

set(get(colbar(1,1),'children'),'CData',(cmap_ind_min(1):cmap_ind_max(1))');
% set(get(colbar(1,2),'children'),'CData',(cmap_ind_min(2):cmap_ind_max(2))');

title('(a) \partial J/\partial u\prime (Jkg^{-1}m^{-1}s)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca, 'color', 'none');
box on; 
set(ax(1,1),'XTick',[],'Ytick',[]);
set(ax(1,2),'YTick',[10 24 38 52 66 80],'Xtick',[-100 -80 -60 -40 -20 0]);
set(ax(1,2),'XTickLabel',{'100W', '80W', '60W', '40W', '20W', '0'},'YtickLabel',{'10N' '24N' '38N' '52N' '66N' '80N'});

set(ax(1,1),'FontSize',fontsize,'FontName','TimesNewRoman')
set(ax(1,2),'FontSize',fontsize,'FontName','TimesNewRoman')


%%%%%% V PLOTS %%%%%%%
ax(2,1) = axes('position',[pos_sp_2(1) pos_sp_2(2) pos_sp_2(3) pos_sp_2(4) ]);
[C1,h1] = contourf(lon(loni1:loni2),lat(lati1:lati2),V_mp_dn(loni1:loni2,lati1:lati2,plotlevel)','LineStyle','none');
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
[C2,h2] = contour(lon(loni1:loni2),lat(lati1:lati2),SLP_1(loni1:loni2,lati1:lati2)','Color','k','ShowText','on');
colbar(2,1) = colorbar('location','WestOutside');

set(h1,'LevelStep',0.5*get(h1,'LevelStep'));

caxis([-max_abs_v_dpdn_vs_mpdn max_abs_v_dpdn_vs_mpdn])

poscbar1 = get(colbar(2,1),'Position');
set(colbar(2,1),'Position',[poscbar1(1)-0.07 poscbar1(2) poscbar1(3) poscbar1(4)])

colbar1_ticks(2,1:length(get(colbar(2,1),'Ytick'))) = get(colbar(2,1),'Ytick');
clim_old_1(:,2) = get(ax(2,1),'CLim');

%%%%%%%
ax(2,2) = axes('position',[pos_sp_2(1) pos_sp_2(2) pos_sp_2(3) pos_sp_2(4) ]);
[C3,h3] = contour(lon(loni1:loni2),lat(lati1:lati2),PRECCON_1(loni1:loni2,lati1:lati2)');
hold on
plot(boxhor,boxver(1)*ones(1,5),'k')
plot(boxhor,boxver(end)*ones(1,5),'k')
plot(boxhor(1)*ones(1,5),boxver,'k')
plot(boxhor(end)*ones(1,5),boxver,'k')
set(h3,'LevelStep',0.75*get(h3,'LevelStep'))
colormap(cmap1)

%DETERMINE AND RESET RANGES OF COLORMAP AND COLORBARS
colormap([cmap2;cmap1]);
cm_length = size(get(gcf,'colormap'),1);
range(1,1:2) = [1 size(cmap1,1)];
range(2,1:2) = [1 size(cmap2,1)]+range(1,2);
cmap_ind_min = range(:,1);
cmap_ind_max = range(:,2);
cdata_range(1,:) = get(ax(2,1),'CLim');
cdata_range(2,:) = get(ax(2,2),'CLim');
color_data_min = cdata_range(:,1);
color_data_max = cdata_range(:,2);

cmin=((cmap_ind_max-1).*color_data_min-(cmap_ind_min-1).*color_data_max)./(cmap_ind_max-cmap_ind_min);
cmax=((cmap_ind_max-1).*color_data_min-(cmap_ind_min-1).*color_data_max+cm_length*(color_data_max-color_data_min))./(cmap_ind_max-cmap_ind_min);

set(ax(2,1),'CLim',[cmin(1),cmax(1)]);
set(ax(2,2),'CLim',[cmin(2),cmax(2)]);

set(get(colbar(2,1),'children'),'CData',(cmap_ind_min(1):cmap_ind_max(1))');

title('(b) \partial J/\partial v\prime (Jkg^{-1}m^{-1}s)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca, 'color', 'none');
box on; 
set(ax(2,1),'XTick',[],'Ytick',[]);
set(ax(2,2),'YTick',[10 24 38 52 66 80],'Xtick',[-100 -80 -60 -40 -20 0]);
set(ax(2,2),'XTickLabel',{'100W', '80W', '60W', '40W', '20W', '0'},'YtickLabel',{'10N' '24N' '38N' '52N' '66N' '80N'});

set(ax(2,1),'FontSize',fontsize,'FontName','TimesNewRoman')
set(ax(2,2),'FontSize',fontsize,'FontName','TimesNewRoman')


%%%%% T PLOTS %%%%%%%
ax(3,1) = axes('position',[pos_sp_3(1) pos_sp_3(2) pos_sp_3(3) pos_sp_3(4) ]);
[C1,h1] = contourf(lon(loni1:loni2),lat(lati1:lati2),T_mp_dn(loni1:loni2,lati1:lati2,plotlevel)','LineStyle','none');
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
[C2,h2] = contour(lon(loni1:loni2),lat(lati1:lati2),SLP_1(loni1:loni2,lati1:lati2)','Color','k','ShowText','on');
colbar(3,1) = colorbar('location','WestOutside');

set(h1,'LevelStep',0.5*get(h1,'LevelStep'));

caxis([-max_abs_t_dpdn_vs_mpdn max_abs_t_dpdn_vs_mpdn])

poscbar1 = get(colbar(3,1),'Position');
set(colbar(3,1),'Position',[poscbar1(1)-0.07 poscbar1(2) poscbar1(3) poscbar1(4)])

colbar1_ticks(3,1:length(get(colbar(3,1),'Ytick'))) = get(colbar(3,1),'Ytick');
clim_old_1(:,3) = get(ax(3,1),'CLim');

colbar1 = get(colbar(3,1),'Ytick');

%%%%%%%
ax(3,2) = axes('position',[pos_sp_3(1) pos_sp_3(2) pos_sp_3(3) pos_sp_3(4) ]);
[C3,h3] = contour(lon(loni1:loni2),lat(lati1:lati2),PRECCON_1(loni1:loni2,lati1:lati2)');
hold on
plot(boxhor,boxver(1)*ones(1,5),'k')
plot(boxhor,boxver(end)*ones(1,5),'k')
plot(boxhor(1)*ones(1,5),boxver,'k')
plot(boxhor(end)*ones(1,5),boxver,'k')
set(h3,'LevelStep',0.75*get(h3,'LevelStep'))
colormap(cmap1)

%DETERMINE AND RESET RANGES OF COLORMAP AND COLORBARS
colormap([cmap2;cmap1]);
cm_length = size(get(gcf,'colormap'),1);
range(1,1:2) = [1 size(cmap1,1)];
range(2,1:2) = [1 size(cmap2,1)]+range(1,2);
cmap_ind_min = range(:,1);
cmap_ind_max = range(:,2);
cdata_range(1,:) = get(ax(3,1),'CLim');
cdata_range(2,:) = get(ax(3,2),'CLim');
color_data_min = cdata_range(:,1);
color_data_max = cdata_range(:,2);

cmin=((cmap_ind_max-1).*color_data_min-(cmap_ind_min-1).*color_data_max)./(cmap_ind_max-cmap_ind_min);
cmax=((cmap_ind_max-1).*color_data_min-(cmap_ind_min-1).*color_data_max+cm_length*(color_data_max-color_data_min))./(cmap_ind_max-cmap_ind_min);

set(ax(3,1),'CLim',[cmin(1),cmax(1)]);
set(ax(3,2),'CLim',[cmin(2),cmax(2)]);

set(get(colbar(3,1),'children'),'CData',(cmap_ind_min(1):cmap_ind_max(1))');

title('(c) \partial J/\partial T_v\prime (Jkg^{-1}K^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca, 'color', 'none');
box on; 
set(ax(3,1),'XTick',[],'Ytick',[]);
set(ax(3,2),'YTick',[10 24 38 52 66 80],'Xtick',[-100 -80 -60 -40 -20 0]);
set(ax(3,2),'XTickLabel',{'100W', '80W', '60W', '40W', '20W', '0'},'YtickLabel',{'10N' '24N' '38N' '52N' '66N' '80N'});

set(ax(3,1),'FontSize',fontsize,'FontName','TimesNewRoman')
set(ax(3,2),'FontSize',fontsize,'FontName','TimesNewRoman')



%%%%%% Q PLOTS %%%%%%%
ax(4,1) = axes('position',[pos_sp_4(1) pos_sp_4(2) pos_sp_4(3) pos_sp_4(4) ]);
[C1,h1] = contourf(lon(loni1:loni2),lat(lati1:lati2),Q_mp_dn(loni1:loni2,lati1:lati2,plotlevel)','LineStyle','none');
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
[C2,h2] = contour(lon(loni1:loni2),lat(lati1:lati2),SLP_1(loni1:loni2,lati1:lati2)','Color','k','ShowText','on');
colbar(4,1) = colorbar('location','WestOutside');

set(h1,'LevelStep',0.5*get(h1,'LevelStep'));

caxis([-max_abs_q_mp_dn max_abs_q_mp_dn])

poscbar1 = get(colbar(4,1),'Position');
set(colbar(4,1),'Position',[poscbar1(1)-0.07 poscbar1(2) poscbar1(3) poscbar1(4)])

colbar1_ticks(4,1:length(get(colbar(4,1),'Ytick'))) = get(colbar(4,1),'Ytick');
clim_old_1(:,4) = get(ax(4,1),'CLim');

%%%%%%%
ax(4,2) = axes('position',[pos_sp_4(1) pos_sp_4(2) pos_sp_4(3) pos_sp_4(4) ]);
[C3,h3] = contour(lon(loni1:loni2),lat(lati1:lati2),PRECCON_1(loni1:loni2,lati1:lati2)');
hold on
plot(boxhor,boxver(1)*ones(1,5),'k')
plot(boxhor,boxver(end)*ones(1,5),'k')
plot(boxhor(1)*ones(1,5),boxver,'k')
plot(boxhor(end)*ones(1,5),boxver,'k')
set(h3,'LevelStep',0.75*get(h3,'LevelStep'))
colormap(cmap1)
colbar(4,2) = colorbar;

colbar2_ticks = get(colbar(4,2),'Ytick');
clim_old_2 = get(ax(4,2),'CLim');

%DETERMINE AND RESET RANGES OF COLORMAP AND COLORBARS
colormap([cmap2;cmap1]);
cm_length = size(get(gcf,'colormap'),1);
range(1,1:2) = [1 size(cmap1,1)];
range(2,1:2) = [1 size(cmap2,1)]+range(1,2);
cmap_ind_min = range(:,1);
cmap_ind_max = range(:,2);
cdata_range(1,:) = get(ax(4,1),'CLim');
cdata_range(2,:) = get(ax(4,2),'CLim');
color_data_min = cdata_range(:,1);
color_data_max = cdata_range(:,2);

cmin=((cmap_ind_max-1).*color_data_min-(cmap_ind_min-1).*color_data_max)./(cmap_ind_max-cmap_ind_min);
cmax=((cmap_ind_max-1).*color_data_min-(cmap_ind_min-1).*color_data_max+cm_length*(color_data_max-color_data_min))./(cmap_ind_max-cmap_ind_min);

set(ax(4,1),'CLim',[cmin(1),cmax(1)]);
set(ax(4,2),'CLim',[cmin(2),cmax(2)]);

set(get(colbar(1,1),'children'),'CData',(cmap_ind_min(1):cmap_ind_max(1))');
% set(get(colbar(1,2),'children'),'CData',(cmap_ind_min(2):cmap_ind_max(2))');

set(get(colbar(2,1),'children'),'CData',(cmap_ind_min(1):cmap_ind_max(1))');
% set(get(colbar(2,2),'children'),'CData',(cmap_ind_min(2):cmap_ind_max(2))');

set(get(colbar(3,1),'children'),'CData',(cmap_ind_min(1):cmap_ind_max(1))');
% set(get(colbar(3,2),'children'),'CData',(cmap_ind_min(2):cmap_ind_max(2))');

set(get(colbar(4,1),'children'),'CData',(cmap_ind_min(1):cmap_ind_max(1))');
set(get(colbar(4,2),'children'),'CData',(cmap_ind_min(2):cmap_ind_max(2))');

%RESET POSITIONS WITH NEW COLORBARS
pos1 = get(ax(4,1),'position');
set(ax(4,2),'position',pos1)

title('(d) \partial J/\partial q\prime (Jkg^{-1}kg^{-1}kg)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca, 'color', 'none');
box on; 
set(ax(4,1),'XTick',[],'Ytick',[]);
set(ax(4,2),'YTick',[10 24 38 52 66 80],'Xtick',[-100 -80 -60 -40 -20 0]);
set(ax(4,2),'XTickLabel',{'100W', '80W', '60W', '40W', '20W', '0'},'YtickLabel',{'10N' '24N' '38N' '52N' '66N' '80N'});

set(ax(4,1),'FontSize',fontsize,'FontName','TimesNewRoman')
set(ax(4,2),'FontSize',fontsize,'FontName','TimesNewRoman')

clim_new_1 = get(ax(1,1),'CLim');
clim_new_1_diff = clim_new_1(2)-clim_new_1(1);
clim_new_1_inc = (colbar1_ticks(1,1:length(get(colbar(1,1),'Ytick')))-clim_old_1(1,1))./(clim_old_1(2,1)-clim_old_1(1,1));
clim_new_1 = clim_new_1(1) + clim_new_1_inc*clim_new_1_diff;
set(colbar(1,1),'YTick',clim_new_1)

set(colbar(1,1),'YtickLabel',{'-4', '-3', '-2', '-1', '0', '1', '2', '3', '4'})
title(colbar(1,1),'\times 10^{-6}','FontSize',fontsize,'FontName','TimesNewRoman')

clim_new_1 = get(ax(2,1),'CLim');
clim_new_1_diff = clim_new_1(2)-clim_new_1(1);
clim_new_1_inc = (colbar1_ticks(2,1:length(get(colbar(2,1),'Ytick')))-clim_old_1(1,2))./(clim_old_1(2,2)-clim_old_1(1,2));
clim_new_1 = clim_new_1(1) + clim_new_1_inc*clim_new_1_diff;
set(colbar(2,1),'YTick',clim_new_1(1:end-1))

set(colbar(2,1),'YtickLabel',{'-2', '-1.5', '-1', '-0.5', '0', '0.5', '1.0', '1.5', '2.0'})
title(colbar(2,1),'\times 10^{-6}','FontSize',fontsize,'FontName','TimesNewRoman')

clim_new_1 = get(ax(3,1),'CLim');
clim_new_1_diff = clim_new_1(2)-clim_new_1(1);
clim_new_1_inc = (colbar1_ticks(3,1:length(colbar1))-clim_old_1(1,3))./(clim_old_1(2,3)-clim_old_1(1,3));
clim_new_1 = clim_new_1(1) + clim_new_1_inc*clim_new_1_diff;
set(colbar(3,1),'YTick',clim_new_1)

set(colbar(3,1),'YtickLabel',{'-1.5', '-1', '-0.5', '0', '0.5', '1.0', '1.5'})
title(colbar(3,1),'\times 10^{-5}','FontSize',fontsize,'FontName','TimesNewRoman')

clim_new_1 = get(ax(4,1),'CLim');
clim_new_1_diff = clim_new_1(2)-clim_new_1(1);
clim_new_1_inc = (colbar1_ticks(4,1:length(get(colbar(4,1),'Ytick')))-clim_old_1(1,4))./(clim_old_1(2,4)-clim_old_1(1,4));
clim_new_1 = clim_new_1(1) + clim_new_1_inc*clim_new_1_diff;
set(colbar(4,1),'YTick',clim_new_1)

set(colbar(4,1),'YtickLabel',{'-0.04', '-0.03', '-0.02', '-0.01', '0', '0.01', '0.02', '0.03', '0.04'})


clim_new_2 = get(ax(4,2),'CLim');
clim_new_2_diff = clim_new_2(2)-clim_new_2(1);
clim_new_2_inc = (colbar2_ticks-clim_old_2(1))./(clim_old_2(2)-clim_old_2(1));
clim_new_2 = clim_new_2(1) + clim_new_2_inc*clim_new_2_diff;
% set(colbar(1,2),'YTick',clim_new_2)
% set(colbar(2,2),'YTick',clim_new_2)
% set(colbar(3,2),'YTick',clim_new_2)
set(colbar(4,2),'YTick',clim_new_2)

poscbar2 = get(colbar(4,2),'Position');
set(colbar(4,2),'Position',[0.93 0.28 poscbar2(3) 0.45])
ylabel(colbar(4,2),'Convective Precipitation Rate (kgm^{-2}s^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')

tickstr = num2str(colbar2_ticks,1);
% set(colbar(1,2),'YtickLabel',{'2', '4', '6', '8', '10', '12'})
% title(colbar(1,2),'\times 10^{-4}','FontSize',fontsize,'FontName','TimesNewRoman')
% 
% set(colbar(2,2),'YtickLabel',{'2', '4', '6', '8', '10', '12'})
% title(colbar(2,2),'\times 10^{-4}','FontSize',fontsize,'FontName','TimesNewRoman')
% 
% 
% set(colbar(3,2),'YtickLabel',{'2', '4', '6', '8', '10', '12'})
% title(colbar(3,2),'\times 10^{-4}','FontSize',fontsize,'FontName','TimesNewRoman')


set(colbar(4,2),'YtickLabel',{'2', '4', '6', '8', '10', '12'})
title(colbar(4,2),'\times 10^{-4}','FontSize',fontsize,'FontName','TimesNewRoman')



saveas(gcf,'storm_sensitivity_mp_mn.eps', 'psc2')



