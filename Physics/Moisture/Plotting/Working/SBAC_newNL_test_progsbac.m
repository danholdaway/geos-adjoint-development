close all
clear
clc

load coast
coast_lat = lat; clear lat
coast_lon = long; clear long

count = 0;

qit_progoper_vert = zeros(72,24);
qlt_progoper_vert = zeros(72,24);
qit_diagsbac_vert = zeros(72,24);
qlt_diagsbac_vert = zeros(72,24);

loc_lon = 191;
loc_lat = 260;

for i = 0%:24
        
    count = count + 1;

    j = i;
    if i == 24
        j = 0;
    end
    
    v = [2013, 1, 7, 0, j, 0];

    a = datestr(v);

    if length(a) < 20
        h = '0000';
    else
        h = [a(16:17) a(19:20)];
    end

    disp(h)
    
    file = ['x0011dh_a.prog.eta.20130107_',h,'z.nc4'];
    if i == 24
        file = ['x0011dh_a.prog.eta.20130108_',h,'z.nc4'];
    end


    cd /discover/nobackup/drholdaw/SBAC/PROG_vs_DIAG/PROGSBAC

    lat = ncread(file,'lat');
    lon = ncread(file,'lon');
    
    qit_progoper = ncread(file,'QIT_old');
    qlt_progoper = ncread(file,'QLT_old');

    qit_diagsbac = ncread(file,'QIT_dh');
    qlt_diagsbac = ncread(file,'QLT_dh');
    

    cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/

    qit_progoper_vert(:,count) = qit_progoper(loc_lon,loc_lat,:);
    qlt_progoper_vert(:,count) = qlt_progoper(loc_lon,loc_lat,:);
    qit_diagsbac_vert(:,count) = qit_diagsbac(loc_lon,loc_lat,:);
    qlt_diagsbac_vert(:,count) = qlt_diagsbac(loc_lon,loc_lat,:);
        
    
end


ylim_top = 40;
cint = 10;

% 
% qit_vert = [qit_progoper_vert; qit_diagsbac_vert];
% qlt_vert = [qlt_progoper_vert; qlt_diagsbac_vert];
% 
% 
% qitv_max = max(qit_vert(:));
% qitv_min = min(qit_vert(:));
% qltv_max = max(qlt_vert(:));
% qltv_min = min(qlt_vert(:));
% 
% 
% cint_i = (qitv_max-qitv_min)/cint;
% cint_l = (qltv_max-qltv_min)/cint;
% 
% cint_i_po = cint_i;
% cint_l_po = cint_l;
% 
% 
% 
% figure
% set(gcf,'position',[1312 123 560 854])
% 
% subplot(2,1,1)
% contour(qit_progoper_vert,'LevelStep',cint_i_po)
% colorbar
% set(gca,'YDir','reverse')
% ylim([ylim_top 72])
% caxis([qitv_min qitv_max])
% 
% subplot(2,1,2)
% contour(qit_diagsbac_vert,'LevelStep',cint_i)
% colorbar
% set(gca,'YDir','reverse')
% ylim([ylim_top 72])
% caxis([qitv_min qitv_max])
% 
% 
% figure
% set(gcf,'position',[1977 85 560 854])
% 
% subplot(2,1,1)
% contour(qlt_progoper_vert,'LevelStep',cint_l_po)
% colorbar
% set(gca,'YDir','reverse')
% ylim([ylim_top 72])
% caxis([qltv_min qltv_max])
% 
% 
% subplot(2,1,2)
% contour(qlt_diagsbac_vert,'LevelStep',cint_l)
% colorbar
% set(gca,'YDir','reverse')
% ylim([ylim_top 72])
% caxis([qltv_min qltv_max])





plot_level_i = 50;
plot_level_l = 60;



qit_lev = [qit_progoper(:,:,plot_level_i) qit_diagsbac(:,:,plot_level_i)];
qlt_lev = [qlt_progoper(:,:,plot_level_l) qlt_diagsbac(:,:,plot_level_l)];


qit_max = max(abs(qit_lev(:)));
qit_min = min(abs(qit_lev(:)));
qlt_max = max(abs(qlt_lev(:)));
qlt_min = min(abs(qlt_lev(:)));



cint_i = (qit_max-qit_min)/cint;
cint_l = (qlt_max-qlt_min)/cint;

cint_i_po = cint_i;
cint_l_po = cint_l;




figure
set(gcf,'position',[128    65   560   854])

subplot(3,1,1)
contour(lon,lat,qit_progoper(:,:,plot_level_i)','LevelStep',cint_i_po)
hold on
plot(coast_lon,coast_lat,'Color',[0.7 0.7 0.7])
colorbar
caxis([qit_min qit_max])
    
subplot(3,1,2)
contour(lon,lat,qit_diagsbac(:,:,plot_level_i)','LevelStep',cint_i)
hold on
plot(coast_lon,coast_lat,'Color',[0.7 0.7 0.7])
colorbar
caxis([qit_min qit_max])

subplot(3,1,3)
contour(lon,lat,qit_diagsbac(:,:,plot_level_i)'-qit_progoper(:,:,plot_level_i)')
hold on
plot(coast_lon,coast_lat,'Color',[0.7 0.7 0.7])
colorbar
% caxis([qit_min qit_max])


figure
set(gcf,'position',[719   65   560   854])

subplot(3,1,1)
contour(lon,lat,qlt_progoper(:,:,plot_level_l)','LevelStep',cint_l_po)
hold on
plot(coast_lon,coast_lat,'Color',[0.7 0.7 0.7])
colorbar
caxis([qlt_min qlt_max])

subplot(3,1,2)
contour(lon,lat,qlt_diagsbac(:,:,plot_level_l)','LevelStep',cint_l)
hold on
plot(coast_lon,coast_lat,'Color',[0.7 0.7 0.7])
colorbar
caxis([qlt_min qlt_max])

subplot(3,1,3)
contour(lon,lat,qlt_diagsbac(:,:,plot_level_i)'-qlt_progoper(:,:,plot_level_i)')
hold on
plot(coast_lon,coast_lat,'Color',[0.7 0.7 0.7])
colorbar
% caxis([qlt_min qlt_max])

