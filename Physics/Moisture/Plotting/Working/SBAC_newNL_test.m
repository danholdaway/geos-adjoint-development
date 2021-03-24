close all
clear
clc

load coast
coast_lat = lat; clear lat
coast_lon = long; clear long

cmap = colormap;
cmap(32,:) = [1 1 1];
cmap(33,:) = [1 1 1];
close

count = 0;

qit_1_vert = zeros(72,24);
qlt_1_vert = zeros(72,24);
qit_2_vert = zeros(72,24);
qlt_2_vert = zeros(72,24);
qit_diagsbac_vert = zeros(72,24);
qlt_diagsbac_vert = zeros(72,24);

loc_lon = 191;
loc_lat = 260;

for i = 18%:24
        
    count = count + 1;

    j = i;
    if i == 24
        j = 0;
    end
    
    v = [2013, 1, 7, 0, j, 0];

    a = datestr(v);

    if length(a) < 20
        h = '00';
    else
        h = [a(16:17)];
    end

    disp(h)
    
    file = ['x0011dh_a.prog.eta.20130116_',h,'z.nc4'];
    if i == 24
        file = ['x0011dh_a.prog.eta.20130107_',h,'z.nc4'];
    end


    cd /discover/nobackup/drholdaw/SBAC/ORIG_vs_NEW/d_ORIG_LSCLOUD/
    cd /discover/nobackup/drholdaw/tmp.22292/

    lat = ncread(file,'lat');
    lon = ncread(file,'lon');
    
    qit_1 = ncread(file,'qitot');
    qlt_1 = ncread(file,'qltot');
    
    cd /discover/nobackup/drholdaw/SBAC/ORIG_vs_NEW/d_SBAC_LSCLOUD/
    cd /discover/nobackup/drholdaw/tmp.33392/
    
    qit_2 = ncread(file,'qitot');
    qlt_2 = ncread(file,'qltot');
    
        
    cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/

            
    
end

[maxival location] = max(abs(qit_1(:)-qit_2(:)));
[rmax,cmax,lmax] = ind2sub(size(qit_1),location);
fprintf('Max qi difference is %f \n',maxival)
fprintf('Located at (%d,%d,%d) \n',rmax,cmax,lmax)
fprintf('qi_a = %f, qi_b = %f \n\n',qit_1(rmax,cmax,lmax),qit_2(rmax,cmax,lmax))

[maxival location] = max(abs(qlt_1(:)-qlt_2(:)));
[rmax,cmax,lmax] = ind2sub(size(qlt_1),location);
fprintf('Max ql difference is %f \n',maxival)
fprintf('Located at (%d,%d,%d) \n',rmax,cmax,lmax)
fprintf('ql_a = %f, ql_b = %f \n\n',qlt_1(rmax,cmax,lmax),qlt_2(rmax,cmax,lmax))

plot_level_i = 50;
plot_level_l = 50;
cint = 30;


qit_1_lev = qit_1(:,:,plot_level_i);
qlt_1_lev = qlt_1(:,:,plot_level_l);

qit_2_lev = qit_2(:,:,plot_level_i);
qlt_2_lev = qlt_2(:,:,plot_level_l);

qit_error = qit_2_lev - qit_1_lev;
qlt_error = qlt_2_lev - qlt_1_lev;

qit_all_lev = [ qit_1_lev qit_2_lev ];
qlt_all_lev = [ qlt_1_lev qlt_2_lev ];


qit_all_lev_max = max(qit_all_lev(:));
qit_all_lev_min = min(qit_all_lev(:));
qlt_all_lev_max = max(qlt_all_lev(:));
qlt_all_lev_min = min(qlt_all_lev(:));

qi_error_max = max(abs(qit_error(:)));
ql_error_max = max(abs(qlt_error(:)));

for i = 1:length(lon);
    for j = 1:length(lat);
        
        if abs(qit_error(i,j)) < 0.001*qi_error_max
            
            qit_error(i,j) = 0;
            
        end
        if abs(qlt_error(i,j)) < 0.001*ql_error_max
            
            qlt_error(i,j) = 0;
            
        end
        
    end
end

cint_i = (qit_all_lev_max - qit_all_lev_min)/cint;
cint_l = (qlt_all_lev_max - qlt_all_lev_min)/cint;



figure
set(gcf,'position',[128    65   560   854])

subplot(3,1,1)
contour(lon,lat,qit_1_lev','LevelStep',cint_i)
hold on
plot(coast_lon,coast_lat,'Color',[0.7 0.7 0.7])
colorbar
caxis([qit_all_lev_min qit_all_lev_max])
    
subplot(3,1,2)
contour(lon,lat,qit_2_lev','LevelStep',cint_i)
hold on
plot(coast_lon,coast_lat,'Color',[0.7 0.7 0.7])
colorbar
caxis([qit_all_lev_min qit_all_lev_max])

subplot(3,1,3)
contour(lon,lat,qit_error')
hold on
title('Middle - Top')
plot(coast_lon,coast_lat,'Color',[0.7 0.7 0.7])
colorbar
caxis([-qi_error_max qi_error_max])



figure
set(gcf,'position',[719   65   560   854])

subplot(3,1,1)
contour(lon,lat,qlt_1_lev','LevelStep',cint_l)
hold on
plot(coast_lon,coast_lat,'Color',[0.7 0.7 0.7])
colorbar
caxis([qlt_all_lev_min qlt_all_lev_max])

subplot(3,1,2)
contour(lon,lat,qlt_2_lev','LevelStep',cint_l)
hold on
plot(coast_lon,coast_lat,'Color',[0.7 0.7 0.7])
colorbar
caxis([qlt_all_lev_min qlt_all_lev_max])

subplot(3,1,3)
contour(lon,lat,qlt_error')
hold on
title('Middle - Top')
plot(coast_lon,coast_lat,'Color',[0.7 0.7 0.7])
colorbar
caxis([-ql_error_max ql_error_max])

