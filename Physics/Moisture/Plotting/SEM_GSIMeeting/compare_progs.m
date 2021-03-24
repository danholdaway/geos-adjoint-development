close all
clear
clc

load coast
coast_lat = lat; clear lat
coast_lon = long; clear long


%Choose model level to plot.
plot_level = 63;

fontsize = 11;

line_wid_cont = 1.0;
line_wid_det = 0.6;

grey = 0.75;

%Choose whether to match the scale.
match_scale = 1;

%Make center of colormap white
cmap = colormap;
cmap(32,:) = [1 1 1];
cmap(33,:) = [1 1 1];
close

cd /discover/nobackup/drholdaw/tmp.22292a/prog/prog_free/
file = 'x0011dh_a.prog.eta.20130117_00z.nc4';

lon = ncread(file,'lon');
lat = ncread(file,'lat');

u_free = ncread(file,'u');
v_free = ncread(file,'v');
t_free = ncread(file,'tv');
q_free = ncread(file,'sphu');
p_free = ncread(file,'delp');
qi_free = ncread(file,'qitot');
ql_free = ncread(file,'qltot');
o3_free = ncread(file,'ozone');

cd /discover/nobackup/drholdaw/tmp.22292a/prog/prog_replay10/
file = 'x0011dh_a.prog.eta.20130117_00z.nc4';

u_r10 = ncread(file,'u');
v_r10 = ncread(file,'v');
t_r10 = ncread(file,'tv');
q_r10 = ncread(file,'sphu');
p_r10 = ncread(file,'delp');
qi_r10 = ncread(file,'qitot');
ql_r10 = ncread(file,'qltot');
o3_r10 = ncread(file,'ozone');

cd /discover/nobackup/drholdaw/tmp.22292a/prog/prog_replay05/
file = 'x0011dh_a.prog.eta.20130117_00z.nc4';

u_r05 = ncread(file,'u');
v_r05 = ncread(file,'v');
t_r05 = ncread(file,'tv');
q_r05 = ncread(file,'sphu');
p_r05 = ncread(file,'delp');
qi_r05 = ncread(file,'qitot');
ql_r05 = ncread(file,'qltot');
o3_r05 = ncread(file,'ozone');

cd /discover/nobackup/drholdaw/tmp.22292a/prog/prog_replay01/
file = 'x0011dh_a.prog.eta.20130117_00z.nc4';

u_r01 = ncread(file,'u');
v_r01 = ncread(file,'v');
t_r01 = ncread(file,'tv');
q_r01 = ncread(file,'sphu');
p_r01 = ncread(file,'delp');
qi_r01 = ncread(file,'qitot');
ql_r01 = ncread(file,'qltot');
o3_r01 = ncread(file,'ozone');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/GSI' Talk'/

u_10 = u_r10 - u_free;
v_10 = v_r10 - v_free;
t_10 = t_r10 - t_free;
q_10 = q_r10 - q_free;
p_10 = p_r10 - p_free;
qi_10 = qi_r10 - qi_free;
ql_10 = ql_r10 - ql_free;
o3_10 = o3_r10 - o3_free;

u_05 = u_r05 - u_free;
v_05 = v_r05 - v_free;
t_05 = t_r05 - t_free;
q_05 = q_r05 - q_free;
p_05 = p_r05 - p_free;
qi_05 = qi_r05 - qi_free;
ql_05 = ql_r05 - ql_free;
o3_05 = o3_r05 - o3_free;

u_01 = u_r01 - u_free;
v_01 = v_r01 - v_free;
t_01 = t_r01 - t_free;
q_01 = q_r01 - q_free;
p_01 = p_r01 - p_free;
qi_01 = qi_r01 - qi_free;
ql_01 = ql_r01 - ql_free;
o3_01 = o3_r01 - o3_free;


u_a_lev = u_10(:,:,plot_level);
v_a_lev = v_10(:,:,plot_level);
t_a_lev = t_10(:,:,plot_level);
q_a_lev = q_10(:,:,plot_level);
p_a_lev = p_10(:,:,plot_level);
qi_a_lev = qi_10(:,:,plot_level);
ql_a_lev = ql_10(:,:,plot_level);
o3_a_lev = o3_10(:,:,plot_level);

u_b_lev = u_05(:,:,plot_level) * 2;
v_b_lev = v_05(:,:,plot_level) * 2;
t_b_lev = t_05(:,:,plot_level) * 2;
q_b_lev = q_05(:,:,plot_level) * 2;
p_b_lev = p_05(:,:,plot_level) * 2;
qi_b_lev = qi_05(:,:,plot_level) * 2;
ql_b_lev = ql_05(:,:,plot_level) * 2;
o3_b_lev = o3_05(:,:,plot_level) * 2;

u_c_lev = u_01(:,:,plot_level) * 10;
v_c_lev = v_01(:,:,plot_level) * 10;
t_c_lev = t_01(:,:,plot_level) * 10;
q_c_lev = q_01(:,:,plot_level) * 10;
p_c_lev = p_01(:,:,plot_level) * 10;
qi_c_lev = qi_01(:,:,plot_level) * 10;
ql_c_lev = ql_01(:,:,plot_level) * 10;
o3_c_lev = o3_01(:,:,plot_level) * 10;

%Get Maximums for Levels
u_lev = [u_a_lev u_b_lev u_c_lev];
maxu = max(abs(u_lev(:)));
if min(u_lev(:)) >= 0
    minu = 0;
else
    minu = -maxu;
end
v_lev = [v_a_lev v_b_lev v_c_lev];
maxv = max(abs(v_lev(:)));
if min(v_lev(:)) >= 0
    minv = min(v_lev(:));
else
    minv = -maxv;
end
t_lev = [t_a_lev t_b_lev t_c_lev];
maxt = max(abs(t_lev(:)));
if min(t_lev(:)) >= 0
    mint = min(t_lev(:));
else
    mint = -maxt;
end
q_lev = [q_a_lev q_b_lev q_c_lev];
maxq = max(abs(q_lev(:)));
if min(q_lev(:)) >= 0
    minq = min(q_lev(:));
else
    minq = -maxq;
end
p_lev = [p_a_lev p_b_lev p_c_lev];
maxp = max(abs(p_lev(:)));
if min(p_lev(:)) >= 0
    minp = min(p_lev(:));
else
    minp = -maxp;
end
qi_lev = [qi_a_lev qi_b_lev qi_c_lev];
maxqi = max(abs(qi_lev(:)));
if min(qi_lev(:)) >= 0
    minqi = min(qi_lev(:));
else
    minqi = -maxqi;
end
ql_lev = [ql_a_lev ql_b_lev ql_c_lev];
maxql = max(abs(ql_lev(:)));
if min(ql_lev(:)) >= 0
    minql = min(ql_lev(:));
else
    minql = -maxql;
end
o3_lev = [o3_a_lev o3_b_lev o3_c_lev];
maxo3 = max(abs(o3_lev(:)));
if min(o3_lev(:)) >= 0
    mino3 = min(o3_lev(:));
else
    mino3 = -maxo3;
end





for i = 1:length(lon)
    for j = 1:length(lat)
        
        if abs(q_b_lev(i,j)) >= max(abs(q_a_lev(i,j)))
            
            q_b_lev(i,j) = max(abs(q_a_lev(i,j))) * abs(q_b_lev(i,j))/q_b_lev(i,j);
            
        end
        
        if abs(q_c_lev(i,j)) >= max(abs(q_a_lev(i,j)))
            
            q_c_lev(i,j) = max(abs(q_a_lev(i,j))) * abs(q_c_lev(i,j))/q_c_lev(i,j);
            
        end
        
    end
end

cint = 2*maxq/30;


figure
set(gcf,'position',[130 250 1149 669])

subplot(2,2,1)

[C,h] = contour(lon,lat,((q_a_lev))');
% caxis([minq maxq])
caxis([-max(abs(q_a_lev(:))) max(abs(q_a_lev(:)))])
title('q')
colormap(cmap)
colbar = colorbar;
% set(h,'LevelStep',cint);
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
title('\delta x = (x_a - x_b)','FontSize',fontsize,'FontName','TimesNewRoman')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])

subplot(2,2,2)

[C,h] = contour(lon,lat,((q_b_lev))');
% caxis([minq maxq])
caxis([-max(abs(q_a_lev(:))) max(abs(q_a_lev(:)))])
title('q')
colormap(cmap)
colbar = colorbar;
% set(h,'LevelStep',cint);
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
title('\delta x = 0.5(x_a - x_b)','FontSize',fontsize,'FontName','TimesNewRoman')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])


subplot(2,2,3)

[C,h] = contour(lon,lat,((q_c_lev))');
% caxis([minq maxq])
caxis([-max(abs(q_a_lev(:))) max(abs(q_a_lev(:)))])
title('q')
colormap(cmap)
colbar = colorbar;
% set(h,'LevelStep',2*max(abs(q_c_lev(:)))/20);
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
title('\delta x = 0.1(x_a - x_b)','FontSize',fontsize,'FontName','TimesNewRoman')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
set(gca,'position',[0.3500    0.1100    0.2781    0.3412])




for i = 1:length(lon)
    for j = 1:length(lat)
        
        if abs(ql_b_lev(i,j)) >= max(abs(ql_a_lev(i,j)))
            
            ql_b_lev(i,j) = max(abs(ql_a_lev(i,j))) * abs(ql_b_lev(i,j))/ql_b_lev(i,j);
            
        end
        
        if abs(ql_c_lev(i,j)) >= max(abs(ql_a_lev(i,j)))
            
            ql_c_lev(i,j) = max(abs(ql_a_lev(i,j))) * abs(ql_c_lev(i,j))/ql_c_lev(i,j);
            
        end
        
    end
end

cint = 2*maxql/30;


figure
set(gcf,'position',[130 250 1149 669])

subplot(2,2,1)

[C,h] = contour(lon,lat,((ql_a_lev))');
% caxis([minq maxq])
caxis([-max(abs(ql_a_lev(:))) max(abs(ql_a_lev(:)))])
title('q')
colormap(cmap)
colbar = colorbar;
% set(h,'LevelStep',cint);
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
title('\delta x = (x_a - x_b)','FontSize',fontsize,'FontName','TimesNewRoman')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])

subplot(2,2,2)

[C,h] = contour(lon,lat,((ql_b_lev))');
% caxis([minq maxq])
caxis([-max(abs(ql_a_lev(:))) max(abs(ql_a_lev(:)))])
title('q')
colormap(cmap)
colbar = colorbar;
% set(h,'LevelStep',cint);
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
title('\delta x = 0.5(x_a - x_b)','FontSize',fontsize,'FontName','TimesNewRoman')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])


subplot(2,2,3)

[C,h] = contour(lon,lat,((ql_c_lev))');
% caxis([minq maxq])
caxis([-max(abs(ql_a_lev(:))) max(abs(ql_a_lev(:)))])
title('q')
colormap(cmap)
colbar = colorbar;
% set(h,'LevelStep',2*max(abs(q_c_lev(:)))/20);
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
title('\delta x = 0.1(x_a - x_b)','FontSize',fontsize,'FontName','TimesNewRoman')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
set(gca,'position',[0.3500    0.1100    0.2781    0.3412])






% for i = 1:length(lon)
%     for j = 1:length(lat)
%         
%         if abs(u_b_lev(i,j)) >= max(abs(u_a_lev(i,j)))
%             
%             u_b_lev(i,j) = max(abs(u_a_lev(i,j))) * abs(u_b_lev(i,j))/u_b_lev(i,j);
%             
%         end
%         
%         if abs(u_c_lev(i,j)) >= max(abs(u_a_lev(i,j)))
%             
%             u_c_lev(i,j) = max(abs(u_a_lev(i,j))) * abs(u_c_lev(i,j))/u_c_lev(i,j);
%             
%         end
%         
%     end
% end

cint = 2*maxql/30;


figure
set(gcf,'position',[130 250 1149 669])

subplot(2,2,1)

[C,h] = contour(lon,lat,((u_a_lev))');
% caxis([minq maxq])
% caxis([-max(abs(u_a_lev(:))) max(abs(u_a_lev(:)))])
title('q')
colormap(cmap)
colbar = colorbar;
% set(h,'LevelStep',cint);
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
title('\delta x = (x_a - x_b)','FontSize',fontsize,'FontName','TimesNewRoman')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])

subplot(2,2,2)

[C,h] = contour(lon,lat,((u_b_lev))');
% caxis([minq maxq])
% caxis([-max(abs(u_a_lev(:))) max(abs(u_a_lev(:)))])
title('q')
colormap(cmap)
colbar = colorbar;
% set(h,'LevelStep',cint);
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
title('\delta x = 0.5(x_a - x_b)','FontSize',fontsize,'FontName','TimesNewRoman')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])


subplot(2,2,3)

[C,h] = contour(lon,lat,((u_c_lev))');
% caxis([minq maxq])
% caxis([-max(abs(u_a_lev(:))) max(abs(u_a_lev(:)))])
title('q')
colormap(cmap)
colbar = colorbar;
% set(h,'LevelStep',2*max(abs(q_c_lev(:)))/20);
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
title('\delta x = 0.1(x_a - x_b)','FontSize',fontsize,'FontName','TimesNewRoman')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
set(gca,'position',[0.3500    0.1100    0.2781    0.3412])

