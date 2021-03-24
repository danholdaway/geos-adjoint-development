close all
clear
clc

cd /discover/nobackup/drholdaw/tmp.31275/sens.20120318.000000/

U_control  = ncread('fsens.eta.moi.control.20120317_06z.nc4','u');
V_control  = ncread('fsens.eta.moi.control.20120317_06z.nc4','v');
TV_control = ncread('fsens.eta.moi.control.20120317_06z.nc4','tv');
DP_control = ncread('fsens.eta.moi.control.20120317_06z.nc4','delp');
QV_control = ncread('fsens.eta.moi.control.20120317_06z.nc4','sphu');
QL_control = ncread('fsens.eta.moi.control.20120317_06z.nc4','qltot');
QI_control = ncread('fsens.eta.moi.control.20120317_06z.nc4','qitot');
O3_control = ncread('fsens.eta.moi.control.20120317_06z.nc4','ozone');

U_new  = ncread('fsens.eta.moi.new.20120317_06z.nc4','u');
V_new  = ncread('fsens.eta.moi.new.20120317_06z.nc4','v');
TV_new = ncread('fsens.eta.moi.new.20120317_06z.nc4','tv');
DP_new = ncread('fsens.eta.moi.new.20120317_06z.nc4','delp');
QV_new = ncread('fsens.eta.moi.new.20120317_06z.nc4','sphu');
QL_new = ncread('fsens.eta.moi.new.20120317_06z.nc4','qltot');
QI_new = ncread('fsens.eta.moi.new.20120317_06z.nc4','qitot');
O3_new = ncread('fsens.eta.moi.new.20120317_06z.nc4','ozone');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/

plot_level = 60;

scrsz = get(0,'ScreenSize');
figure('visible','on','Position',[0.3*scrsz(3) scrsz(4) (0.3)*scrsz(3) (0.5)*scrsz(4)])

subplot(3,2,1)
contourf(U_control(:,:,plot_level)')
colorbar

subplot(3,2,2)
contourf(V_control(:,:,plot_level)')
colorbar

subplot(3,2,3)
contourf(TV_control(:,:,plot_level)')
colorbar

subplot(3,2,4)
contourf(QV_control(:,:,plot_level)')
colorbar

subplot(3,2,5)
contourf(DP_control(:,:,plot_level)')
colorbar

scrsz = get(0,'ScreenSize');
figure('visible','on','Position',[0.3*scrsz(3) scrsz(4) (0.3)*scrsz(3) (0.5)*scrsz(4)])

subplot(3,2,1)
contourf(U_new(:,:,plot_level)')
colorbar

subplot(3,2,2)
contourf(V_new(:,:,plot_level)')
colorbar

subplot(3,2,3)
contourf(TV_new(:,:,plot_level)')
colorbar

subplot(3,2,4)
contourf(QV_new(:,:,plot_level)')
colorbar

subplot(3,2,5)
contourf(DP_new(:,:,plot_level)')
colorbar

scrsz = get(0,'ScreenSize');
figure('visible','on','Position',[0.3*scrsz(3) scrsz(4) (0.3)*scrsz(3) (0.5)*scrsz(4)])

subplot(3,2,1)
contourf(U_control(:,:,plot_level)' - U_new(:,:,plot_level)')
colorbar

subplot(3,2,2)
contourf(V_control(:,:,plot_level)' - V_new(:,:,plot_level)')
colorbar

subplot(3,2,3)
contourf(TV_control(:,:,plot_level)' - TV_new(:,:,plot_level)')
colorbar

subplot(3,2,4)
contourf(QV_control(:,:,plot_level)' - QV_new(:,:,plot_level)')
colorbar

subplot(3,2,5)
contourf(DP_control(:,:,plot_level)' - DP_new(:,:,plot_level)')
colorbar

U_rel_error  = (U_control(:,:,:)  - U_new(:,:,:))./ U_new(:,:,:);
V_rel_error  = (V_control(:,:,:)  - V_new(:,:,:))./ V_new(:,:,:);
TV_rel_error = (TV_control(:,:,:) - TV_new(:,:,:))./TV_new(:,:,:);
QV_rel_error = (QV_control(:,:,:) - QV_new(:,:,:))./QV_new(:,:,:);
DP_rel_error = (DP_control(:,:,:) - DP_new(:,:,:))./DP_new(:,:,:);

[num1 idx] = max(U_rel_error(:));
[x_u y_u z_u] = ind2sub(size(U_rel_error),idx);
disp([num1 x_u y_u z_u])

[num2 idx] = max(V_rel_error(:));
[x_v y_v z_v] = ind2sub(size(V_rel_error),idx);
disp([num2 x_v y_v z_v])

[num3 idx] = max(TV_rel_error(:));
[x_t y_t z_t] = ind2sub(size(TV_rel_error),idx);
disp([num3 x_t y_t z_t])

[num4 idx] = max(QV_rel_error(:));
[x_q y_q z_q] = ind2sub(size(QV_rel_error),idx);
disp([num4 x_q y_q z_q])

[num5 idx] = max(DP_rel_error(:));
[x_d y_d z_d] = ind2sub(size(DP_rel_error),idx);
disp([num5 x_d y_d z_d])
