close all
clear
clc

cd /discover/nobackup/drholdaw/tmp.iau590/tmp.8823/sens.20120318.000000/


cd /discover/nobackup/drholdaw/tmp.iau590/tmp.8823/nlm_runs/free/

U_nlm_free  = ncread('iau590.prog.eta.20120317_12z.nc4','u');
V_nlm_free  = ncread('iau590.prog.eta.20120317_12z.nc4','v');   
TV_nlm_free = ncread('iau590.prog.eta.20120317_12z.nc4','tv');
DP_nlm_free = ncread('iau590.prog.eta.20120317_12z.nc4','delp');
QV_nlm_free = ncread('iau590.prog.eta.20120317_12z.nc4','sphu');
QL_nlm_free = ncread('iau590.prog.eta.20120317_12z.nc4','qltot');
QI_nlm_free = ncread('iau590.prog.eta.20120317_12z.nc4','qitot');
O3_nlm_free = ncread('iau590.prog.eta.20120317_12z.nc4','ozone');

cd /discover/nobackup/drholdaw/tmp.iau590/tmp.8823/nlm_runs/del_0.01/

U_nlm_del_0_01  = ncread('iau590.prog.eta.20120317_12z.nc4','u');
V_nlm_del_0_01  = ncread('iau590.prog.eta.20120317_12z.nc4','v');   
TV_nlm_del_0_01 = ncread('iau590.prog.eta.20120317_12z.nc4','tv');
DP_nlm_del_0_01 = ncread('iau590.prog.eta.20120317_12z.nc4','delp');
QV_nlm_del_0_01 = ncread('iau590.prog.eta.20120317_12z.nc4','sphu');
QL_nlm_del_0_01 = ncread('iau590.prog.eta.20120317_12z.nc4','qltot');
QI_nlm_del_0_01 = ncread('iau590.prog.eta.20120317_12z.nc4','qitot');
O3_nlm_del_0_01 = ncread('iau590.prog.eta.20120317_12z.nc4','ozone');

cd /discover/nobackup/drholdaw/tmp.iau590/tmp.8823/nlm_runs/del_0.1/

U_nlm_del_0_1  = ncread('iau590.prog.eta.20120317_12z.nc4','u');
V_nlm_del_0_1  = ncread('iau590.prog.eta.20120317_12z.nc4','v');   
TV_nlm_del_0_1 = ncread('iau590.prog.eta.20120317_12z.nc4','tv');
DP_nlm_del_0_1 = ncread('iau590.prog.eta.20120317_12z.nc4','delp');
QV_nlm_del_0_1 = ncread('iau590.prog.eta.20120317_12z.nc4','sphu');
QL_nlm_del_0_1 = ncread('iau590.prog.eta.20120317_12z.nc4','qltot');
QI_nlm_del_0_1 = ncread('iau590.prog.eta.20120317_12z.nc4','qitot');
O3_nlm_del_0_1 = ncread('iau590.prog.eta.20120317_12z.nc4','ozone');

cd /discover/nobackup/drholdaw/tmp.iau590/tmp.8823/nlm_runs/del_0.2/

U_nlm_del_0_2  = ncread('iau590.prog.eta.20120317_12z.nc4','u');
V_nlm_del_0_2  = ncread('iau590.prog.eta.20120317_12z.nc4','v');   
TV_nlm_del_0_2 = ncread('iau590.prog.eta.20120317_12z.nc4','tv');
DP_nlm_del_0_2 = ncread('iau590.prog.eta.20120317_12z.nc4','delp');
QV_nlm_del_0_2 = ncread('iau590.prog.eta.20120317_12z.nc4','sphu');
QL_nlm_del_0_2 = ncread('iau590.prog.eta.20120317_12z.nc4','qltot');
QI_nlm_del_0_2 = ncread('iau590.prog.eta.20120317_12z.nc4','qitot');
O3_nlm_del_0_2 = ncread('iau590.prog.eta.20120317_12z.nc4','ozone');

cd /discover/nobackup/drholdaw/tmp.iau590/tmp.8823/nlm_runs/del_0.5/

U_nlm_del_0_5  = ncread('iau590.prog.eta.20120317_12z.nc4','u');
V_nlm_del_0_5  = ncread('iau590.prog.eta.20120317_12z.nc4','v');   
TV_nlm_del_0_5 = ncread('iau590.prog.eta.20120317_12z.nc4','tv');
DP_nlm_del_0_5 = ncread('iau590.prog.eta.20120317_12z.nc4','delp');
QV_nlm_del_0_5 = ncread('iau590.prog.eta.20120317_12z.nc4','sphu');
QL_nlm_del_0_5 = ncread('iau590.prog.eta.20120317_12z.nc4','qltot');
QI_nlm_del_0_5 = ncread('iau590.prog.eta.20120317_12z.nc4','qitot');
O3_nlm_del_0_5 = ncread('iau590.prog.eta.20120317_12z.nc4','ozone');

cd /discover/nobackup/drholdaw/tmp.iau590/tmp.8823/nlm_runs/del_1.0/

U_nlm_del_1_0  = ncread('iau590.prog.eta.20120317_12z.nc4','u');
V_nlm_del_1_0  = ncread('iau590.prog.eta.20120317_12z.nc4','v');   
TV_nlm_del_1_0 = ncread('iau590.prog.eta.20120317_12z.nc4','tv');
DP_nlm_del_1_0 = ncread('iau590.prog.eta.20120317_12z.nc4','delp');
QV_nlm_del_1_0 = ncread('iau590.prog.eta.20120317_12z.nc4','sphu');
QL_nlm_del_1_0 = ncread('iau590.prog.eta.20120317_12z.nc4','qltot');
QI_nlm_del_1_0 = ncread('iau590.prog.eta.20120317_12z.nc4','qitot');
O3_nlm_del_1_0 = ncread('iau590.prog.eta.20120317_12z.nc4','ozone');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/

plot_level = 50;

figure
subplot(3,2,1)
contourf(TV_nlm_del_0_01(:,:,plot_level)' - TV_nlm_free(:,:,plot_level)','LineStyle','none')
title('F[x+0.01(x_a-x_b),12h] - F[x,12h]')
colorbar

subplot(3,2,2)
contourf(TV_nlm_del_0_1(:,:,plot_level)' - TV_nlm_free(:,:,plot_level)','LineStyle','none')
title('F[x+0.1(x_a-x_b),12h] - F[x,12h]')
colorbar

subplot(3,2,3)
contourf(TV_nlm_del_0_2(:,:,plot_level)' - TV_nlm_free(:,:,plot_level)','LineStyle','none')
title('F[x+0.2(x_a-x_b),12h] - F[x,12h]')
colorbar

subplot(3,2,4)
contourf(TV_nlm_del_0_5(:,:,plot_level)' - TV_nlm_free(:,:,plot_level)','LineStyle','none')
title('F[x+0.5(x_a-x_b),12h] - F[x,12h]')
colorbar

subplot(3,2,5)
contourf(TV_nlm_del_1_0(:,:,plot_level)' - TV_nlm_free(:,:,plot_level)','LineStyle','none')
title('F[x+1.0(x_a-x_b),12h] - F[x,12h]')
colorbar