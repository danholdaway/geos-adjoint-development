close all
clear
clc

% cd /discover/nobackup/drholdaw/ExperimentData/New_akbk/
% 
% u_old = ncread('x0009b.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00zBACK.nc4','u');
% v_old = ncread('x0009b.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00zBACK.nc4','v');
% t_old = ncread('x0009b.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00zBACK.nc4','tv');
% q_old = ncread('x0009b.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00zBACK.nc4','sphu');
% p_old = ncread('x0009b.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00zBACK.nc4','delp');
% oz_old = ncread('x0009b.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00zBACK.nc4','ozone');
% qi_old = ncread('x0009b.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00zBACK.nc4','qitot');
% ql_old = ncread('x0009b.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00zBACK.nc4','qltot');
% 
% u_new = ncread('x0009b.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','u');
% v_new = ncread('x0009b.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','v');
% t_new = ncread('x0009b.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','tv');
% q_new = ncread('x0009b.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','sphu');
% p_new = ncread('x0009b.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','delp');
% oz_new = ncread('x0009b.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','ozone');
% qi_new = ncread('x0009b.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','qitot');
% ql_new = ncread('x0009b.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','qltot');
% 
% cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/
% 
% cd /discover/nobackup/drholdaw/x0009a/asens
% 
% u_old = ncread('x0009a.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00zBACK.nc4','u');
% v_old = ncread('x0009a.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00zBACK.nc4','v');
% t_old = ncread('x0009a.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00zBACK.nc4','tv');
% q_old = ncread('x0009a.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00zBACK.nc4','sphu');
% p_old = ncread('x0009a.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00zBACK.nc4','delp');
% oz_old = ncread('x0009a.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00zBACK.nc4','ozone');
% qi_old = ncread('x0009a.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00zBACK.nc4','qitot');
% ql_old = ncread('x0009a.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00zBACK.nc4','qltot');

% cd /discover/nobackup/drholdaw/x0009b/asens
% 
% u_old = ncread('x0009b.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00zBACK.nc4','u');
% v_old = ncread('x0009b.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00zBACK.nc4','v');
% t_old = ncread('x0009b.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00zBACK.nc4','tv');
% q_old = ncread('x0009b.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00zBACK.nc4','sphu');
% p_old = ncread('x0009b.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00zBACK.nc4','delp');
% oz_old = ncread('x0009b.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00zBACK.nc4','ozone');
% qi_old = ncread('x0009b.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00zBACK.nc4','qitot');
% ql_old = ncread('x0009b.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00zBACK.nc4','qltot');
% 
cd /discover/nobackup/drholdaw/tmp.23681/sens.20120702.000000

% u_new = ncread('fsens.eta.20120702_00z_dry.nc4','u');
% v_new = ncread('fsens.eta.20120702_00z_dry.nc4','v');
% t_new = ncread('fsens.eta.20120702_00z_dry.nc4','tv');
% q_new = ncread('fsens.eta.20120702_00z_dry.nc4','sphu');
% p_new = ncread('fsens.eta.20120702_00z_dry.nc4','delp');
% oz_new = ncread('fsens.eta.20120702_00z_dry.nc4','ozone');
% qi_new = ncread('fsens.eta.20120702_00z_dry.nc4','qitot');
% ql_new = ncread('fsens.eta.20120702_00z_dry.nc4','qltot');

u_old = ncread('fsens.eta.20120701_00z_moist.nc4','u');
v_old = ncread('fsens.eta.20120701_00z_moist.nc4','v');
t_old = ncread('fsens.eta.20120701_00z_moist.nc4','tv');
q_old = ncread('fsens.eta.20120701_00z_moist.nc4','sphu');
p_old = ncread('fsens.eta.20120701_00z_moist.nc4','delp');
oz_old = ncread('fsens.eta.20120701_00z_moist.nc4','ozone');
qi_old = ncread('fsens.eta.20120701_00z_moist.nc4','qitot');
ql_old = ncread('fsens.eta.20120701_00z_moist.nc4','qltot');

u_new = ncread('fsens.eta.20120701_00z_moist3.nc4','u');
v_new = ncread('fsens.eta.20120701_00z_moist3.nc4','v');
t_new = ncread('fsens.eta.20120701_00z_moist3.nc4','tv');
q_new = ncread('fsens.eta.20120701_00z_moist3.nc4','sphu');
p_new = ncread('fsens.eta.20120701_00z_moist3.nc4','delp');
oz_new = ncread('fsens.eta.20120701_00z_moist3.nc4','ozone');
qi_new = ncread('fsens.eta.20120701_00z_moist3.nc4','qitot');
ql_new = ncread('fsens.eta.20120701_00z_moist3.nc4','qltot');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/


u_diff = (u_new - u_old);%./u_new;
v_diff = (v_new - v_old);%./v_new;
t_diff = (t_new - t_old);%./t_new;
q_diff = (q_new - q_old);%./q_new;
p_diff = (p_new - p_old);%./p_new;
oz_diff = (oz_new - oz_old);%./oz_new;
qi_diff = (qi_new - qi_old);%./qi_new;
ql_diff = (ql_new - ql_old);%./ql_new;

max(u_diff(:))
max(v_diff(:))
max(t_diff(:))
max(q_diff(:))
max(p_diff(:))
max(oz_diff(:))
max(qi_diff(:))
max(ql_diff(:))

plot_level = 50;

contourf(u_diff(:,:,plot_level)')
colorbar

figure
contourf(v_diff(:,:,plot_level)')
colorbar

figure
contourf(t_diff(:,:,plot_level)')
colorbar

figure
contourf(q_diff(:,:,plot_level)')
colorbar

figure
contourf(p_diff(:,:,plot_level)')
colorbar
