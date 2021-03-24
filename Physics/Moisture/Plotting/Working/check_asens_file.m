close all
clear
clc

cd /discover/nobackup/drholdaw/iau590b_mp/asens/jong/

% u_twe_moi    = ncread('iau590b.fsens_twe.eta.20120316_21z+20120318_00z-20120317_00z.nc4','u');
% u_txe_moi    = ncread('iau590b_mp.fsens_twe.eta.20120412_21z+20120414_00z-20120413_00z.nc4','u');
% v_twe_moi    = ncread('iau590b.fsens_twe.eta.20120316_21z+20120318_00z-20120317_00z.nc4','v');
% v_txe_moi    = ncread('iau590b_mp.fsens_twe.eta.20120412_21z+20120414_00z-20120413_00z.nc4','v');
% tv_twe_moi   = ncread('iau590b.fsens_twe.eta.20120316_21z+20120318_00z-20120317_00z.nc4','tv');
% tv_txe_moi   = ncread('iau590b_mp.fsens_twe.eta.20120412_21z+20120414_00z-20120413_00z.nc4','tv');
% q_twe_moi    = ncread('iau590b.fsens_twe.eta.20120316_21z+20120318_00z-20120317_00z.nc4','sphu');
% q_txe_moi    = ncread('iau590b_mp.fsens_twe.eta.20120412_21z+20120414_00z-20120413_00z.nc4','sphu');
% delp_twe_moi = ncread('iau590b.fsens_twe.eta.20120316_21z+20120318_00z-20120317_00z.nc4','delp');
% delp_txe_moi = ncread('iau590b_mp.fsens_twe.eta.20120412_21z+20120414_00z-20120413_00z.nc4','delp');

% cd /discover/nobackup/drholdaw/iau590b/asens/fsens_data_dryphys/

% u_twe_dry = ncread('iau590.fsens_twe.eta.20120316_21z+20120318_00z-20120317_00z.nc4','u');
u_txe_dry = ncread('iau590b_mp.fsens_twe.eta.20120412_21z+20120414_00z-20120413_00z.nc4','u');
% v_twe_dry = ncread('iau590.fsens_twe.eta.20120316_21z+20120318_00z-20120317_00z.nc4','v');
v_txe_dry = ncread('iau590b_mp.fsens_twe.eta.20120412_21z+20120414_00z-20120413_00z.nc4','v');
% tv_twe_dry = ncread('iau590.fsens_twe.eta.20120316_21z+20120318_00z-20120317_00z.nc4','tv');
tv_txe_dry = ncread('iau590b_mp.fsens_twe.eta.20120412_21z+20120414_00z-20120413_00z.nc4','tv');
% q_twe_dry = ncread('iau590.fsens_twe.eta.20120316_21z+20120318_00z-20120317_00z.nc4','sphu');
q_txe_dry = ncread('iau590b_mp.fsens_twe.eta.20120412_21z+20120414_00z-20120413_00z.nc4','sphu');
% delp_twe_dry = ncread('iau590.fsens_twe.eta.20120316_21z+20120318_00z-20120317_00z.nc4','delp');
delp_txe_dry = ncread('iau590b_mp.fsens_twe.eta.20120412_21z+20120414_00z-20120413_00z.nc4','delp');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/




plot_level = 3;


figure
subplot(2,3,1)
contourf(u_txe_dry(:,:,plot_level)')
colorbar

subplot(2,3,2)
contourf(v_txe_dry(:,:,plot_level)')
colorbar

subplot(2,3,3)
contourf(tv_txe_dry(:,:,plot_level)')
colorbar

subplot(2,3,4)
contourf(q_txe_dry(:,:,plot_level)')
colorbar

subplot(2,3,5)
contourf(delp_txe_dry(:,:,plot_level)')
colorbar
% 
% 
% 
% 
% figure
% subplot(2,3,1)
% contourf(u_twe_dry(:,:,plot_level)')
% colorbar
% 
% subplot(2,3,2)
% contourf(v_twe_dry(:,:,plot_level)')
% colorbar
% 
% subplot(2,3,3)
% contourf(tv_twe_dry(:,:,plot_level)')
% colorbar
% 
% subplot(2,3,4)
% contourf(q_twe_dry(:,:,plot_level)')
% colorbar
% 
% subplot(2,3,5)
% contourf(delp_twe_dry(:,:,plot_level)')
% colorbar






% figure
% subplot(2,3,1)
% contourf(u_txe_moi(:,:,plot_level)')
% colorbar
% 
% subplot(2,3,2)
% contourf(v_txe_moi(:,:,plot_level)')
% colorbar
% 
% subplot(2,3,3)
% contourf(tv_txe_moi(:,:,plot_level)')
% colorbar
% 
% subplot(2,3,4)
% contourf(q_txe_moi(:,:,plot_level)')
% colorbar
% 
% subplot(2,3,5)
% contourf(delp_txe_moi(:,:,plot_level)')
% colorbar
% 
% 
% 
% 
% figure
% subplot(2,3,1)
% contourf(u_twe_moi(:,:,plot_level)')
% colorbar
% 
% subplot(2,3,2)
% contourf(v_twe_moi(:,:,plot_level)')
% colorbar
% 
% subplot(2,3,3)
% contourf(tv_twe_moi(:,:,plot_level)')
% colorbar
% 
% subplot(2,3,4)
% contourf(q_twe_moi(:,:,plot_level)')
% colorbar
% 
% subplot(2,3,5)
% contourf(delp_twe_moi(:,:,plot_level)')
% colorbar