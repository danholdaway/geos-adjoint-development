close all
clear
clc

cd /discover/nobackup/drholdaw/tmp.31275/sens.20120318.000000

time = ncread('JgradfX.eta.nc4','time');
lon = ncread('JgradfX.eta.nc4','lon');
lat = ncread('JgradfX.eta.nc4','lat');
lev = ncread('JgradfX.eta.nc4','lev');

% U_drynorm = ncread('Jgradf_txe.eta.nc4','u');
% V_drynorm = ncread('Jgradf_txe.eta.nc4','v');
% T_drynorm = ncread('Jgradf_txe.eta.nc4','tv');
% Q_drynorm = ncread('Jgradf_txe.eta.nc4','sphu');
% P_drynorm = ncread('Jgradf_txe.eta.nc4','delp');
% 
% cd /discover/nobackup/drholdaw/tmp.21867/sens.20120318.000000/
% 
% U_wetnorm = ncread('Jgradf_twe.eta.nc4','u');
% V_wetnorm = ncread('Jgradf_twe.eta.nc4','v');
% T_wetnorm = ncread('Jgradf_twe.eta.nc4','tv');
% Q_wetnorm = ncread('Jgradf_twe.eta.nc4','sphu');
% P_wetnorm = ncread('Jgradf_twe.eta.nc4','delp');

U_X = ncread('JgradfX.eta.nc4','U');
V_X = ncread('JgradfX.eta.nc4','V');
T_X = ncread('JgradfX.eta.nc4','TV');
Q_X = ncread('JgradfX.eta.nc4','QV');
P_X = ncread('JgradfX.eta.nc4','DP');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/


plotlevel = 65;

% subplot(2,3,1)
% contour(lon,lat,(U_drynorm(:,:,plotlevel))')
% title('U \prime')
% colorbar
% 
% subplot(2,3,2)
% contour(lon,lat,V_drynorm(:,:,plotlevel)')
% title('V \prime')
% colorbar
% 
% subplot(2,3,3)
% contour(lon,lat,T_drynorm(:,:,plotlevel)')
% title('T \prime')
% colorbar
% 
% subplot(2,3,4)
% contour(lon,lat,Q_drynorm(:,:,plotlevel)')
% title('Q \prime')
% colorbar
% 
% subplot(2,3,5)
% contour(lon,lat,P_drynorm(:,:,plotlevel)')
% title('P \prime')
% colorbar
% 
% figure
% 
% subplot(2,3,1)
% contour(lon,lat,(U_wetnorm(:,:,plotlevel))')
% title('U \prime')
% colorbar
% 
% subplot(2,3,2)
% contour(lon,lat,V_wetnorm(:,:,plotlevel)')
% title('V \prime')
% colorbar
% 
% subplot(2,3,3)
% contour(lon,lat,T_wetnorm(:,:,plotlevel)')
% title('T \prime')
% colorbar
% 
% subplot(2,3,4)
% contour(lon,lat,Q_wetnorm(:,:,plotlevel)')
% title('Q \prime')
% colorbar
% 
% subplot(2,3,5)
% contour(lon,lat,P_wetnorm(:,:,plotlevel)')
% title('P \prime')
% colorbar

figure

subplot(2,3,1)
contour(lon,lat,(U_X(:,:,plotlevel))')
title('U \prime')
colorbar

subplot(2,3,2)
contour(lon,lat,V_X(:,:,plotlevel)')
title('V \prime')
colorbar

subplot(2,3,3)
contour(lon,lat,T_X(:,:,plotlevel)')
title('T \prime')
colorbar

subplot(2,3,4)
contour(lon,lat,Q_X(:,:,plotlevel)')
title('Q \prime')
colorbar

subplot(2,3,5)
contour(lon,lat,P_X(:,:,plotlevel)')
title('P \prime')
colorbar

