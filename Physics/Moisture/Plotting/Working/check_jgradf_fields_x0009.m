close all
clear
clc

load coast
coast_lat = lat; clear lat
coast_lon = long; clear long

cd /archive/u/drholdaw/x0009a/prog/Y2012/M06/D30/H21

lon = ncread('x0009a.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','lon');
lat = ncread('x0009a.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','lat');

ua_twe = ncread('x0009a.Jgradf_twe.eta.20120630_21z+20120702_00z.nc4','u');
va_twe = ncread('x0009a.Jgradf_twe.eta.20120630_21z+20120702_00z.nc4','v');
ta_twe = ncread('x0009a.Jgradf_twe.eta.20120630_21z+20120702_00z.nc4','tv');
qa_twe = ncread('x0009a.Jgradf_twe.eta.20120630_21z+20120702_00z.nc4','sphu');
pa_twe = ncread('x0009a.Jgradf_twe.eta.20120630_21z+20120702_00z.nc4','delp');
o3a_twe = ncread('x0009a.Jgradf_twe.eta.20120630_21z+20120702_00z.nc4','ozone');
qia_twe = ncread('x0009a.Jgradf_twe.eta.20120630_21z+20120702_00z.nc4','qitot');
qla_twe = ncread('x0009a.Jgradf_twe.eta.20120630_21z+20120702_00z.nc4','qltot');


cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/

plotlevel = 63;
grey = 0.5;

upa_twe = ua_twe(:,:,plotlevel);
vpa_twe = va_twe(:,:,plotlevel);
tpa_twe = ta_twe(:,:,plotlevel);
qpa_twe = qa_twe(:,:,plotlevel);
ppa_twe = pa_twe(:,:,plotlevel);
o3pa_twe = o3a_twe(:,:,plotlevel);
qipa_twe = qia_twe(:,:,plotlevel);
qlpa_twe = qla_twe(:,:,plotlevel);


figure('position',[4 30 1152 889])
subplot(4,2,1)
[C,h] = contourf(lon,lat,upa_twe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(upa_twe(:))) max(abs(upa_twe(:)))])
title('u')

subplot(4,2,2)
[C,h] = contourf(lon,lat,vpa_twe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(vpa_twe(:))) max(abs(vpa_twe(:)))])
title('v')

subplot(4,2,3)
[C,h] = contourf(lon,lat,tpa_twe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('t')
caxis([-max(abs(tpa_twe(:))) max(abs(tpa_twe(:)))])

subplot(4,2,4)
[C,h] = contourf(lon,lat,qpa_twe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(qpa_twe(:))) max(abs(qpa_twe(:)))])
title('q')

subplot(4,2,5)
[C,h] = contourf(lon,lat,ppa_twe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(ppa_twe(:))) max(abs(ppa_twe(:)))])
title('p')

subplot(4,2,6)
[C,h] = contourf(lon,lat,o3pa_twe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(ppa_twe(:))) max(abs(ppa_twe(:)))])
title('o_3')

subplot(4,2,7)
[C,h] = contourf(lon,lat,qipa_twe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(ppa_twe(:))) max(abs(ppa_twe(:)))])
title('q_i')

subplot(4,2,8)
[C,h] = contourf(lon,lat,qlpa_twe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(ppa_twe(:))) max(abs(ppa_twe(:)))])
title('q_l')
