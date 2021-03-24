close all
clear
clc

load coast
coast_lat = lat; clear lat
coast_lon = long; clear long

cd /archive/u/drholdaw/x0009a/prog/Y2012/M06/D30/H21

lon = ncread('x0009a.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','lon');
lat = ncread('x0009a.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','lat');

ua_twe = ncread('x0009a.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','u');
va_twe = ncread('x0009a.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','v');
ta_twe = ncread('x0009a.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','tv');
qa_twe = ncread('x0009a.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','sphu');
pa_twe = ncread('x0009a.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','delp');
o3a_twe = ncread('x0009a.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','ozone');
qia_twe = ncread('x0009a.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','qitot');
qla_twe = ncread('x0009a.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','qltot');

cd /archive/u/drholdaw/x0009b/prog/Y2012/M06/D30/H21

ub_twe = ncread('x0009b.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','u');
vb_twe = ncread('x0009b.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','v');
tb_twe = ncread('x0009b.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','tv');
qb_twe = ncread('x0009b.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','sphu');
pb_twe = ncread('x0009b.fsens_txe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','delp');

cd /archive/u/drholdaw/x0009c/prog/Y2012/M06/D30/H21

uc_twe = ncread('x0009c.fsens_twe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','u');
vc_twe = ncread('x0009c.fsens_twe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','v');
tc_twe = ncread('x0009c.fsens_twe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','tv');
qc_twe = ncread('x0009c.fsens_twe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','sphu');
pc_twe = ncread('x0009c.fsens_twe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','delp');

cd /archive/u/drholdaw/x0009d/prog/Y2012/M06/D30/H21

ud_twe = ncread('x0009d.fsens_twe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','u');
vd_twe = ncread('x0009d.fsens_twe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','v');
td_twe = ncread('x0009d.fsens_twe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','tv');
qd_twe = ncread('x0009d.fsens_twe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','sphu');
pd_twe = ncread('x0009d.fsens_twe.eta.20120630_21z+20120702_00z-20120701_00z.nc4','delp');

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

upb_twe = ub_twe(:,:,plotlevel);
vpb_twe = vb_twe(:,:,plotlevel);
tpb_twe = tb_twe(:,:,plotlevel);
qpb_twe = qb_twe(:,:,plotlevel);
ppb_twe = pb_twe(:,:,plotlevel);

upc_twe = uc_twe(:,:,plotlevel);
vpc_twe = vc_twe(:,:,plotlevel);
tpc_twe = tc_twe(:,:,plotlevel);
qpc_twe = qc_twe(:,:,plotlevel);
ppc_twe = pc_twe(:,:,plotlevel);

upd_twe = ud_twe(:,:,plotlevel);
vpd_twe = vd_twe(:,:,plotlevel);
tpd_twe = td_twe(:,:,plotlevel);
qpd_twe = qd_twe(:,:,plotlevel);
ppd_twe = pd_twe(:,:,plotlevel);


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


figure('position',[707 30 1152 889])
subplot(4,2,1)
[C,h] = contourf(lon,lat,upb_twe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(upb_twe(:))) max(abs(upb_twe(:)))])
title('u')

subplot(4,2,2)
[C,h] = contourf(lon,lat,vpb_twe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(vpb_twe(:))) max(abs(vpb_twe(:)))])
title('v')

subplot(4,2,3)
[C,h] = contourf(lon,lat,tpb_twe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(tpb_twe(:))) max(abs(tpb_twe(:)))])
title('t')

subplot(4,2,4)
[C,h] = contourf(lon,lat,qpb_twe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(qpb_twe(:))) max(abs(qpb_twe(:)))])
title('q')

subplot(4,2,5)
[C,h] = contourf(lon,lat,ppb_twe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(ppb_twe(:))) max(abs(ppb_twe(:)))])
title('p')

asd
figure('position',[1407 47 1152 889])
subplot(4,2,1)
[C,h] = contourf(lon,lat,upc_twe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(upc_twe(:))) max(abs(upc_twe(:)))])
title('u')

subplot(4,2,2)
[C,h] = contourf(lon,lat,vpc_twe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(vpc_twe(:))) max(abs(vpc_twe(:)))])
title('v')

subplot(4,2,3)
[C,h] = contourf(lon,lat,tpc_twe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(tpc_twe(:))) max(abs(tpc_twe(:)))])
title('t')

subplot(4,2,4)
[C,h] = contourf(lon,lat,qpc_twe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(qpc_twe(:))) max(abs(qpc_twe(:)))])
title('q')

subplot(4,2,5)
[C,h] = contourf(lon,lat,ppc_twe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(ppc_twe(:))) max(abs(ppc_twe(:)))])
title('p')


figure('position',[1407 47 1152 889])
subplot(4,2,1)
[C,h] = contourf(lon,lat,upd_twe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(upd_twe(:))) max(abs(upd_twe(:)))])
title('u')

subplot(4,2,2)
[C,h] = contourf(lon,lat,vpd_twe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(vpd_twe(:))) max(abs(vpd_twe(:)))])
title('v')

subplot(4,2,3)
[C,h] = contourf(lon,lat,tpd_twe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(tpd_twe(:))) max(abs(tpd_twe(:)))])
title('t')

subplot(4,2,4)
[C,h] = contourf(lon,lat,qpd_twe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(qpd_twe(:))) max(abs(qpd_twe(:)))])
title('q')

subplot(4,2,5)
[C,h] = contourf(lon,lat,ppd_twe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(ppd_twe(:))) max(abs(ppd_twe(:)))])
title('p')