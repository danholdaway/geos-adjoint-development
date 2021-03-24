close all
clear
clc


cd /discover/nobackup/drholdaw/atmp.22292/;
file = 'x0011dh_a.prog.eta.20130116_0100z.nc4';

lon = ncread(file,'lon');
lat = ncread(file,'lat');
lev = ncread(file,'lev');

qi_free = ncread(file,'qitot');
ql_free = ncread(file,'qltot');
qils_free = ncread(file,'QILS');
qlls_free = ncread(file,'QLLS');
qicn_free = ncread(file,'QICN');
qlcn_free = ncread(file,'QLCN');


cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

sumql = zeros(1,72);
sumqi = zeros(1,72);
sumqlls = zeros(1,72);
sumqils = zeros(1,72);
sumqlcn = zeros(1,72);
sumqicn = zeros(1,72);

for k = 1:72
    sumql(k) = sum(sum(ql_free(:,:,k)));
    sumqi(k) = sum(sum(qi_free(:,:,k)));
    sumqlls(k) = sum(sum(qlls_free(:,:,k)));
    sumqlcn(k) = sum(sum(qlcn_free(:,:,k)));
    sumqils(k) = sum(sum(qils_free(:,:,k)));
    sumqicn(k) = sum(sum(qicn_free(:,:,k)));
end


plot(sumql,1:72)
hold on
plot(sumqlls,1:72,'r')
plot(sumqlcn,1:72,'g')
set(gca,'YDir','reverse')
ylim([30 72])

figure
plot(sumqi,1:72)
hold on
plot(sumqils,1:72,'r')
plot(sumqicn,1:72,'g')
set(gca,'YDir','reverse')
ylim([30 72])


figure
set(gcf,'position',[719    30   560   889])
subplot(3,1,1)
contour(lon,lat,ql_free(:,:,63)')
colorbar

subplot(3,1,2)
contour(lon,lat,qlls_free(:,:,63)')
colorbar

subplot(3,1,3)
contour(lon,lat,qlcn_free(:,:,63)')
colorbar


figure
set(gcf,'position',[719    30   560   889])
subplot(3,1,1)
contour(lon,lat,qi_free(:,:,45)')
colorbar

subplot(3,1,2)
contour(lon,lat,qils_free(:,:,45)')
colorbar

subplot(3,1,3)
contour(lon,lat,qicn_free(:,:,45)')
colorbar