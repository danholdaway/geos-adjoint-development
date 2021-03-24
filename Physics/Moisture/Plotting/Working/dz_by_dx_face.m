close all
clear
clc

load coast
coast_lat = lat; clear lat
coast_lon = long; clear long
grey = 0.0;

cd /discover/nobackup/drholdaw/tmp.4975

phis1 = ncread('iau590a.traj.lcv.20120317_0000z.nc4','phis');

lons1 = ncread('iau590a.traj.M2.lcv.20120317_0000z.nc4','LONS_cube');
lats1 = ncread('iau590a.traj.M2.lcv.20120317_0000z.nc4','LATS_cube');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/

lons1 = radtodeg(lons1);
lats1 = radtodeg(lats1);

phis1 = phis1/9.81;

PHIS = zeros(6,90,90);
LATS = zeros(6,90,90);
LONS = zeros(6,90,90);

Dz = zeros(6,90,90);
Dx = zeros(6,90,90);

PHIS(1,:,:) = phis1(:,1:90);
PHIS(2,:,:) = phis1(:,91:180);
PHIS(3,:,:) = phis1(:,181:270);
PHIS(4,:,:) = phis1(:,271:360);
PHIS(5,:,:) = phis1(:,361:450);
PHIS(6,:,:) = phis1(:,451:540);

LATS(1,:,:) = lats1(:,1:90);
LATS(2,:,:) = lats1(:,91:180);
LATS(3,:,:) = lats1(:,181:270);
LATS(4,:,:) = lats1(:,271:360);
LATS(5,:,:) = lats1(:,361:450);
LATS(6,:,:) = lats1(:,451:540);

LONS(1,:,:) = lons1(:,1:90);
LONS(2,:,:) = lons1(:,91:180);
LONS(3,:,:) = lons1(:,181:270);
LONS(4,:,:) = lons1(:,271:360);
LONS(5,:,:) = lons1(:,361:450);
LONS(6,:,:) = lons1(:,451:540);

LATSr = deg2rad(LATS);
LONSr = deg2rad(LONS);

[f,e,n] = size(LATS);

LAND = zeros(f,e,n);
for k = 1:f
    for i = 1:e
        for j = 1:n
            if PHIS(k,i,j) <= 10
                LAND(k,i,j) = 1;
            end
        end
    end
end

for k = 1:f
    for i = 2:e-1
        for j = 2:n-1

            dz_i(1) = (abs(PHIS(k,i,j) - PHIS(k,i+1,j  ))) ;
            dz_i(2) = (abs(PHIS(k,i,j) - PHIS(k,i-1,j  ))) ;
            dz_i(3) = (abs(PHIS(k,i,j) - PHIS(k,i  ,j+1))) ;
            dz_i(4) = (abs(PHIS(k,i,j) - PHIS(k,i  ,j-1))) ;
            Dz(k,i,j) = max(dz_i);
            
            d1 = sin(0.5*(LATSr(k,i,j)-LATSr(k,i+1,j)));
            d2 = sin(0.5*(LONSr(k,i,j)-LONSr(k,i+1,j)));
            d3 = cos(LATSr(k,i,j))*cos(LATSr(k,i+1,j));
            d = (2*earthRadius)*asin(sqrt(d1*d1+d2*d2*d3));
            du_i(1) = abs(d);

            d1 = sin(0.5*(LATSr(k,i,j)-LATSr(k,i-1,j)));
            d2 = sin(0.5*(LONSr(k,i,j)-LONSr(k,i-1,j)));
            d3 = cos(LATSr(k,i,j))*cos(LATSr(k,i-1,j));
            d = (2*earthRadius)*asin(sqrt(d1*d1+d2*d2*d3));
            du_i(2) = abs(d);

            d1 = sin(0.5*(LATSr(k,i,j)-LATSr(k,i,j+1)));
            d2 = sin(0.5*(LONSr(k,i,j)-LONSr(k,i,j+1)));
            d3 = cos(LATSr(k,i,j))*cos(LATSr(k,i,j+1));
            d = (2*earthRadius)*asin(sqrt(d1*d1+d2*d2*d3));
            du_i(3) = abs(d);

            d1 = sin(0.5*(LATSr(k,i,j)-LATSr(k,i,j-1)));
            d2 = sin(0.5*(LONSr(k,i,j)-LONSr(k,i,j-1)));
            d3 = cos(LATSr(k,i,j))*cos(LATSr(k,i,j-1));
            d = (2*earthRadius)*asin(sqrt(d1*d1+d2*d2*d3));
            du_i(4) = abs(d);

            Dx(k,i,j) = max(du_i);
            
        end
    end
end

Dz1(:,:) = Dz(1,:,:);
Dz2(:,:) = Dz(2,:,:);
Dz3(:,:) = Dz(3,:,:);
Dz4(:,:) = Dz(4,:,:);
Dz5(:,:) = Dz(5,:,:);
Dz6(:,:) = Dz(6,:,:);

Dz_Npole = zeros(45,360);
Dz_Spole = zeros(45,360);

Dz_Npole(:,1:90) = rot90(Dz3(:,46:end),3);
Dz_Npole(:,91:180) = rot90(Dz3(1:45,:),2);
Dz_Npole(:,181:270) = rot90(Dz3(:,1:45));
Dz_Npole(:,271:end) = Dz3(46:end,:);
Dz_Spole(:,1:90) = Dz6(1:45,:);
Dz_Spole(:,91:180) = rot90(Dz6(:,46:end),1);
Dz_Spole(:,181:270) = rot90(Dz6(46:end,:),2);
Dz_Spole(:,271:end) = rot90(Dz6(:,1:45),3);

dz = flipud([Dz_Npole; [ Dz5(:,:) rot90(Dz1(:,:)) rot90(Dz2(:,:)) Dz4(:,:)]; Dz_Spole]);

Dx1(:,:) = Dx(1,:,:);
Dx2(:,:) = Dx(2,:,:);
Dx3(:,:) = Dx(3,:,:);
Dx4(:,:) = Dx(4,:,:);
Dx5(:,:) = Dx(5,:,:);
Dx6(:,:) = Dx(6,:,:);

Dx_Npole = zeros(45,360);
Dx_Spole = zeros(45,360);

Dx_Npole(:,1:90) = rot90(Dx3(:,46:end),3);
Dx_Npole(:,91:180) = rot90(Dx3(1:45,:),2);
Dx_Npole(:,181:270) = rot90(Dx3(:,1:45));
Dx_Npole(:,271:end) = Dx3(46:end,:);
Dx_Spole(:,1:90) = Dx6(1:45,:);
Dx_Spole(:,91:180) = rot90(Dx6(:,46:end),1);
Dx_Spole(:,181:270) = rot90(Dx6(46:end,:),2);
Dx_Spole(:,271:end) = rot90(Dx6(:,1:45),3);

dx = flipud([Dx_Npole; [ Dx5(:,:) rot90(Dx1(:,:)) rot90(Dx2(:,:)) Dx4(:,:)]; Dx_Spole]);

LAND1(:,:) = LAND(1,:,:);
LAND2(:,:) = LAND(2,:,:);
LAND3(:,:) = LAND(3,:,:);
LAND4(:,:) = LAND(4,:,:);
LAND5(:,:) = LAND(5,:,:);
LAND6(:,:) = LAND(6,:,:);

LAND_Npole = zeros(45,360);
LAND_Spole = zeros(45,360);

LAND_Npole(:,1:90) = rot90(LAND3(:,46:end),3);
LAND_Npole(:,91:180) = rot90(LAND3(1:45,:),2);
LAND_Npole(:,181:270) = rot90(LAND3(:,1:45));
LAND_Npole(:,271:end) = LAND3(46:end,:);
LAND_Spole(:,1:90) = LAND6(1:45,:);
LAND_Spole(:,91:180) = rot90(LAND6(:,46:end),1);
LAND_Spole(:,181:270) = rot90(LAND6(46:end,:),2);
LAND_Spole(:,271:end) = rot90(LAND6(:,1:45),3);

land = flipud([LAND_Npole; [ Dx5(:,:) rot90(LAND1(:,:)) rot90(LAND2(:,:)) LAND4(:,:)]; LAND_Spole]);

LATS1(:,:) = LATS(1,:,:);
LATS2(:,:) = LATS(2,:,:);
LATS3(:,:) = LATS(3,:,:);
LATS4(:,:) = LATS(4,:,:);
LATS5(:,:) = LATS(5,:,:);
LATS6(:,:) = LATS(6,:,:);

LATS_Npole = zeros(45,360);
LATS_Spole = zeros(45,360);

LATS_Npole(:,1:90) = rot90(LATS3(:,46:end),3);
LATS_Npole(:,91:180) = rot90(LATS3(1:45,:),2);
LATS_Npole(:,181:270) = rot90(LATS3(:,1:45));
LATS_Npole(:,271:end) = LATS3(46:end,:);
LATS_Spole(:,1:90) = LATS6(1:45,:);
LATS_Spole(:,91:180) = rot90(LATS6(:,46:end),1);
LATS_Spole(:,181:270) = rot90(LATS6(46:end,:),2);
LATS_Spole(:,271:end) = rot90(LATS6(:,1:45),3);

lats = flipud([LATS_Npole; [ LATS5(:,:) rot90(LATS1(:,:)) rot90(LATS2(:,:)) LATS4(:,:)]; LATS_Spole]);

LONS1(:,:) = LONS(1,:,:);
LONS2(:,:) = LONS(2,:,:);
LONS3(:,:) = LONS(3,:,:);
LONS4(:,:) = LONS(4,:,:);
LONS5(:,:) = LONS(5,:,:);
LONS6(:,:) = LONS(6,:,:);

LONS_Npole = zeros(45,360);
LONS_Spole = zeros(45,360);

LONS_Npole(:,1:90) = rot90(LONS3(:,46:end),3);
LONS_Npole(:,91:180) = rot90(LONS3(1:45,:),2);
LONS_Npole(:,181:270) = rot90(LONS3(:,1:45));
LONS_Npole(:,271:end) = LONS3(46:end,:);
LONS_Spole(:,1:90) = LONS6(1:45,:);
LONS_Spole(:,91:180) = rot90(LONS6(:,46:end),1);
LONS_Spole(:,181:270) = rot90(LONS6(46:end,:),2);
LONS_Spole(:,271:end) = rot90(LONS6(:,1:45),3);

lons = flipud([LONS_Npole; [ LONS5(:,:) rot90(LONS1(:,:)) rot90(LONS2(:,:)) LONS4(:,:)]; LONS_Spole]);


figure('position',[151 236 1128 683])
contourf(lons)
hold on
contour(land,[0 0],'Color',[1 1 1])
colorbar

figure('position',[151 236 1128 683])
contourf(lats)
hold on
contour(land,[0 0],'Color',[1 1 1])
colorbar

figure('position',[151 236 1128 683])
contourf(dx)
hold on
contour(land,[0 0],'Color',[1 1 1])
colorbar

figure('position',[151 236 1128 683])
contourf(dz)
hold on
contour(land,[0 0],'Color',[1 1 1])
colorbar

figure('position',[151 236 1128 683])
contourf(dz./dx)
hold on
contour(land,[0 0],'Color',[1 1 1])
colorbar