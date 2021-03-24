close all
clear
clc

load coast
coast_lat = lat; clear lat
coast_lon = long; clear long
grey = 0.0;

cd /discover/nobackup/drholdaw/tmp.iau590/tmp.31275

phis1 = ncread('iau590.prog.sfc.20120317_00z.nc4','PHIS');

lons1 = ncread('iau590.prog.sfc.20120317_00z.nc4','lon');
lats1 = ncread('iau590.prog.sfc.20120317_00z.nc4','lat');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/

lons1 = radtodeg(lons1);
lats1 = radtodeg(lats1);

phis1 = phis1/9.81;

PHIS1 = phis1(:,1:90);
PHIS2 = phis1(:,91:180);
PHIS3 = phis1(:,181:270);
PHIS4 = phis1(:,271:360);
PHIS5 = phis1(:,361:450);
PHIS6 = phis1(:,451:540);

DPHIS1 = zeros(90,90);
DPHIS2 = zeros(90,90);
DPHIS3 = zeros(90,90);
DPHIS4 = zeros(90,90);
DPHIS5 = zeros(90,90);
DPHIS6 = zeros(90,90);
for i= 2:89
    for j = 2:89

        dzdu_i_gc(1) = (abs(PHIS1(i,j) - PHIS1(i+1,j  ))) ;
        dzdu_i_gc(2) = (abs(PHIS1(i,j) - PHIS1(i-1,j  ))) ;
        dzdu_i_gc(3) = (abs(PHIS1(i,j) - PHIS1(i  ,j+1))) ;
        dzdu_i_gc(4) = (abs(PHIS1(i,j) - PHIS1(i  ,j-1))) ;
        DPHIS1(i,j) = max(dzdu_i_gc);
        dzdu_i_gc(1) = (abs(PHIS2(i,j) - PHIS2(i+1,j  ))) ;
        dzdu_i_gc(2) = (abs(PHIS2(i,j) - PHIS2(i-1,j  ))) ;
        dzdu_i_gc(3) = (abs(PHIS2(i,j) - PHIS2(i  ,j+1))) ;
        dzdu_i_gc(4) = (abs(PHIS2(i,j) - PHIS2(i  ,j-1))) ;
        DPHIS2(i,j) = max(dzdu_i_gc);
        dzdu_i_gc(1) = (abs(PHIS3(i,j) - PHIS3(i+1,j  ))) ;
        dzdu_i_gc(2) = (abs(PHIS3(i,j) - PHIS3(i-1,j  ))) ;
        dzdu_i_gc(3) = (abs(PHIS3(i,j) - PHIS3(i  ,j+1))) ;
        dzdu_i_gc(4) = (abs(PHIS3(i,j) - PHIS3(i  ,j-1))) ;
        DPHIS3(i,j) = max(dzdu_i_gc);
        dzdu_i_gc(1) = (abs(PHIS4(i,j) - PHIS4(i+1,j  ))) ;
        dzdu_i_gc(2) = (abs(PHIS4(i,j) - PHIS4(i-1,j  ))) ;
        dzdu_i_gc(3) = (abs(PHIS4(i,j) - PHIS4(i  ,j+1))) ;
        dzdu_i_gc(4) = (abs(PHIS4(i,j) - PHIS4(i  ,j-1))) ;
        DPHIS4(i,j) = max(dzdu_i_gc);
        dzdu_i_gc(1) = (abs(PHIS5(i,j) - PHIS5(i+1,j  ))) ;
        dzdu_i_gc(2) = (abs(PHIS5(i,j) - PHIS5(i-1,j  ))) ;
        dzdu_i_gc(3) = (abs(PHIS5(i,j) - PHIS5(i  ,j+1))) ;
        dzdu_i_gc(4) = (abs(PHIS5(i,j) - PHIS5(i  ,j-1))) ;
        DPHIS5(i,j) = max(dzdu_i_gc);
        dzdu_i_gc(1) = (abs(PHIS6(i,j) - PHIS6(i+1,j  ))) ;
        dzdu_i_gc(2) = (abs(PHIS6(i,j) - PHIS6(i-1,j  ))) ;
        dzdu_i_gc(3) = (abs(PHIS6(i,j) - PHIS6(i  ,j+1))) ;
        dzdu_i_gc(4) = (abs(PHIS6(i,j) - PHIS6(i  ,j-1))) ;
        DPHIS6(i,j) = max(dzdu_i_gc);
        
    end
end


PHIS_Npole = zeros(45,360);
PHIS_Spole = zeros(45,360);
PHIS_Npole(:,1:90) = rot90(PHIS3(:,46:end),3);
PHIS_Npole(:,91:180) = rot90(PHIS3(1:45,:),2);
PHIS_Npole(:,181:270) = rot90(PHIS3(:,1:45));
PHIS_Npole(:,271:end) = PHIS3(46:end,:);

% PHIS_Npole(:,1:90) = rot90(PHIS3(:,1:45),3);
% PHIS_Npole(:,91:180) = PHIS3(1:45,:);
% PHIS_Npole(:,181:270) = rot90(PHIS3(:,46:end),1);
% PHIS_Npole(:,271:end) = rot90(PHIS3(46:end,:),2);

PHIS_Spole(:,1:90) = PHIS6(1:45,:);
PHIS_Spole(:,91:180) = rot90(PHIS6(:,46:end),1);
PHIS_Spole(:,181:270) = rot90(PHIS6(46:end,:),2);
PHIS_Spole(:,271:end) = rot90(PHIS6(:,1:45),3);

phis = flipud([PHIS_Npole; [ PHIS5(:,:) rot90(PHIS1(:,:)) rot90(PHIS2(:,:)) PHIS4(:,:)]; PHIS_Spole]);
% phis = [phis(:,end-36:end) phis(:,1:end-37)];

DPHIS_Npole = zeros(45,360);
DPHIS_Spole = zeros(45,360);
DPHIS_Npole(:,1:90) = rot90(DPHIS3(:,46:end),3);
DPHIS_Npole(:,91:180) = rot90(DPHIS3(1:45,:),2);
DPHIS_Npole(:,181:270) = rot90(DPHIS3(:,1:45));
DPHIS_Npole(:,271:end) = DPHIS3(46:end,:);
DPHIS_Spole(:,1:90) = DPHIS6(1:45,:);
DPHIS_Spole(:,91:180) = rot90(DPHIS6(:,46:end),1);
DPHIS_Spole(:,181:270) = rot90(DPHIS6(46:end,:),2);
DPHIS_Spole(:,271:end) = rot90(DPHIS6(:,1:45),3);

dphis = flipud([DPHIS_Npole; [ DPHIS5(:,:) rot90(DPHIS1(:,:)) rot90(DPHIS2(:,:)) DPHIS4(:,:)]; DPHIS_Spole]);
% dphis = [dphis(:,end-36:end) dphis(:,1:end-37)];

LONS1 = lons1(:,1:90);
LONS2 = lons1(:,91:180);
LONS3 = lons1(:,181:270);
LONS4 = lons1(:,271:360);
LONS5 = lons1(:,361:450);
LONS6 = lons1(:,451:540);

LONS3 = fliplr(LONS3);

LONS_Npole = zeros(45,360);
LONS_Spole = zeros(45,360);
LONS_Npole(:,1:90) = rot90(LONS3(:,46:end),3);
LONS_Npole(:,91:180) = rot90(LONS3(1:45,:),2);
LONS_Npole(:,181:270) = rot90(LONS3(:,1:45));
LONS_Npole(:,271:end) = LONS3(46:end,:);

LONS_Npole(:,1:90) = rot90(LONS3(:,1:45),3);
LONS_Npole(:,91:180) = LONS3(1:45,:);
LONS_Npole(:,181:270) = rot90(LONS3(:,46:end),1);
LONS_Npole(:,271:end) = rot90(LONS3(46:end,:),2);

LONS_Spole(:,1:90) = LONS6(1:45,:);
LONS_Spole(:,91:180) = rot90(LONS6(:,46:end),1);
LONS_Spole(:,181:270) = rot90(LONS6(46:end,:),2);
LONS_Spole(:,271:end) = rot90(LONS6(:,1:45),3);

LONS_FACE = flipud([ LONS5(:,:) rot90(LONS1(:,:)) rot90(LONS2(:,:)) LONS4(:,:)]);
LONS_Npole = flipud(LONS_Npole);
LONS_Spole = flipud(LONS_Spole);

LONS_Npole1 = flipud(LONS_Npole);
LONS_Spole1 = flipud(LONS_Spole);

for i = 1:45
    
    LONS_Npole(i,:) = LONS_FACE(end-i+1,:);
    LONS_Spole(i,:) = LONS_FACE(1+1-1,:);
    
end

lons = [LONS_Spole; LONS_FACE ; LONS_Npole];
% lons = [lons(:,end-36:end) lons(:,1:end-37)];

lons1 = [LONS_Spole1; LONS_FACE ; LONS_Npole1];
% lons1 = [lons1(:,end-36:end) lons1(:,1:end-37)];

% lons(:,1:181) = lons(:,1:181)-360;
% lons1(:,1:181) = lons1(:,1:181)-360;

LATS1 = lats1(:,1:90);
LATS2 = lats1(:,91:180);
LATS3 = lats1(:,181:270);
LATS4 = lats1(:,271:360);
LATS5 = lats1(:,361:450);
LATS6 = lats1(:,451:540);

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

LATS_FACE = flipud([ LATS5(:,:) rot90(LATS1(:,:)) rot90(LATS2(:,:)) LATS4(:,:)]);
LATS_Npole = flipud(LATS_Npole);
LATS_Spole = flipud(LATS_Spole);

lats = [LATS_Spole; LATS_FACE ; LATS_Npole];
% lats = [lats(:,end-36:end) lats(:,1:end-37)];


lon_plot = -180:179;
lat_plot = -90:89;


[n,e] = size(lats);

land = zeros(n,e);
for i = 1:n
    for j = 1:e
        if phis(i,j) <= 10
            land(i,j) = 1;
        end
    end
end

figure('position',[151 236 1128 683])
contourf(phis)
hold on
contour(land,[0 0],'Color',[1 1 1])
colorbar
title('u')

figure('position',[151 236 1128 683])
contourf(lats)
hold on
contour(land,[0 0],'Color',[1 1 1])
colorbar
title('lats')

figure('position',[151 236 1128 683])
contourf(lons1,'ShowText','on')
hold on
contour(land,[0 0],'Color',[1 1 1])
colorbar
title('lats')

figure('position',[151 236 1128 683])
contourf(lons,'ShowText','on')
hold on
contour(land,[0 0],'Color',[1 1 1])
colorbar
title('lons')


latsr = deg2rad(lats);
lonsr = deg2rad(lons);

dz = zeros(n,e);
du = zeros(n,e);
dzdu_1 = zeros(n,e);
dzdu_2 = zeros(n,e);


for i = 2:n-1
    for j = 2:e-1

        dz_i(1) = (abs(phis(i,j) - phis(i+1,j  ))) ;
        dz_i(2) = (abs(phis(i,j) - phis(i-1,j  ))) ;
        dz_i(3) = (abs(phis(i,j) - phis(i  ,j+1))) ;
        dz_i(4) = (abs(phis(i,j) - phis(i  ,j-1))) ;
        
        dz(i,j) = max(dz_i);

        d1 = sin(0.5*(latsr(i,j)-latsr(i+1,j)));
        d2 = sin(0.5*(lonsr(i,j)-lonsr(i+1,j)));
        d3 = cos(latsr(i,j))*cos(latsr(i+1,j));
        d = (2*earthRadius)*asin(sqrt(d1*d1+d2*d2*d3));
        du_i(1) = abs(d);
        
        d1 = sin(0.5*(latsr(i,j)-latsr(i-1,j)));
        d2 = sin(0.5*(lonsr(i,j)-lonsr(i-1,j)));
        d3 = cos(latsr(i,j))*cos(latsr(i-1,j));
        d = (2*earthRadius)*asin(sqrt(d1*d1+d2*d2*d3));
        du_i(2) = abs(d);
        
        d1 = sin(0.5*(latsr(i,j)-latsr(i,j+1)));
        d2 = sin(0.5*(lonsr(i,j)-lonsr(i,j+1)));
        d3 = cos(latsr(i,j))*cos(latsr(i,j+1));
        d = (2*earthRadius)*asin(sqrt(d1*d1+d2*d2*d3));
        du_i(3) = abs(d);
        
        d1 = sin(0.5*(latsr(i,j)-latsr(i,j-1)));
        d2 = sin(0.5*(lonsr(i,j)-lonsr(i,j-1)));
        d3 = cos(latsr(i,j))*cos(latsr(i,j-1));
        d = (2*earthRadius)*asin(sqrt(d1*d1+d2*d2*d3));
        du_i(4) = abs(d);
        
        du(i,j) = max(du_i);
        
        dzdu_i_1(1) = dphis(i,j) / du_i(1);
        dzdu_i_1(2) = dphis(i,j) / du_i(2);
        dzdu_i_1(3) = dphis(i,j) / du_i(3);
        dzdu_i_1(4) = dphis(i,j) / du_i(4);               
        dzdu_1(i,j) = max(dzdu_i_1);
               
        dzdu_i_2(1) = (abs(phis(i,j) - phis(i+1,j  ))) / du_i(1);
        dzdu_i_2(2) = (abs(phis(i,j) - phis(i-1,j  ))) / du_i(2);
        dzdu_i_2(3) = (abs(phis(i,j) - phis(i  ,j+1))) / du_i(3);
        dzdu_i_2(4) = (abs(phis(i,j) - phis(i  ,j-1))) / du_i(4); 
        dzdu_2(i,j) = max(dzdu_i_2);
        
        
    end
end

lat_min_i = 50;
lat_max_i = 160;

figure('position',[151 236 1128 683])
contourf(du(2:end-1,2:end-1))
hold on
contour(land(2:end-1,2:end-1),[0 0],'Color',[1 1 1])
colorbar

figure('position',[151 236 1128 683])
contourf(dz(2:end-1,2:end-1))
hold on
contour(land(2:end-1,2:end-1),[0 0],'Color',[1 1 1])
colorbar

figure('position',[151 236 1128 683])
contourf(dphis(2:end-1,2:end-1))
hold on
contour(land(2:end-1,2:end-1),[0 0],'Color',[1 1 1])
colorbar

figure('position',[151 236 1128 683])
contourf(dzdu_1(lat_min_i:lat_max_i,2:end-1))
hold on
contour(land(lat_min_i:lat_max_i,2:end-1),[0 0],'Color',[1 1 1])
colorbar

figure('position',[151 236 1128 683])
contourf(dzdu_2(lat_min_i:lat_max_i,2:end-1))
hold on
contour(land(lat_min_i:lat_max_i,2:end-1),[0 0],'Color',[1 1 1])
colorbar


A = dzdu_1(lat_min_i:lat_max_i,2:end-1);
[num idx] = max(A(:))
[x y] = find(dzdu_1 == num)
lats(x,y)
lons(x,y)
lons1(x,y)

