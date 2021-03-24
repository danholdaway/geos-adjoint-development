close all
clear
clc


load coast
coast_lat = lat; clear lat
coast_lon = long; clear long
grey = 0.0;

cd /discover/nobackup/drholdaw/tmp.11265/

lat = ncread('iau590a.traj.lcv.20120410_1300z.nc4','lat');
lon = ncread('iau590a.traj.lcv.20120410_1300z.nc4','lon');

u1 = ncread('iau590a.traj.lcv.20120410_1300z.nc4','U');
v1 = ncread('iau590a.traj.lcv.20120410_1300z.nc4','V');
t1 = ncread('iau590a.traj.lcv.20120410_1300z.nc4','PT');
q1 = ncread('iau590a.traj.lcv.20120410_1300z.nc4','QV');
p2 = ncread('iau590a.traj.lcv.20120410_1300z.nc4','DP');

p1 = zeros(90,540,72);
p1(:,:,1) = 1;
for i = 2:72
    p1(:,:,i) = p1(:,:,i-1)+p2(:,:,i-1);
end

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/

plotlevel = 72;
lon_plot = -180:179;
lat_plot = -90:89;

% figure
% contourf(u1(:,:,plotlevel),'LineStyle','none')
% colorbar
% title('u')
% 
% asd

U1 = u1(:,1:90,plotlevel);
U2 = u1(:,91:180,plotlevel);
U3 = u1(:,181:270,plotlevel);
U4 = u1(:,271:360,plotlevel);
U5 = u1(:,361:450,plotlevel);
U6 = u1(:,451:540,plotlevel);

V1 = v1(:,1:90,plotlevel);
V2 = v1(:,91:180,plotlevel);
V3 = v1(:,181:270,plotlevel);
V4 = v1(:,271:360,plotlevel);
V5 = v1(:,361:450,plotlevel);
V6 = v1(:,451:540,plotlevel);

T1 = t1(:,1:90,plotlevel);
T2 = t1(:,91:180,plotlevel);
T3 = t1(:,181:270,plotlevel);
T4 = t1(:,271:360,plotlevel);
T5 = t1(:,361:450,plotlevel);
T6 = t1(:,451:540,plotlevel);

Q1 = q1(:,1:90,plotlevel);
Q2 = q1(:,91:180,plotlevel);
Q3 = q1(:,181:270,plotlevel);
Q4 = q1(:,271:360,plotlevel);
Q5 = q1(:,361:450,plotlevel);
Q6 = q1(:,451:540,plotlevel);

P1 = p1(:,1:90,plotlevel);
P2 = p1(:,91:180,plotlevel);
P3 = p1(:,181:270,plotlevel); %North Pole
P4 = p1(:,271:360,plotlevel);
P5 = p1(:,361:450,plotlevel);
P6 = p1(:,451:540,plotlevel); %South Pole

U_Npole = zeros(45,360);
U_Spole = zeros(45,360);
U_Npole(:,1:90) = rot90(U3(:,46:end),3);
U_Npole(:,91:180) = rot90(U3(1:45,:),2);
U_Npole(:,181:270) = rot90(U3(:,1:45));
U_Npole(:,271:end) = U3(46:end,:);
U_Spole(:,1:90) = U6(1:45,:);
U_Spole(:,91:180) = rot90(U6(:,46:end),1);
U_Spole(:,181:270) = rot90(U6(46:end,:),2);
U_Spole(:,271:end) = rot90(U6(:,1:45),3);

V_Npole = zeros(45,360);
V_Spole = zeros(45,360);
V_Npole(:,1:90) = rot90(V3(:,46:end),3);
V_Npole(:,91:180) = rot90(V3(1:45,:),2);
V_Npole(:,181:270) = rot90(V3(:,1:45));
V_Npole(:,271:end) = V3(46:end,:);
V_Spole(:,1:90) = V6(1:45,:);
V_Spole(:,91:180) = rot90(V6(:,46:end),1);
V_Spole(:,181:270) = rot90(V6(46:end,:),2);
V_Spole(:,271:end) = rot90(V6(:,1:45),3);

T_Npole = zeros(45,360);
T_Spole = zeros(45,360);
T_Npole(:,1:90) = rot90(T3(:,46:end),3);
T_Npole(:,91:180) = rot90(T3(1:45,:),2);
T_Npole(:,181:270) = rot90(T3(:,1:45));
T_Npole(:,271:end) = T3(46:end,:);
T_Spole(:,1:90) = T6(1:45,:);
T_Spole(:,91:180) = rot90(T6(:,46:end),1);
T_Spole(:,181:270) = rot90(T6(46:end,:),2);
T_Spole(:,271:end) = rot90(T6(:,1:45),3);

Q_Npole = zeros(45,360);
Q_Spole = zeros(45,360);
Q_Npole(:,1:90) = rot90(Q3(:,46:end),3);
Q_Npole(:,91:180) = rot90(Q3(1:45,:),2);
Q_Npole(:,181:270) = rot90(Q3(:,1:45));
Q_Npole(:,271:end) = Q3(46:end,:);
Q_Spole(:,1:90) = Q6(1:45,:);
Q_Spole(:,91:180) = rot90(Q6(:,46:end),1);
Q_Spole(:,181:270) = rot90(Q6(46:end,:),2);
Q_Spole(:,271:end) = rot90(Q6(:,1:45),3);

P_Npole = zeros(45,360);
P_Spole = zeros(45,360);
P_Npole(:,1:90) = rot90(P3(:,46:end),3);
P_Npole(:,91:180) = rot90(P3(1:45,:),2);
P_Npole(:,181:270) = rot90(P3(:,1:45));
P_Npole(:,271:end) = P3(46:end,:);
P_Spole(:,1:90) = P6(1:45,:);
P_Spole(:,91:180) = rot90(P6(:,46:end),1);
P_Spole(:,181:270) = rot90(P6(46:end,:),2);
P_Spole(:,271:end) = rot90(P6(:,1:45),3);

u = flipud([U_Npole; [ U5(:,:) rot90(U1(:,:)) rot90(U2(:,:)) U4(:,:)]; U_Spole]);
v = flipud([V_Npole; [ V5(:,:) rot90(V1(:,:)) rot90(V2(:,:)) V4(:,:)]; V_Spole]);
t = flipud([T_Npole; [ T5(:,:) rot90(T1(:,:)) rot90(T2(:,:)) T4(:,:)]; T_Spole]);
q = flipud([Q_Npole; [ Q5(:,:) rot90(Q1(:,:)) rot90(Q2(:,:)) Q4(:,:)]; Q_Spole]);
p = 0.01*flipud([P_Npole; [ P5(:,:) rot90(P1(:,:)) rot90(P2(:,:)) P4(:,:)]; P_Spole]);

u = [u(:,end-36:end) u(:,1:end-37)];
v = [v(:,end-36:end) v(:,1:end-37)];
t = [t(:,end-36:end) t(:,1:end-37)];
q = [q(:,end-36:end) q(:,1:end-37)];
p = [p(:,end-36:end) p(:,1:end-37)];

figure('position',[151 236 1128 683])
subplot(2,2,1)
contourf(lon_plot,lat_plot,u(:,:),'LineStyle','none')
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('u')

subplot(2,2,2)
contourf(lon_plot,lat_plot,v(:,:),'LineStyle','none')
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('v')

subplot(2,2,3)
contourf(lon_plot,lat_plot,t(:,:),'LineStyle','none')
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('t')

subplot(2,2,4)
contourf(lon_plot,lat_plot,q(:,:),'LineStyle','none')
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('q')

figure
contourf(lon_plot,lat_plot,u(:,:),'LineStyle','none')
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('u')

figure
contourf(lon_plot,lat_plot,p(:,:),'LineStyle','none','LevelStep',5)
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('p')