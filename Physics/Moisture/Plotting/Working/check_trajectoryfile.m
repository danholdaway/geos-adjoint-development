close all
clear
clc

cd /discover/nobackup/drholdaw/tmp.22292/traj1

file = 'x0011dh_a.traj.lcv.20130116_0000z.nc4';


time = ncread(file,'time');
lon = ncread(file,'lon');
lat = ncread(file,'lat');
lev = ncread(file,'lev');

U = ncread(file,'U');
V = ncread(file,'V');
TH = ncread(file,'PT');
Q = ncread(file,'QV');
DP = ncread(file,'DP');
QI = ncread(file,'QI');
QL = ncread(file,'QL');
CFLS = ncread(file,'CFLS');
CFCN = ncread(file,'CFCN');

% timemoist = ncread('iau590.traj.M3.lcv.20120317_0000z.nc4','time');
% lonmoist = ncread('iau590.traj.M3.lcv.20120317_0000z.nc4','lon');
% latmoist = ncread('iau590.traj.M3.lcv.20120317_0000z.nc4','lat');
% levmoist = ncread('iau590.traj.M3.lcv.20120317_0000z.nc4','lev');
% 
% Umoist00 = ncread('iau590.traj.M3.lcv.20120317_0000z.nc4','Umoist');
% Vmoist00 = ncread('iau590.traj.M3.lcv.20120317_0000z.nc4','Vmoist');
% THmoist00 = ncread('iau590.traj.M3.lcv.20120317_0000z.nc4','THmoist');
% Qmoist00 = ncread('iau590.traj.M3.lcv.20120317_0000z.nc4','Qmoist');
% KCBLmoist00 = ncread('iau590.traj.M2.lcv.20120317_0000z.nc4','KCBL');
% TSmoist00 = ncread('iau590.traj.M2.lcv.20120317_0000z.nc4','TS');
% FRLAND00 = ncread('iau590.traj.M2.lcv.20120317_0000z.nc4','FRLAND');
% CTOPmoist00 = ncread('iau590.traj.M2.lcv.20120317_0000z.nc4','CTOP');

% Umoist30 = ncread('iau590.traj.M3.lcv.20120317_0030z.nc4','Umoist');
% Vmoist30 = ncread('iau590.traj.M3.lcv.20120317_0030z.nc4','Vmoist');
% THmoist30 = ncread('iau590.traj.M3.lcv.20120317_0030z.nc4','THmoist');
% Qmoist30 = ncread('iau590.traj.M3.lcv.20120317_0030z.nc4','Qmoist');
% KCBLmoist30 = ncread('iau590.traj.M2.lcv.20120317_0030z.nc4','KCBL');
% TSmoist30 = ncread('iau590.traj.M2.lcv.20120317_0030z.nc4','TS');
% FRLAND30 = ncread('iau590.traj.M2.lcv.20120317_0030z.nc4','FRLAND');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/

%P thickness to P
% P = zeros(90,540,73);
% P(:,:,1) = 1;
% for i = 2:73
%     P(:,:,i) = P(:,:,i-1) + DP(:,:,i-1);
% end
% P = P*0.01;

P = U;

% %COUNT NUMBER OF CONVECTIVE POINTS
% CON_DEP = KCBLmoist00 - CTOPmoist00;
% 
% num_convec_points = 0;
% for i = 1:length(lon)
%     for j = 1:length(lat)
%         
%         if CTOPmoist00(i,j) > 0 && CON_DEP(i,j) > 20
%             
%             num_convec_points = num_convec_points + 1;
%             
%         end
%     end
% end
% fprintf('Number of Convective Points = %g \n',num_convec_points)
% fprintf('Number of Convective Points = %g \n',length(lon)*length(lat) - num_convec_points)


plotlevel = 60;

scrsz = get(0,'ScreenSize');
figure('visible','on','Position',[0 scrsz(4) scrsz(3)/2.5 scrsz(4)/2])

% P1 = P(:,1:90,plotlevel);
% P2 = P(:,91:180,plotlevel);
% P3a = P(:,181:225,plotlevel);
% P3b = P(:,226:270,plotlevel);
% P4 = P(:,271:360,plotlevel);
% P5 = P(:,361:450,plotlevel);
% P6 = P(:,451:540,plotlevel);

P1 = P(:,1:180,plotlevel);
P2 = P(:,181:360,plotlevel);
P3a = P(:,361:450,plotlevel);
P3b = P(:,451:540,plotlevel);
P4 = P(:,541:720,plotlevel);
P5 = P(:,721:900,plotlevel);
P6 = P(:,901:1080,plotlevel);

Pnew = [ P5(:,:) rot90(P1(:,:)) rot90(P2(:,:)) P4(:,:)];

[C,h] = contourf(Pnew(:,:));
% [C,h] = contour(P(:,:,plotlevel));
% set(h,'ShowText','off','LevelStep',2);
colorbar
set(gca,'YDir','reverse')

% subplot(2,1,2)
% contourf(lon,lat,Q(:,:,plotlevel)')
% title('ln (Q_{tot} \prime)')
% % xlim([xmin xmax])
% % ylim([ymin ymax])
% colorbar


%column of theta
TH_col = zeros(72,1);

TH_col(:) = TH(50,50,:)




% figure
% contour(TSmoist00)
% colorbar
% 
% figure
% contour(Qmoist00(:,:,plotlevel))
% colorbar
% 
% 
% [a,b] = find(Qmoist00(:,:,70) <= 3.891122e-3  & Qmoist00(:,:,70) >= 3.891120e-3)
% Umoist00(a,b,70)
% Vmoist00(a,b,70)
% THmoist00(a,b,70)
% Qmoist00(a,b,70)
% KCBLmoist00(a,b)
% TSmoist00(a,b)
% FRLAND00(a,b)
% 
% [a,b] = find(Qmoist30(:,:,70) <= 4.516725e-3 & Qmoist30(:,:,70) >= 4.516723e-3)
% 
% Umoist30(a,b,70)
% Vmoist30(a,b,70)
% THmoist30(a,b,70)
% Qmoist30(a,b,70)
% KCBLmoist30(a,b)
% TSmoist30(a,b)
% FRLAND30(a,b)

%u   1.997661     v  -9.247960     TH   281.3824     q  3.8911218E-03 
%kcbl  64.00000     TS   290.6362     frland  2.1308039E-03 

%u -0.2616079     v  -2.511851     TH   289.1126     q  4.5167240E-03
%kcbl  64.00000     TS   289.2204     frland  0.9946468
