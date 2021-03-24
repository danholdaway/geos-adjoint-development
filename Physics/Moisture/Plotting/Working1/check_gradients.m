close all
clear
clc

date = datenum(2013, 01, 03, 00, 00, 00);

% start = '15';
start = '21';
% start = '00';
lead = 24;

plot_lev = 65;

plot_var = 2;

%Get required dates
fprintf('Date is %s \n', datestr(date, 'yyyymmddHHMM'))

dates = datestr(date, 'yyyymmddHHMM');
if strcmp(start,'00') == 1
    dater = datestr(date, 'yyyymmddHHMM');
else
    dater = datestr(date - (24 - str2double(start))/24, 'yyyymmddHHMM');    
end
datee  = datestr(date + lead/24, 'yyyymmddHHMM');


cd /home/drholdaw/Lin_Moist_Physics/Inputs/
pref
    

dir = ['/archive/u/drholdaw/u001_C180_PH0/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H',start,'/'];
cd(dir)
file = ['u001_C180_PH0.fsens_txe.eta.',dater(1:8),'_',dater(9:10),'z+',datee(1:8),'_00z-',dates(1:8),'_00z.nc4'];

lon = ncread(file,'lon');
lat = ncread(file,'lat');
lev = ncread(file,'lev');

u_ph0 = ncread(file,'u');
v_ph0 = ncread(file,'v');
t_ph0 = ncread(file,'tv');
q_ph0 = ncread(file,'sphu');

dir = ['/archive/u/drholdaw/u001_C180_PH1/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H',start,'/'];
cd(dir)
file = ['u001_C180_PH1.fsens_txe.eta.',dater(1:8),'_',dater(9:10),'z+',datee(1:8),'_00z-',dates(1:8),'_00z.nc4'];

u_ph1 = ncread(file,'u');
v_ph1 = ncread(file,'v');
t_ph1 = ncread(file,'tv');
q_ph1 = ncread(file,'sphu');

dir = ['/archive/u/drholdaw/u001_C180_PH2/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H',start,'/'];
cd(dir)
file = ['u001_C180_PH2.fsens_txe.eta.',dater(1:8),'_',dater(9:10),'z+',datee(1:8),'_00z-',dates(1:8),'_00z.nc4'];

u_ph2 = ncread(file,'u');
v_ph2 = ncread(file,'v');
t_ph2 = ncread(file,'tv');
q_ph2 = ncread(file,'sphu');

dir = ['/archive/u/drholdaw/u001_C180_PH2w/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H',start,'/'];
cd(dir)
file = ['u001_C180_PH2w.fsens_twe.eta.',dater(1:8),'_',dater(9:10),'z+',datee(1:8),'_00z-',dates(1:8),'_00z.nc4'];

u_ph2w = ncread(file,'u');
v_ph2w = ncread(file,'v');
t_ph2w = ncread(file,'tv');
q_ph2w = ncread(file,'sphu');


cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/


if plot_var == 1
    ph0 = u_ph0(:,:,plot_lev);
    ph1 = u_ph1(:,:,plot_lev);
    ph2 = u_ph2(:,:,plot_lev);
    ph2w = u_ph2w(:,:,plot_lev);
elseif plot_var == 2
    ph0 = v_ph0(:,:,plot_lev);
    ph1 = v_ph1(:,:,plot_lev);
    ph2 = v_ph2(:,:,plot_lev);
    ph2w = v_ph2w(:,:,plot_lev);
elseif plot_var == 3
    ph0 = t_ph0(:,:,plot_lev);
    ph1 = t_ph1(:,:,plot_lev);
    ph2 = t_ph2(:,:,plot_lev);
    ph2w = t_ph2w(:,:,plot_lev);
elseif plot_var == 4
    ph0 = q_ph0(:,:,plot_lev);
    ph1 = q_ph1(:,:,plot_lev);
    ph2 = q_ph2(:,:,plot_lev);
    ph2w = q_ph2w(:,:,plot_lev);
end


figure
set(gcf,'position',[71 58 1173 837])

subplot(2,2,1)
contourf(lon,lat,ph0','LineStyle','none')
colorbar
caxis([-max(abs(ph0(:))) max(abs(ph0(:)))])

subplot(2,2,2)
contourf(lon,lat,ph1','LineStyle','none')
colorbar
caxis([-max(abs(ph1(:))) max(abs(ph1(:)))])

subplot(2,2,3)
contourf(lon,lat,ph2','LineStyle','none')
colorbar
caxis([-max(abs(ph2(:))) max(abs(ph2(:)))])

subplot(2,2,4)
contourf(lon,lat,ph2w','LineStyle','none')
colorbar
caxis([-max(abs(ph2w(:))) max(abs(ph2w(:)))])


