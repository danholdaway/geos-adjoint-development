close all
clear
clc

cd /discover/nobackup/drholdaw/tmp.iau590/tmp.8823/sens.20120318.000000/

lon = ncread('fvpert.eta.nc4','lon');
lat = ncread('fvpert.eta.nc4','lat');
lev = ncread('fvpert.eta.nc4','lev');

U_dry  = ncread('fvpert.eta.dry.20120317_06z.nc4','u');
V_dry  = ncread('fvpert.eta.dry.20120317_06z.nc4','v');
TV_dry = ncread('fvpert.eta.dry.20120317_06z.nc4','tv');
DP_dry = ncread('fvpert.eta.dry.20120317_06z.nc4','delp');
QV_dry = ncread('fvpert.eta.dry.20120317_06z.nc4','sphu');
QL_dry = ncread('fvpert.eta.dry.20120317_06z.nc4','qltot');
QI_dry = ncread('fvpert.eta.dry.20120317_06z.nc4','qitot');
O3_dry = ncread('fvpert.eta.dry.20120317_06z.nc4','ozone');

cd /discover/nobackup/drholdaw/tmp.iau590/tmp.31275/sens.20120318.000000/

U_moi  = ncread('fvpert.eta.moi.1.20120317_06z.nc4','u');
V_moi  = ncread('fvpert.eta.moi.1.20120317_06z.nc4','v');
TV_moi = ncread('fvpert.eta.moi.1.20120317_06z.nc4','tv');
DP_moi = ncread('fvpert.eta.moi.1.20120317_06z.nc4','delp');
QV_moi = ncread('fvpert.eta.moi.1.20120317_06z.nc4','sphu');
QL_moi = ncread('fvpert.eta.moi.1.20120317_06z.nc4','qltot');
QI_moi = ncread('fvpert.eta.moi.1.20120317_06z.nc4','qitot');
O3_moi = ncread('fvpert.eta.moi.1.20120317_06z.nc4','ozone');

cd /discover/nobackup/drholdaw/tmp.iau590/tmp.8823/nlm_runs/free/

U_nlm_free  = ncread('iau590.prog.eta.20120317_06z.nc4','u');
V_nlm_free  = ncread('iau590.prog.eta.20120317_06z.nc4','v');   
TV_nlm_free = ncread('iau590.prog.eta.20120317_06z.nc4','tv');
DP_nlm_free = ncread('iau590.prog.eta.20120317_06z.nc4','delp');
QV_nlm_free = ncread('iau590.prog.eta.20120317_06z.nc4','sphu');
QL_nlm_free = ncread('iau590.prog.eta.20120317_06z.nc4','qltot');
QI_nlm_free = ncread('iau590.prog.eta.20120317_06z.nc4','qitot');
O3_nlm_free = ncread('iau590.prog.eta.20120317_06z.nc4','ozone');

cd /discover/nobackup/drholdaw/tmp.iau590/tmp.8823/nlm_runs/del_1.0/

U_nlm_del_0_2  = ncread('iau590.prog.eta.20120317_06z.nc4','u');
V_nlm_del_0_2  = ncread('iau590.prog.eta.20120317_06z.nc4','v');   
TV_nlm_del_0_2 = ncread('iau590.prog.eta.20120317_06z.nc4','tv');
DP_nlm_del_0_2 = ncread('iau590.prog.eta.20120317_06z.nc4','delp');
QV_nlm_del_0_2 = ncread('iau590.prog.eta.20120317_06z.nc4','sphu');
QL_nlm_del_0_2 = ncread('iau590.prog.eta.20120317_06z.nc4','qltot');
QI_nlm_del_0_2 = ncread('iau590.prog.eta.20120317_06z.nc4','qitot');
O3_nlm_del_0_2 = ncread('iau590.prog.eta.20120317_06z.nc4','ozone');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/

U_nlm = U_nlm_del_0_2 - U_nlm_free;
V_nlm = V_nlm_del_0_2 - V_nlm_free;
TV_nlm = TV_nlm_del_0_2 - TV_nlm_free;
DP_nlm = DP_nlm_del_0_2 - DP_nlm_free;
QV_nlm = QV_nlm_del_0_2 - QV_nlm_free;
QL_nlm = QL_nlm_del_0_2 - QL_nlm_free;
QI_nlm = QI_nlm_del_0_2 - QI_nlm_free;
O3_nlm = O3_nlm_del_0_2 - O3_nlm_free;

im = length(lon);
jm = length(lat);
lm = length(lev);


%Weighting
r_earth = 6378.1;
lat_pole = lat+lat(1)*-1;
A = zeros(length(lon),length(lat));
A(1,:) = (2*pi^2*r_earth^2/(length(lon)*length(lat))) * sind(lat_pole);
for i = 2:length(lon)
    A(i,:) = A(1,:);
end

dA = 0.0;

cov_dry = zeros(lm,4);
cov_moi = zeros(lm,4);
var_nlm = zeros(lm,4);
var_dry = zeros(lm,4);
var_moi = zeros(lm,4);
corr_dry = zeros(lm,4);
corr_moi = zeros(lm,4);

lon_min = 1;
lon_max = 288;
lat_min = 64;
lat_max = 128;


for k = 1:length(lev)

    dA = 0.0;
    
    for i = lon_min:lon_max
        for j = lat_min:lat_max

            dA = dA + A(i,j);

            cov_dry(k,1) = cov_dry(k,1) + U_nlm(i,j,k)*U_dry(i,j,k)*A(i,j);
            cov_moi(k,1) = cov_moi(k,1) + U_nlm(i,j,k)*U_moi(i,j,k)*A(i,j);

            var_nlm(k,1) = var_nlm(k,1) + U_nlm(i,j,k)*U_nlm(i,j,k)*A(i,j);
            var_dry(k,1) = var_dry(k,1) + U_dry(i,j,k)*U_dry(i,j,k)*A(i,j);
            var_moi(k,1) = var_moi(k,1) + U_moi(i,j,k)*U_moi(i,j,k)*A(i,j);
            
            cov_dry(k,2) = cov_dry(k,2) + V_nlm(i,j,k)*V_dry(i,j,k)*A(i,j);
            cov_moi(k,2) = cov_moi(k,2) + V_nlm(i,j,k)*V_moi(i,j,k)*A(i,j);

            var_nlm(k,2) = var_nlm(k,2) + V_nlm(i,j,k)*V_nlm(i,j,k)*A(i,j);
            var_dry(k,2) = var_dry(k,2) + V_dry(i,j,k)*V_dry(i,j,k)*A(i,j);
            var_moi(k,2) = var_moi(k,2) + V_moi(i,j,k)*V_moi(i,j,k)*A(i,j);
            
            cov_dry(k,3) = cov_dry(k,3) + TV_nlm(i,j,k)*TV_dry(i,j,k)*A(i,j);
            cov_moi(k,3) = cov_moi(k,3) + TV_nlm(i,j,k)*TV_moi(i,j,k)*A(i,j);

            var_nlm(k,3) = var_nlm(k,3) + TV_nlm(i,j,k)*TV_nlm(i,j,k)*A(i,j);
            var_dry(k,3) = var_dry(k,3) + TV_dry(i,j,k)*TV_dry(i,j,k)*A(i,j);
            var_moi(k,3) = var_moi(k,3) + TV_moi(i,j,k)*TV_moi(i,j,k)*A(i,j);
            
            cov_dry(k,4) = cov_dry(k,4) + QV_nlm(i,j,k)*QV_dry(i,j,k)*A(i,j);
            cov_moi(k,4) = cov_moi(k,4) + QV_nlm(i,j,k)*QV_moi(i,j,k)*A(i,j);

            var_nlm(k,4) = var_nlm(k,4) + QV_nlm(i,j,k)*QV_nlm(i,j,k)*A(i,j);
            var_dry(k,4) = var_dry(k,4) + QV_dry(i,j,k)*QV_dry(i,j,k)*A(i,j);
            var_moi(k,4) = var_moi(k,4) + QV_moi(i,j,k)*QV_moi(i,j,k)*A(i,j);

        end
    end
    
    cov_dry(k,1) = cov_dry(k,1)/dA;
    cov_moi(k,1) = cov_moi(k,1)/dA;
    var_nlm(k,1) = var_nlm(k,1)/dA;
    var_dry(k,1) = var_dry(k,1)/dA;
    var_moi(k,1) = var_moi(k,1)/dA;
        
    cov_dry(k,2) = cov_dry(k,2)/dA;
    cov_moi(k,2) = cov_moi(k,2)/dA;
    var_nlm(k,2) = var_nlm(k,2)/dA;
    var_dry(k,2) = var_dry(k,2)/dA;
    var_moi(k,2) = var_moi(k,2)/dA;
    
    cov_dry(k,3) = cov_dry(k,3)/dA;
    cov_moi(k,3) = cov_moi(k,3)/dA;
    var_nlm(k,3) = var_nlm(k,3)/dA;
    var_dry(k,3) = var_dry(k,3)/dA;
    var_moi(k,3) = var_moi(k,3)/dA;
    
    cov_dry(k,4) = cov_dry(k,4)/dA;
    cov_moi(k,4) = cov_moi(k,4)/dA;
    var_nlm(k,4) = var_nlm(k,4)/dA;
    var_dry(k,4) = var_dry(k,4)/dA;
    var_moi(k,4) = var_moi(k,4)/dA;
    
    corr_dry(k,1) = cov_dry(k,1)/sqrt(var_nlm(k,1)*var_dry(k,1));
    corr_moi(k,1) = cov_moi(k,1)/sqrt(var_nlm(k,1)*var_moi(k,1));
    
    corr_dry(k,2) = cov_dry(k,2)/sqrt(var_nlm(k,2)*var_dry(k,2));
    corr_moi(k,2) = cov_moi(k,2)/sqrt(var_nlm(k,2)*var_moi(k,2));
    
    corr_dry(k,3) = cov_dry(k,3)/sqrt(var_nlm(k,3)*var_dry(k,3));
    corr_moi(k,3) = cov_moi(k,3)/sqrt(var_nlm(k,3)*var_moi(k,3));
    
    corr_dry(k,4) = cov_dry(k,4)/sqrt(var_nlm(k,4)*var_dry(k,4));
    corr_moi(k,4) = cov_moi(k,4)/sqrt(var_nlm(k,4)*var_moi(k,4));
    
end

format short

disp([lev corr_dry(:,1) corr_moi(:,1) corr_dry(:,1)-corr_moi(:,1)])
disp([lev corr_dry(:,2) corr_moi(:,2) corr_dry(:,2)-corr_moi(:,2)])
disp([lev corr_dry(:,3) corr_moi(:,3) corr_dry(:,3)-corr_moi(:,3)])
disp([lev corr_dry(:,4) corr_moi(:,4) corr_dry(:,4)-corr_moi(:,4)])


disp([(1:72)' lev corr_dry(:,1)-corr_moi(:,1) corr_dry(:,2)-corr_moi(:,2) corr_dry(:,3)-corr_moi(:,3) corr_dry(:,4)-corr_moi(:,4)])

mean_dry_corr_qv = mean(corr_dry(37:72,3))
mean_moi_corr_qv = mean(corr_moi(37:72,3))


mean_dry_corr_qv = mean(corr_dry(37:72,4))
mean_moi_corr_qv = mean(corr_moi(37:72,4))


