close all
clear
clc

cd /discover/nobackup/drholdaw/tmp.7775/

lat = ncread('x0011dh_a.prog.eta.20130105_00z.nc4','lat');
lon = ncread('x0011dh_a.prog.eta.20130105_00z.nc4','lon');

t = ncread('x0011dh_a.prog.eta.20130105_00z.nc4','tv');
dp = ncread('x0011dh_a.prog.eta.20130105_00z.nc4','delp');

lat_1 = 120;
lat_2 = 240;
lon_1 = 400;
lon_2 = 576;

% lat_1 = 1;
% lat_2 = 361;
% lon_1 = 1;
% lon_2 = 576;

t_trop = t;

p = zeros(576,361,73);
p(:,:,1) = 1;
for k = 2:73
    p(:,:,k) = p(:,:,k-1) + dp(:,:,k-1);
end
p = p*0.01;

ph = 0.5*(p(:,:,1:72) + p(:,:,2:73));

t_mean = zeros(1,72);
t_col = zeros(1,72);
p_mean = zeros(1,72);
p_col = zeros(1,72);

k = 0;
% for i = 1:576
%     for j = 1:361
for i = lon_1:lon_2
    for j = lat_1:lat_2
        
        t_col(1:72) = t(i,j,:);
        t_mean = t_col + t_mean;
        
        p_col(1:72) = ph(i,j,:);
        p_mean = p_col + p_mean;
        
        k = k + 1;
        
        
    end
end

t_mean = t_mean/k;
p_mean = p_mean/k;

plot(t_mean,1:72)
set(gca,'YDir','reverse')
ylim([1 72])

figure
plot(p_mean,1:72)
set(gca,'YDir','reverse')
ylim([1 72])


F3 = zeros(1,72);

for k = 1:72
        
    if ( p_mean(k) >= 775.  && t_mean(k) <=  275. )
     F3(k) = max(-0.016 * p_mean(k) + 13.4, 0.2);
    end
    if ( p_mean(k) >= 825.  && t_mean(k) <=  282. )
     F3(k) = max(0.11 * t_mean(k) - 30.02, 0.2);
    end
    if ( p_mean(k) >= 775.  && p_mean(k) < 825. && t_mean(k) <=  282. && t_mean(k) > 275.)
     F3(k) = min(max(-0.016*p_mean(k) + 0.11 * t_mean(k) - 16.85, 0.2),1.);
    end
    if ( p_mean(k) >= 825.  && t_mean(k) <=  275. )
     F3(k) = 0.2;
    end
    if ( p_mean(k) <= 775.  || t_mean(k) >  282. )
     F3(k) = 1.0;
    end

end

figure
plot(F3,1:72)
set(gca,'YDir','reverse')
ylim([1 72])

F3a = zeros(1,72);

for k = 1:72
    
    if ( p_mean(k) >= 950.  && t_mean(k) >=  285. )         
        F3a(k) = min(0.2 * t_mean(k) - 56, 2.);
    end
    if ( p_mean(k) >= 925.  && t_mean(k) >=  290. ) 
        F3a(k) = min(0.04 * p_mean(k) - 36., 2.);
    end
    if ( p_mean(k) >= 925.  && p_mean(k) < 950. && t_mean(k) >  285. && t_mean(k) < 290.)
        F3a(k) = max(min(0.04*p_mean(k) + 0.2 * t_mean(k) - 94., 2.),1.);
    end
    if ( p_mean(k) >= 950.  && t_mean(k) >=  290. )
        F3a(k) = 2.0;
    end
    
end

figure
plot(F3a,1:72)
set(gca,'YDir','reverse')
ylim([1 72])

