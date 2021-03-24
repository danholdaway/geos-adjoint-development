close all
clear
clc

fprintf('DRY DRY vs MOIST DRY \n');

aircraft(1) = -0.23732850894;
aircraft(2) = -0.243593189812;
fprintf('Aircraft & %5.2f \\\\ \n', aircraft(1)/aircraft(2)*100)

amsua(1) = -0.419498295773;
amsua(2) = -0.44613092998;
fprintf('AMSUA & %5.2f \\\\ \n', amsua(1)/amsua(2)*100)

aqua_airs(1) = -0.202647488792;
aqua_airs(2) = -0.21398466592;
fprintf('Aqua AIRS & %5.2f \\\\ \n', aqua_airs(1)/aqua_airs(2)*100)

ascat_wind(1) = -0.0100655663125;
ascat_wind(2) = -0.0135975421875;
fprintf('ASCAT Wind & %5.2f \\\\ \n', ascat_wind(1)/ascat_wind(2)*100)

dropsondes(1) = -0.00120289841703;
dropsondes(2) = -0.000881743796658;
fprintf('Dropsonde & %5.2f \\\\ \n', dropsondes(1)/dropsondes(2)*100)

gpsro(1) = -0.0773009581708;
gpsro(2) = -0.0896491298544;
fprintf('GPSRO & %5.2f \\\\ \n', gpsro(1)/gpsro(2)*100)

hirs(1) = -0.0331718994097;
hirs(2) = -0.0366425337472;
fprintf('HIRS & %5.2f \\\\ \n', hirs(1)/hirs(2)*100)

iasi(1) = -0.228633066939;
iasi(2) = -0.234282538936;
fprintf('IASI & %5.2f \\\\ \n', iasi(1)/iasi(2)*100)

land_sfc(1) = -0.0255188297813;
land_sfc(2) = -0.033798884875;
fprintf('Land-Surface & %5.2f \\\\ \n', land_sfc(1)/land_sfc(2)*100)

marine_sfc(1) = -0.029556969879;
marine_sfc(2) = -0.0341383216356;
fprintf('Marine-Surface & %5.2f \\\\ \n', marine_sfc(1)/marine_sfc(2)*100)

mhs(1) = -0.00221482517156;
mhs(2) = -0.00450690739469;
fprintf('MHS & %5.2f \\\\ \n', mhs(1)/mhs(2)*100)

modis_wind(1) = -0.000169086425281;
modis_wind(2) = -0.000452767265772;
fprintf('MODIS Wind & %5.2f \\\\ \n', modis_wind(1)/modis_wind(2)*100)

nexrad_wind(1) = -0.00241851854715;
nexrad_wind(2) = -0.00216787510566;
fprintf('NEXRAD Wind & %5.2f \\\\ \n', nexrad_wind(1)/nexrad_wind(2)*100)

pibal(1) = -0.00825011785769;
pibal(2) = -0.00938577889941;
fprintf('PIBAL & %5.2f \\\\ \n', pibal(1)/pibal(2)*100)

profiler_wind(1) = -0.00122863183546;
profiler_wind(2) = -0.00110551360747;
fprintf('Profiler Wind & %5.2f \\\\ \n', profiler_wind(1)/profiler_wind(2)*100)

radiosondes(1) = -0.40962654201;
radiosondes(2) = -0.427633866883;
fprintf('Radiosonde & %5.2f \\\\ \n', radiosondes(1)/radiosondes(2)*100)

satellite_wind(1) = -0.172065570701;
satellite_wind(2) = -0.18471370266;
fprintf('Satellite Wind & %5.2f \\\\ \n', satellite_wind(1)/satellite_wind(2)*100)

tmi_rain_rate(1) = 6.39276289655e-05;
tmi_rain_rate(2) = 0.00014934726;
fprintf('TMI Rain Rate & %5.2f \\\\ \n', tmi_rain_rate(1)/tmi_rain_rate(2)*100)

windsat_wind(1) = -0.00557256496875;
windsat_wind(2) = -0.0070739305;
fprintf('WINDSAT Wind & %5.2f \\\\ \n', windsat_wind(1)/windsat_wind(2)*100)

total_dry = aircraft(1)+amsua(1)+aqua_airs(1)+ascat_wind(1)+dropsondes(1)+gpsro(1)+hirs(1)+iasi(1)+land_sfc(1)+marine_sfc(1)+mhs(1)+modis_wind(1)+nexrad_wind(1)+pibal(1)+profiler_wind(1)+radiosondes(1)+satellite_wind(1)+tmi_rain_rate(1)+windsat_wind(1);
total_wet = aircraft(2)+amsua(2)+aqua_airs(2)+ascat_wind(2)+dropsondes(2)+gpsro(2)+hirs(2)+iasi(2)+land_sfc(2)+marine_sfc(2)+mhs(2)+modis_wind(2)+nexrad_wind(2)+pibal(2)+profiler_wind(2)+radiosondes(2)+satellite_wind(2)+tmi_rain_rate(2)+windsat_wind(2);

fprintf('TOTAL & %5.2f \\\\ \n \n \n \n', total_dry/total_wet*100)


%%%%% DRY DRY vs DRY MOIST %%%%%%
fprintf('DRY DRY vs DRY MOIST \n');


aircraft(1) = -0.23732850894;
aircraft(2) = -0.252019137916;
amsua(1) = -0.419498295773;
amsua(2) = -0.464194324135;
aqua_airs(1) = -0.202647488792;
aqua_airs(2) = -0.282047920302;
ascat_wind(1) = -0.0100655663125;
ascat_wind(2) = -0.0119284119062;
dropsondes(1) = -0.00120289841703;
dropsondes(2) = -0.00153453611784;
gpsro(1) = -0.0773009581708;
gpsro(2) = -0.0887280744042;
hirs(1) = -0.0331718994097;
hirs(2) = -0.0519835615256;
iasi(1) = -0.228633066939;
iasi(2) = -0.284707587227;
land_sfc(1) = -0.0255188297813;
land_sfc(2) = -0.0279189185531;
marine_sfc(1) = -0.029556969879;
marine_sfc(2) = -0.0319866527344;
mhs(1) = -0.00221482517156;
mhs(2) = -0.0136134302297;
modis_wind(1) = -0.000169086425281;
modis_wind(2) = -0.000214427279375;
nexrad_wind(1) = -0.00241851854715;
nexrad_wind(2) = -0.00273942782823;
pibal(1) = -0.00825011785769;
pibal(2) = -0.0131241391731;
profiler_wind(1) = -0.00122863183546;
profiler_wind(2) = -0.00144743675181;
radiosondes(1) = -0.40962654201;
radiosondes(2) = -0.458410797254;
satellite_wind(1) = -0.172065570701;
satellite_wind(2) = -0.190107135343;
tmi_rain_rate(1) = 6.39276289655e-05;
tmi_rain_rate(2) = 0.0001355401;
windsat_wind(1) = -0.00557256496875;
windsat_wind(2) = -0.00658377926875;

total_dry = aircraft(1)+amsua(1)+aqua_airs(1)+ascat_wind(1)+dropsondes(1)+gpsro(1)+hirs(1)+iasi(1)+land_sfc(1)+marine_sfc(1)+mhs(1)+modis_wind(1)+nexrad_wind(1)+pibal(1)+profiler_wind(1)+radiosondes(1)+satellite_wind(1)+tmi_rain_rate(1)+windsat_wind(1);
total_wet = aircraft(2)+amsua(2)+aqua_airs(2)+ascat_wind(2)+dropsondes(2)+gpsro(2)+hirs(2)+iasi(2)+land_sfc(2)+marine_sfc(2)+mhs(2)+modis_wind(2)+nexrad_wind(2)+pibal(2)+profiler_wind(2)+radiosondes(2)+satellite_wind(2)+tmi_rain_rate(2)+windsat_wind(2);

fprintf('Aircraft & %5.2f \\\\ \n', aircraft(1)/aircraft(2)*100)

fprintf('AMSUA & %5.2f \\\\ \n', amsua(1)/amsua(2)*100)

fprintf('Aqua AIRS & %5.2f \\\\ \n', aqua_airs(1)/aqua_airs(2)*100)

fprintf('ASCAT Wind & %5.2f \\\\ \n', ascat_wind(1)/ascat_wind(2)*100)

fprintf('Dropsonde & %5.2f \\\\ \n', dropsondes(1)/dropsondes(2)*100)

fprintf('GPSRO & %5.2f \\\\ \n', gpsro(1)/gpsro(2)*100)

fprintf('HIRS & %5.2f \\\\ \n', hirs(1)/hirs(2)*100)

fprintf('IASI & %5.2f \\\\ \n', iasi(1)/iasi(2)*100)

fprintf('Land-Surface & %5.2f \\\\ \n', land_sfc(1)/land_sfc(2)*100)

fprintf('Marine-Surface & %5.2f \\\\ \n', marine_sfc(1)/marine_sfc(2)*100)

fprintf('MHS & %5.2f \\\\ \n', mhs(1)/mhs(2)*100)

fprintf('MODIS Wind & %5.2f \\\\ \n', modis_wind(1)/modis_wind(2)*100)

fprintf('NEXRAD Wind & %5.2f \\\\ \n', nexrad_wind(1)/nexrad_wind(2)*100)

fprintf('PIBAL & %5.2f \\\\ \n', pibal(1)/pibal(2)*100)

fprintf('Profiler Wind & %5.2f \\\\ \n', profiler_wind(1)/profiler_wind(2)*100)

fprintf('Radiosonde & %5.2f \\\\ \n', radiosondes(1)/radiosondes(2)*100)

fprintf('Satellite Wind & %5.2f \\\\ \n', satellite_wind(1)/satellite_wind(2)*100)

fprintf('TMI Rain Rate & %5.2f \\\\ \n', tmi_rain_rate(1)/tmi_rain_rate(2)*100)

fprintf('WINDSAT Wind & %5.2f \\\\ \n', windsat_wind(1)/windsat_wind(2)*100)

fprintf('TOTAL & %5.2f \\\\ \n \n \n', total_dry/total_wet*100)


fprintf('DRY DRY vs MOIST MOIST \n');

aircraft(1) = -0.23732850894;
aircraft(2) = -0.258654494891;
amsua(1) = -0.419498295773;
amsua(2) = -0.480946766403;
aqua_airs(1) = -0.202647488792;
aqua_airs(2) = -0.269659417722;
ascat_wind(1) = -0.0100655663125;
ascat_wind(2) = -0.0156813567188;
dropsondes(1) = -0.00120289841703;
dropsondes(2) = -0.00128869377819;
gpsro(1) = -0.0773009581708;
gpsro(2) = -0.104160689572;
hirs(1) = -0.0331718994097;
hirs(2) = -0.050446745695;
iasi(1) = -0.228633066939;
iasi(2) = -0.273445181351;
land_sfc(1) = -0.0255188297813;
land_sfc(2) = -0.0372624484062;
marine_sfc(1) = -0.029556969879;
marine_sfc(2) = -0.0367790518798;
mhs(1) = -0.00221482517156;
mhs(2) = -0.0133067356819;
modis_wind(1) = -0.000169086425281;
modis_wind(2) = -0.000496473717469;
nexrad_wind(1) = -0.00241851854715;
nexrad_wind(2) = -0.00242035623978;
pibal(1) = -0.00825011785769;
pibal(2) = -0.0138858388028;
profiler_wind(1) = -0.00122863183546;
profiler_wind(2) = -0.00130691235188;
radiosondes(1) = -0.40962654201;
radiosondes(2) = -0.484167585014;
satellite_wind(1) = -0.172065570701;
satellite_wind(2) = -0.203889845549;
tmi_rain_rate(1) = 6.39276289655e-05;
tmi_rain_rate(2) = 0.000177236876667;
windsat_wind(1) = -0.00557256496875;
windsat_wind(2) = -0.00810228900313;

total_dry = aircraft(1)+amsua(1)+aqua_airs(1)+ascat_wind(1)+dropsondes(1)+gpsro(1)+hirs(1)+iasi(1)+land_sfc(1)+marine_sfc(1)+mhs(1)+modis_wind(1)+nexrad_wind(1)+pibal(1)+profiler_wind(1)+radiosondes(1)+satellite_wind(1)+tmi_rain_rate(1)+windsat_wind(1);
total_wet = aircraft(2)+amsua(2)+aqua_airs(2)+ascat_wind(2)+dropsondes(2)+gpsro(2)+hirs(2)+iasi(2)+land_sfc(2)+marine_sfc(2)+mhs(2)+modis_wind(2)+nexrad_wind(2)+pibal(2)+profiler_wind(2)+radiosondes(2)+satellite_wind(2)+tmi_rain_rate(2)+windsat_wind(2);


fprintf('Aircraft & %5.2f \\\\ \n', aircraft(1)/aircraft(2)*100)

fprintf('AMSUA & %5.2f \\\\ \n', amsua(1)/amsua(2)*100)

fprintf('Aqua AIRS & %5.2f \\\\ \n', aqua_airs(1)/aqua_airs(2)*100)

fprintf('ASCAT Wind & %5.2f \\\\ \n', ascat_wind(1)/ascat_wind(2)*100)

fprintf('Dropsonde & %5.2f \\\\ \n', dropsondes(1)/dropsondes(2)*100)

fprintf('GPSRO & %5.2f \\\\ \n', gpsro(1)/gpsro(2)*100)

fprintf('HIRS & %5.2f \\\\ \n', hirs(1)/hirs(2)*100)

fprintf('IASI & %5.2f \\\\ \n', iasi(1)/iasi(2)*100)

fprintf('Land-Surface & %5.2f \\\\ \n', land_sfc(1)/land_sfc(2)*100)

fprintf('Marine-Surface & %5.2f \\\\ \n', marine_sfc(1)/marine_sfc(2)*100)

fprintf('MHS & %5.2f \\\\ \n', mhs(1)/mhs(2)*100)

fprintf('MODIS Wind & %5.2f \\\\ \n', modis_wind(1)/modis_wind(2)*100)

fprintf('NEXRAD Wind & %5.2f \\\\ \n', nexrad_wind(1)/nexrad_wind(2)*100)

fprintf('PIBAL & %5.2f \\\\ \n', pibal(1)/pibal(2)*100)

fprintf('Profiler Wind & %5.2f \\\\ \n', profiler_wind(1)/profiler_wind(2)*100)

fprintf('Radiosonde & %5.2f \\\\ \n', radiosondes(1)/radiosondes(2)*100)

fprintf('Satellite Wind & %5.2f \\\\ \n', satellite_wind(1)/satellite_wind(2)*100)

fprintf('TMI Rain Rate & %5.2f \\\\ \n', tmi_rain_rate(1)/tmi_rain_rate(2)*100)

fprintf('WINDSAT Wind & %5.2f \\\\ \n', windsat_wind(1)/windsat_wind(2)*100)

fprintf('TOTAL & %5.2f \\\\ \n \n \n \n', total_dry/total_wet*100)










%%%%% DRY MOIST vs MOIST MOIST %%%%%%

fprintf('DRY MOIST vs MOIST MOIST \n');

aircraft(1) = -0.252019137916;
aircraft(2) = -0.258654494891;
amsua(1) = -0.464194324135;
amsua(2) = -0.480946766403;
aqua_airs(1) = -0.282047920302;
aqua_airs(2) = -0.269659417722;
ascat_wind(1) = -0.0119284119062;
ascat_wind(2) = -0.0156813567188;
dropsondes(1) = -0.00153453611784;
dropsondes(2) = -0.00128869377819;
gpsro(1) = -0.0887280744042;
gpsro(2) = -0.104160689572;
hirs(1) = -0.0519835615256;
hirs(2) = -0.050446745695;
iasi(1) = -0.284707587227;
iasi(2) = -0.273445181351;
land_sfc(1) = -0.0279189185531;
land_sfc(2) = -0.0372624484062;
marine_sfc(1) = -0.0319866527344;
marine_sfc(2) = -0.0367790518798;
mhs(1) = -0.0136134302297;
mhs(2) = -0.0133067356819;
modis_wind(1) = -0.000214427279375;
modis_wind(2) = -0.000496473717469;
nexrad_wind(1) = -0.00273942782823;
nexrad_wind(2) = -0.00242035623978;
pibal(1) = -0.0131241391731;
pibal(2) = -0.0138858388028;
profiler_wind(1) = -0.00144743675181;
profiler_wind(2) = -0.00130691235188;
radiosondes(1) = -0.458410797254;
radiosondes(2) = -0.484167585014;
satellite_wind(1) = -0.190107135343;
satellite_wind(2) = -0.203889845549;
tmi_rain_rate(1) = 0.0001355401;
tmi_rain_rate(2) = 0.000177236876667;
windsat_wind(1) = -0.00658377926875;
windsat_wind(2) = -0.00810228900313;

total_dry = aircraft(1)+amsua(1)+aqua_airs(1)+ascat_wind(1)+dropsondes(1)+gpsro(1)+hirs(1)+iasi(1)+land_sfc(1)+marine_sfc(1)+mhs(1)+modis_wind(1)+nexrad_wind(1)+pibal(1)+profiler_wind(1)+radiosondes(1)+satellite_wind(1)+tmi_rain_rate(1)+windsat_wind(1);
total_wet = aircraft(2)+amsua(2)+aqua_airs(2)+ascat_wind(2)+dropsondes(2)+gpsro(2)+hirs(2)+iasi(2)+land_sfc(2)+marine_sfc(2)+mhs(2)+modis_wind(2)+nexrad_wind(2)+pibal(2)+profiler_wind(2)+radiosondes(2)+satellite_wind(2)+tmi_rain_rate(2)+windsat_wind(2);

fprintf('Aircraft & %5.2f \\\\ \n', aircraft(1)/aircraft(2)*100)

fprintf('AMSUA & %5.2f \\\\ \n', amsua(1)/amsua(2)*100)

fprintf('Aqua AIRS & %5.2f \\\\ \n', aqua_airs(1)/aqua_airs(2)*100)

fprintf('ASCAT Wind & %5.2f \\\\ \n', ascat_wind(1)/ascat_wind(2)*100)

fprintf('Dropsonde & %5.2f \\\\ \n', dropsondes(1)/dropsondes(2)*100)

fprintf('GPSRO & %5.2f \\\\ \n', gpsro(1)/gpsro(2)*100)

fprintf('HIRS & %5.2f \\\\ \n', hirs(1)/hirs(2)*100)

fprintf('IASI & %5.2f \\\\ \n', iasi(1)/iasi(2)*100)

fprintf('Land-Surface & %5.2f \\\\ \n', land_sfc(1)/land_sfc(2)*100)

fprintf('Marine-Surface & %5.2f \\\\ \n', marine_sfc(1)/marine_sfc(2)*100)

fprintf('MHS & %5.2f \\\\ \n', mhs(1)/mhs(2)*100)

fprintf('MODIS Wind & %5.2f \\\\ \n', modis_wind(1)/modis_wind(2)*100)

fprintf('NEXRAD Wind & %5.2f \\\\ \n', nexrad_wind(1)/nexrad_wind(2)*100)

fprintf('PIBAL & %5.2f \\\\ \n', pibal(1)/pibal(2)*100)

fprintf('Profiler Wind & %5.2f \\\\ \n', profiler_wind(1)/profiler_wind(2)*100)

fprintf('Radiosonde & %5.2f \\\\ \n', radiosondes(1)/radiosondes(2)*100)

fprintf('Satellite Wind & %5.2f \\\\ \n', satellite_wind(1)/satellite_wind(2)*100)

fprintf('TMI Rain Rate & %5.2f \\\\ \n', tmi_rain_rate(1)/tmi_rain_rate(2)*100)

fprintf('WINDSAT Wind & %5.2f \\\\ \n', windsat_wind(1)/windsat_wind(2)*100)

fprintf('TOTAL & %5.2f \\\\ \n \n \n ', total_dry/total_wet*100)
