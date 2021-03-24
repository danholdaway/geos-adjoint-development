close all
clear
clc


DateNumber_Start = datenum(2013, 01, 16, 06, 00, 00);
DateNumber_End   = datenum(2013, 01, 17, 00, 00, 00);
DateNumber_Inc1   = datenum(2013, 01, 16, 00, 07, 00) - datenum(2013, 01, 16, 00, 00, 00);
DateNumber_Inc2   = datenum(2013, 01, 16, 00, 15, 00) - datenum(2013, 01, 16, 00, 07, 00);

jmax = round((DateNumber_End - DateNumber_Start)/(0.5*(DateNumber_Inc1 + DateNumber_Inc2))) ;

date_num_vec = zeros(1,jmax);
date_num_vec(1) = DateNumber_Start;
for j = 2:2:jmax;
    
    date_num_vec(j)   = date_num_vec(j-1) + DateNumber_Inc1;
    date_num_vec(j+1) = date_num_vec(j-1) + DateNumber_Inc1 + DateNumber_Inc2;
    
end


dates = datestr(date_num_vec, 'yyyymmddHHMM')

cd /discover/nobackup/drholdaw/tmp.22292/

for i = 1:jmax;
    
    file1 = ['x0011dh_a.traj.lcv.',dates(i,1:8),'_',dates(i,9:12),'z.nc4']
    file2 = ['x0011dh_a.traj.lcv.',dates(i+1,1:8),'_',dates(i+1,9:12),'z.nc4']
    
    % Read File 2 and overwrite file 1 with contents of file 2.
    
    PTM = ncread(file2,'PTM');
    QVM = ncread(file2,'QVM');
    KCBL = ncread(file2,'KCBL');
    TS = ncread(file2,'TS');
    CTOP = ncread(file2,'CTOP');
    SEEDRAS = ncread(file2,'SEEDRAS');
    
%     PTM1 = ncread(file1,'PTM');
%     QVM1 = ncread(file1,'QVM');
%     KCBL1 = ncread(file1,'KCBL');
%     TS1 = ncread(file1,'TS');
%     CTOP1 = ncread(file1,'CTOP');
%     SEEDRAS1 = ncread(file1,'SEEDRAS');
    
    ncwrite(file1,'PTM',PTM);
    ncwrite(file1,'QVM',QVM);
    ncwrite(file1,'KCBL',KCBL);
    ncwrite(file1,'TS',TS);
    ncwrite(file1,'CTOP',CTOP);
    ncwrite(file1,'SEEDRAS',SEEDRAS);
    
   fprintf('WRITE DONE \n\n\n\n')    

end




cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/