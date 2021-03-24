close all
clear
clc


DateNumber_Start = datenum(2013, 01, 16, 00, 00, 00);
DateNumber_End   = datenum(2013, 01, 17, 00, 00, 00);
DateNumber_Inc   = datenum(2013, 01, 16, 00, 15, 00) - datenum(2013, 01, 16, 00, 00, 00);

jmax = round((DateNumber_End - DateNumber_Start)/DateNumber_Inc);

date_num_vec = zeros(1,jmax);
date_num_vec(1) = DateNumber_Start;
for j = 2:jmax+1;
    
    date_num_vec(j)   = date_num_vec(j-1) + DateNumber_Inc;
    
end

dates = datestr(date_num_vec, 'yyyymmddHHMM');

cd /discover/nobackup/drholdaw/tmp.22292/

for i = 1:jmax+1;
    
    file = ['x0011dh_a.traj.lcv.',dates(i,1:8),'_',dates(i,9:12),'z.nc4'];
    disp(file);
    
    % Read files and average, write to repl.
    
    cd /discover/nobackup/drholdaw/tmp.22292/traj1
        
    U_a = ncread(file,'U');
    V_a = ncread(file,'V');
    PT_a = ncread(file,'PT');
    QV_a = ncread(file,'QV');
    DP_a = ncread(file,'DP');
    QI_a = ncread(file,'QI');
    QL_a = ncread(file,'QL');
    O3_a = ncread(file,'O3');
    QILS_a = ncread(file,'QILS');
    QICN_a = ncread(file,'QICN');
    QLLS_a = ncread(file,'QLLS');
    QLCN_a = ncread(file,'QLCN');
    CFLS_a = ncread(file,'CFLS');
    CFCN_a = ncread(file,'CFCN');
    PS_a = ncread(file,'PS');
    USTAR_a = ncread(file,'USTAR');
    BSTAR_a = ncread(file,'BSTAR');
    ZPBL_a = ncread(file,'ZPBL');
    CM_a = ncread(file,'CM');
    CT_a = ncread(file,'CT');
    CQ_a = ncread(file,'CQ');
    KCBL_a = ncread(file,'KCBL');
    TS_a = ncread(file,'TS');
    
    cd /discover/nobackup/drholdaw/tmp.22292/traj2
        
    U_b = ncread(file,'U');
    V_b = ncread(file,'V');
    PT_b = ncread(file,'PT');
    QV_b = ncread(file,'QV');
    DP_b = ncread(file,'DP');
    QI_b = ncread(file,'QI');
    QL_b = ncread(file,'QL');
    O3_b = ncread(file,'O3');
    QILS_b = ncread(file,'QILS');
    QICN_b = ncread(file,'QICN');
    QLLS_b = ncread(file,'QLLS');
    QLCN_b = ncread(file,'QLCN');
    CFLS_b = ncread(file,'CFLS');
    CFCN_b = ncread(file,'CFCN');
    PS_b = ncread(file,'PS');
    USTAR_b = ncread(file,'USTAR');
    BSTAR_b = ncread(file,'BSTAR');
    ZPBL_b = ncread(file,'ZPBL');
    CM_b = ncread(file,'CM');
    CT_b = ncread(file,'CT');
    CQ_b = ncread(file,'CQ');
    KCBL_b = ncread(file,'KCBL');
    TS_b = ncread(file,'TS');
   
    %AVERAGE THE TRAJECTORIES
    U_c = 0.5*(U_a + U_b);
    V_c = 0.5*(V_a + V_b);
    PT_c = 0.5*(PT_a + PT_b);
    QV_c = 0.5*(QV_a + QV_b);
    DP_c = 0.5*(DP_a + DP_b);
    QI_c = 0.5*(QI_a + QI_b);
    QL_c = 0.5*(QL_a + QL_b);
    O3_c = 0.5*(O3_a + O3_b);
    QILS_c = 0.5*(QILS_a + QILS_b);
    QICN_c = 0.5*(QICN_a + QICN_b);
    QLLS_c = 0.5*(QLLS_a + QLLS_b);
    QLCN_c = 0.5*(QLCN_a + QLCN_b);
    CFLS_c = 0.5*(CFLS_a + CFLS_b);
    CFCN_c = 0.5*(CFCN_a + CFCN_b);
    PS_c = 0.5*(PS_a + PS_b);
    USTAR_c = 0.5*(USTAR_a + USTAR_b);
    BSTAR_c = 0.5*(BSTAR_a + BSTAR_b);
    ZPBL_c = 0.5*(ZPBL_a + ZPBL_b);
    CM_c = 0.5*(CM_a + CM_b);
    CT_c = 0.5*(CT_a + CT_b);
    CQ_c = 0.5*(CQ_a + CQ_b);
    KCBL_c = 0.5*(KCBL_a + KCBL_b);
    TS_c = 0.5*(TS_a + TS_b);
    
    cd /discover/nobackup/drholdaw/tmp.22292/trajq
    ncwrite(file,'U',U_c);
    ncwrite(file,'V',V_c);
    ncwrite(file,'PT',PT_c);
    ncwrite(file,'QV',QV_c);
    ncwrite(file,'DP',DP_c);
    ncwrite(file,'QI',QI_c);
    ncwrite(file,'QL',QL_c);
    ncwrite(file,'O3',O3_c);
    ncwrite(file,'QILS',QILS_c);
    ncwrite(file,'QICN',QICN_c);
    ncwrite(file,'QLLS',QLLS_c);
    ncwrite(file,'QLCN',QLCN_c);
    ncwrite(file,'CFLS',CFLS_c);
    ncwrite(file,'CFCN',CFCN_c);
    ncwrite(file,'PS',PS_c);
    ncwrite(file,'USTAR',USTAR_c);
    ncwrite(file,'BSTAR',BSTAR_c);
    ncwrite(file,'ZPBL',ZPBL_c);
    ncwrite(file,'CM',CM_c);
    ncwrite(file,'CT',CT_c);
    ncwrite(file,'CQ',CQ_c);
    ncwrite(file,'KCBL',KCBL_c);
    ncwrite(file,'TS',TS_c);
    
    fprintf('WRITE DONE \n\n\n\n')    

end

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

