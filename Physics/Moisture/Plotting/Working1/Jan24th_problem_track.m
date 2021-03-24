close all
clear
clc

% for i = 1:6;
%     
%     h = num2str(i);

    im = 576;
    jm = 361;
    lm = 72;

    LevMin = 30;
    LevMax = lm;
    LonMin = 1;
    LonMax = im;
    LatMin = 1;
    LatMax = jm;
        
    %Load second TLM state to compare.
    cd /discover/nobackup/drholdaw/tmp.2497/sens.20130126.000000/
    file = ['u001_C180.fvpert.eta.20130124_0430z_nocloud.nc4'];

    u_tl2 = ncread(file,'u',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    v_tl2 = ncread(file,'v',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    t_tl2 = ncread(file,'tv',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    q_tl2 = ncread(file,'sphu',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

    cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

    [ max(abs(u_tl2(:))) max(abs(v_tl2(:))) max(abs(t_tl2(:))) max(abs(q_tl2(:)))]
    
    
    %Load second TLM state to compare.
    cd /discover/nobackup/drholdaw/tmp.2497/sens.20130126.000000/
    file = ['u001_C180.fvpert.eta.20130124_0445z_nocloud.nc4'];

    u_tl2 = ncread(file,'u',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    v_tl2 = ncread(file,'v',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    t_tl2 = ncread(file,'tv',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    q_tl2 = ncread(file,'sphu',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

    cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

    [ max(abs(u_tl2(:))) max(abs(v_tl2(:))) max(abs(t_tl2(:))) max(abs(q_tl2(:)))]
    
    
    %Load second TLM state to compare.
    cd /discover/nobackup/drholdaw/tmp.2497/sens.20130126.000000/
    file = ['u001_C180.fvpert.eta.20130124_0500z_nocloud.nc4'];

    u_tl2 = ncread(file,'u',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    v_tl2 = ncread(file,'v',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    t_tl2 = ncread(file,'tv',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    q_tl2 = ncread(file,'sphu',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

    cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

    [ max(abs(u_tl2(:))) max(abs(v_tl2(:))) max(abs(t_tl2(:))) max(abs(q_tl2(:)))]
    
        %Load second TLM state to compare.
    cd /discover/nobackup/drholdaw/tmp.2497/sens.20130126.000000/
    file = ['u001_C180.fvpert.eta.20130124_0515z_nocloud.nc4'];

    u_tl2 = ncread(file,'u',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    v_tl2 = ncread(file,'v',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    t_tl2 = ncread(file,'tv',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    q_tl2 = ncread(file,'sphu',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

    cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

    [ max(abs(u_tl2(:))) max(abs(v_tl2(:))) max(abs(t_tl2(:))) max(abs(q_tl2(:)))]
    
        %Load second TLM state to compare.
    cd /discover/nobackup/drholdaw/tmp.2497/sens.20130126.000000/
    file = ['u001_C180.fvpert.eta.20130124_0530z_nocloud.nc4'];

    u_tl2 = ncread(file,'u',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    v_tl2 = ncread(file,'v',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    t_tl2 = ncread(file,'tv',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    q_tl2 = ncread(file,'sphu',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

    cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

    [ max(abs(u_tl2(:))) max(abs(v_tl2(:))) max(abs(t_tl2(:))) max(abs(q_tl2(:)))]
    
        %Load second TLM state to compare.
    cd /discover/nobackup/drholdaw/tmp.2497/sens.20130126.000000/
    file = ['u001_C180.fvpert.eta.20130124_0545z_nocloud.nc4'];

    u_tl2 = ncread(file,'u',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    v_tl2 = ncread(file,'v',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    t_tl2 = ncread(file,'tv',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    q_tl2 = ncread(file,'sphu',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

    cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

    [ max(abs(u_tl2(:))) max(abs(v_tl2(:))) max(abs(t_tl2(:))) max(abs(q_tl2(:)))]
    
        %Load second TLM state to compare.
    cd /discover/nobackup/drholdaw/tmp.2497/sens.20130126.000000/
    file = ['u001_C180.fvpert.eta.20130124_0600z_nocloud.nc4'];

    u_tl2 = ncread(file,'u',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    v_tl2 = ncread(file,'v',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    t_tl2 = ncread(file,'tv',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    q_tl2 = ncread(file,'sphu',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

    cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

    [ max(abs(u_tl2(:))) max(abs(v_tl2(:))) max(abs(t_tl2(:))) max(abs(q_tl2(:)))]
    
% end