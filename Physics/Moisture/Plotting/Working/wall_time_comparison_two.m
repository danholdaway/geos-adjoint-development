close all
clear
clc

tlm_6_dry = mean([111.3343 - 16.9482]);
tlm_6_moi_nofilt = mean([111.3343 - 7.9987]);
tlm_6_moi = mean([111.3343]);

adm_6_dry = mean([194.9025 - 35.5921]);
adm_6_moi_nofilt = mean([194.9025 - 8.0460]);
adm_6_moi = mean([194.9025]);

% adm_6_dry = mean([221.7672 - 79.0600]);
% adm_6_moi_nofilt = mean([221.7672 - 19.3218]);
% adm_6_moi = mean([221.7672]);

tlm_24_dry = mean([512 512]);
tlm_24_moi = mean([637 706]);
tlm_24_moi_nlmfilt = mean([0 0]);
tlm_24_moi_tlmfilt = mean([650 650]);

adm_24_dry = mean([666.8083 667.9832 666.8995]);
adm_24_moi = mean([774.9695 727.3828 723.3162]);
adm_24_moi_nlmfilt = mean([512 512]);
adm_24_moi_tlmfilt = mean([512 512]);

tlm_percent = 100*((tlm_6_moi_nofilt/tlm_6_dry)-1)
adm_percent = 100*((adm_6_moi_nofilt/adm_6_dry)-1)

tlm_percent = 100*((tlm_6_moi/tlm_6_dry)-1)
adm_percent = 100*((adm_6_moi/adm_6_dry)-1)


%TLM: DRY, MOIST_NOFILT, MOIST
walltime(1,1) =   tlm_6_dry; 
walltime(1,2) =   tlm_6_moi_nofilt; 
walltime(1,3) =   tlm_6_moi; 

%ADM: DRY, MOIST_NOFILT, MOIST
walltime(2,1) =   adm_6_dry;
walltime(2,2) =   adm_6_moi_nofilt;
walltime(2,3) =   adm_6_moi;

% %TLM: DRY, MOIST_NOFILT, MOIST
% walltime(1,1) =   tlm_24_dry; 
% walltime(1,2) =   tlm_24_moi; 
% walltime(1,3) =   tlm_24_moi_nlmfilt; 
% walltime(1,4) =   tlm_24_moi_tlmfilt; 
% 
% %ADM: DRY, MOIST_NOFILT, MOIST
% walltime(2,1) =   adm_24_dry;
% walltime(2,2) =   adm_24_moi;
% walltime(2,3) =   adm_24_moi_nlmfilt;
% walltime(2,4) =   adm_24_moi_tlmfilt;

times = [6 12 18];


bar(walltime(:,:))
colormap summer

% xlim([0 3.6])

fontsize = 14;

title('6 Hour Forecasts Computational Time','FontSize',fontsize,'FontName','TimesNewRoman')

% xlabel('Window Length','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Time (Seconds)','FontSize',fontsize,'FontName','TimesNewRoman')

legend('Dry','Moist','With Filtering','Location','NorthWest')

set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

saveas(gcf,'processor_dist.eps','psc2')

set(gca,'XtickLabel',['Tan Lin Model';'Adjoint Model'])
