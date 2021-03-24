close all
clear
clc

tlm_6_dry = mean([110.8190 - 16.6329]);
tlm_6_moi = mean([110.8190]);

adm_6_dry = mean([193.1074 - 33.6018]);
adm_6_moi = mean([193.1074]);

tlm_24_dry = mean([383.9399 - 74.0506]);
tlm_24_moi = mean([383.9399]);


adm_24_dry = mean([666.8083 667.9832 666.8995]);
adm_24_moi = mean([774.9695 727.3828 723.3162]);


tlm_percent = 100*((tlm_6_moi/tlm_6_dry)-1)
adm_percent = 100*((adm_6_moi/adm_6_dry)-1)


%TLM: DRY, MOIST_NOFILT, MOIST
% walltime(1,1) =   tlm_6_dry; 
% walltime(1,2) =   tlm_6_moi; 

walltime(1,1) =   tlm_24_dry; 
walltime(1,2) =   tlm_24_moi; 

%ADM: DRY, MOIST_NOFILT, MOIST
walltime(2,1) =   adm_6_dry;
walltime(2,2) =   adm_6_moi;


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
