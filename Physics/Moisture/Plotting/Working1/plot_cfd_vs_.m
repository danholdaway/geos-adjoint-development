close all
clear
clc



cd /discover/nobackup/drholdaw/ExperimentData/Journal_Articles/bacmeister_development/

load cf_vs_JacTandq_dt.mat
load qi_vs_JacTandq_dt.mat
load ql_vs_JacTandq_dt.mat
clear i

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/


[cf_T_sort, ind_cf] = sort(abs(cf_T),'ascend');
cf_sort_T = cf(ind_cf);
cf_sort_q_tmp = cf_q(ind_cf);

[ql_T_sort, ind_ql] = sort(abs(ql_T),'ascend');
ql_sort_T = ql(ind_ql);
ql_sort_q_tmp = ql_q(ind_ql);

[qi_T_sort, ind_qi] = sort(abs(qi_T),'ascend');
qi_sort_T = qi(ind_qi);
qi_sort_q_tmp = qi_q(ind_qi);

cut = 1;
space = 10;

%Remove from cf
i_cf = length(cf_sort_T);
cf_ind_new = [(1:cut-1)' ; (cut:space:i_cf-cut)' ; (i_cf-cut-1:i_cf)' ];

cf_sort_T = cf_sort_T(cf_ind_new);
cf_T_sort = cf_T_sort(cf_ind_new);

cf_sort_q_tmp = cf_sort_q_tmp(cf_ind_new);

[cf_q_sort, ind_cf] = sort(abs(cf_sort_q_tmp),'ascend');
cf_sort_q = cf_sort_T(ind_cf);

%Remove from ql
i_ql = length(ql_sort_T);
ql_ind_new = [(1:cut-1)' ; (cut:space:i_ql-cut)' ; (i_ql-cut-1:i_ql)' ];

ql_sort_T = ql_sort_T(ql_ind_new);
ql_T_sort = ql_T_sort(ql_ind_new);

ql_sort_q_tmp = ql_sort_q_tmp(ql_ind_new);

[ql_q_sort, ind_ql] = sort(abs(ql_sort_q_tmp),'ascend');
ql_sort_q = ql_sort_T(ind_ql);

%Remove from qi
i_qi = length(qi_sort_T);
qi_ind_new = [(1:cut-1)' ; (cut:space:i_qi-cut)' ; (i_qi-cut-1:i_qi)' ];

qi_sort_T = qi_sort_T(qi_ind_new);
qi_T_sort = qi_T_sort(qi_ind_new);

qi_sort_q_tmp = qi_sort_q_tmp(qi_ind_new);

[qi_q_sort, ind_qi] = sort(abs(qi_sort_q_tmp),'ascend');
qi_sort_q = qi_sort_T(ind_qi);





cfd_max = max(abs([cf_sort_T cf_sort_q]));
qld_max = max(abs([ql_sort_T ql_sort_q]));
qid_max = max(abs([qi_sort_T qi_sort_q]));



a = 90;
b = 80;
c = 60;

line_wid = 1.0;
fontsize = 12;

i_ql = length(ql_sort_T);
i_ql_a = round(a*i_ql/100);
i_ql_b = round(b*i_ql/100);
i_ql_c = round(c*i_ql/100);
i_qi = length(qi_sort_T);
i_qi_a = round(a*i_qi/100);
i_qi_b = round(b*i_qi/100);
i_qi_c = round(c*i_qi/100);
i_cf = length(cf_sort_T);
i_cf_a = round(a*i_cf/100);
i_cf_b = round(b*i_cf/100);
i_cf_c = round(c*i_cf/100);


figure
set(gcf,'position',[129 308 1150 611])

subplot(2,3,1)
scatter(ql_sort_T(i_ql_a+1:end),abs(ql_T_sort(i_ql_a+1:end)),'rx')
hold on
scatter(ql_sort_T(i_ql_b+1:i_ql_a),abs(ql_T_sort(i_ql_b+1:i_ql_a)),'bx')
scatter(ql_sort_T(i_ql_c+1:i_ql_b),abs(ql_T_sort(i_ql_c+1:i_ql_b)),'gx')
scatter(ql_sort_T(1:i_ql_c),abs(ql_T_sort(1:i_ql_c)),'kx')
set(gca,'YScale','log')
ylabel('|\partial H / \partial T|','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('\partial q_{l,LS}\prime / \partial t','FontSize',fontsize,'FontName','TimesNewRoman')
title('(a)','FontSize',fontsize,'FontName','TimesNewRoman')
box on
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([-1.1*qld_max 1.1*qld_max])
% ylim([0.5*min(abs(ql_T_sort)) 2*max(abs(ql_T_sort))])
ylim([10^-12 2*max(abs(ql_T_sort))])
legend('> 90%','80% - 90%','60% - 80%','< 60%','location','SouthEast')

subplot(2,3,4)
scatter(ql_sort_q(i_ql_a+1:end),abs(ql_q_sort(i_ql_a+1:end)),'rx')
hold on
scatter(ql_sort_q(i_ql_b+1:i_ql_a),abs(ql_q_sort(i_ql_b+1:i_ql_a)),'bx')
scatter(ql_sort_T(i_ql_c+1:i_ql_b),abs(ql_q_sort(i_ql_c+1:i_ql_b)),'gx')
scatter(ql_sort_T(1:i_ql_c),abs(ql_q_sort(1:i_ql_c)),'kx')
set(gca,'YScale','log')
ylabel('|\partial H / \partial q|','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('\partial q_{l,LS}\prime / \partial t','FontSize',fontsize,'FontName','TimesNewRoman')
title('(d)','FontSize',fontsize,'FontName','TimesNewRoman')
box on
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([-1.1*qld_max 1.1*qld_max])
% ylim([0.5*min(abs(ql_q_sort)) 2*max(abs(ql_q_sort))])
ylim([10^-8 2*max(abs(ql_q_sort))])



subplot(2,3,2)
scatter(qi_sort_T(i_ql_a+1:end),abs(qi_T_sort(i_ql_a+1:end)),'rx')
hold on
scatter(qi_sort_T(i_ql_b+1:i_ql_a),abs(qi_T_sort(i_ql_b+1:i_ql_a)),'bx')
scatter(qi_sort_T(i_ql_c+1:i_ql_b),abs(qi_T_sort(i_ql_c+1:i_ql_b)),'gx')
scatter(qi_sort_T(1:i_ql_c),abs(qi_T_sort(1:i_ql_c)),'kx')
set(gca,'YScale','log')
ylabel('|\partial I / \partial T|','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('\partial q_{i,LS}\prime / \partial t','FontSize',fontsize,'FontName','TimesNewRoman')
title('(b)','FontSize',fontsize,'FontName','TimesNewRoman')
box on
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([-1.1*qid_max 1.1*qid_max])
ylim([0.5*min(abs(qi_T_sort)) 2*max(abs(qi_T_sort))])

subplot(2,3,5)
scatter(qi_sort_q(i_ql_a+1:end),abs(qi_q_sort(i_ql_a+1:end)),'rx')
hold on
scatter(qi_sort_q(i_ql_b+1:i_ql_a),abs(qi_q_sort(i_ql_b+1:i_ql_a)),'bx')
scatter(qi_sort_q(i_ql_c+1:i_ql_b),abs(qi_q_sort(i_ql_c+1:i_ql_b)),'gx')
scatter(qi_sort_q(1:i_ql_c),abs(qi_q_sort(1:i_ql_c)),'kx')
set(gca,'YScale','log')
ylabel('|\partial I / \partial q|','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('\partial q_{l,LS}\prime / \partial t','FontSize',fontsize,'FontName','TimesNewRoman')
title('(e)','FontSize',fontsize,'FontName','TimesNewRoman')
box on
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([-1.1*qid_max 1.1*qid_max])
ylim([0.5*min(abs(qi_q_sort)) 2*max(abs(qi_q_sort))])


subplot(2,3,3)
scatter(cf_sort_T(i_cf_a+1:end),abs(cf_T_sort(i_cf_a+1:end)),'rx')
hold on
scatter(cf_sort_T(i_cf_b+1:i_cf_a),abs(cf_T_sort(i_cf_b+1:i_cf_a)),'bx')
scatter(cf_sort_T(i_cf_c+1:i_cf_b),abs(cf_T_sort(i_cf_c+1:i_cf_b)),'gx')
scatter(cf_sort_T(1:i_cf_c),abs(cf_T_sort(1:i_cf_c)),'kx')
set(gca,'YScale','log')
ylabel('|\partial J/ \partial T|','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('\partial C_{LS}\prime / \partial t','FontSize',fontsize,'FontName','TimesNewRoman')
title('(c)','FontSize',fontsize,'FontName','TimesNewRoman')
box on
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([-1.1*cfd_max 1.1*cfd_max])
ylim([0.5*min(abs(cf_T_sort)) 2*max(abs(cf_T_sort))])

subplot(2,3,6)
scatter(cf_sort_q(i_cf_a+1:end),abs(cf_q_sort(i_cf_a+1:end)),'rx')
hold on
scatter(cf_sort_q(i_cf_b+1:i_cf_a),abs(cf_q_sort(i_cf_b+1:i_cf_a)),'bx')
scatter(cf_sort_q(i_cf_c+1:i_cf_b),abs(cf_q_sort(i_cf_c+1:i_cf_b)),'gx')
scatter(cf_sort_q(1:i_cf_c),abs(cf_q_sort(1:i_cf_c)),'kx')
set(gca,'YScale','log')
ylabel('|\partial J/ \partial q|','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('\partial C_{LS}\prime / \partial t','FontSize',fontsize,'FontName','TimesNewRoman')
title('(f)','FontSize',fontsize,'FontName','TimesNewRoman')
box on
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([-1.1*cfd_max 1.1*cfd_max])
ylim([0.5*min(abs(cf_q_sort)) 2*max(abs(cf_q_sort))])


