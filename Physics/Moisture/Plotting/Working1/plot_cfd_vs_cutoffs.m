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

a = 1:99.0000001/length(cf_T_sort):100;
b = 1:99.0000001/length(ql_T_sort):100;
c = 1:99.0000001/length(qi_T_sort):100;

subplot(1,2,1)
semilogy(a,abs(cf_T_sort))
subplot(1,2,2)
plot(a,abs(cf_T_sort))

figure
subplot(1,2,1)
semilogy(b,abs(ql_T_sort))
subplot(1,2,2)
plot(b,abs(ql_T_sort))

figure
subplot(1,2,1)
semilogy(c,abs(qi_T_sort))
subplot(1,2,2)
plot(c,abs(qi_T_sort))

i_ql = length(ql_sort_T);
i_ql_90 = round(90*i_ql/100);
i_ql_80 = round(80*i_ql/100);
i_ql_70 = round(70*i_ql/100);
i_ql_60 = round(60*i_ql/100);
i_ql_50 = round(50*i_ql/100);
i_ql_25 = round(25*i_ql/100);

fprintf('ql \n')
disp(ql_T_sort(i_ql_90))
disp(ql_T_sort(i_ql_80))
disp(ql_T_sort(i_ql_70))
disp(ql_T_sort(i_ql_60))
disp(ql_T_sort(i_ql_50))
disp(ql_T_sort(i_ql_25))

i_qi = length(qi_sort_T);
i_qi_90 = round(90*i_qi/100);
i_qi_80 = round(80*i_qi/100);
i_qi_70 = round(70*i_qi/100);
i_qi_60 = round(60*i_qi/100);
i_qi_50 = round(50*i_qi/100);
i_qi_25 = round(25*i_qi/100);

fprintf('qi \n')
disp(qi_T_sort(i_qi_90))
disp(qi_T_sort(i_qi_80))
disp(qi_T_sort(i_qi_70))
disp(qi_T_sort(i_qi_60))
disp(qi_T_sort(i_qi_50))
disp(qi_T_sort(i_qi_25))

i_cf = length(cf_sort_T);
i_cf_90 = round(90*i_cf/100);
i_cf_80 = round(80*i_cf/100);
i_cf_70 = round(70*i_cf/100);
i_cf_60 = round(60*i_cf/100);
i_cf_50 = round(50*i_cf/100);
i_cf_25 = round(25*i_cf/100);

fprintf('cf \n')
disp(cf_T_sort(i_cf_90))
disp(cf_T_sort(i_cf_80))
disp(cf_T_sort(i_cf_70))
disp(cf_T_sort(i_cf_60))
disp(cf_T_sort(i_cf_50))
disp(cf_T_sort(i_qi_25))

