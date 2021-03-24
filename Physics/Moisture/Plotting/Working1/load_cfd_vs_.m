close all
clear
clc


% cfd = zeros(1,3607040);
% dq  = zeros(1,3607040);
% i = 1;
% 
% cfd_vs_pdfwidth
% dq_cf = dq; clear dq;
% 
% qcld = zeros(1,3417650);
% dq  = zeros(1,3417650);
% i = 1;
% 
% qcld_vs_pdfwidth
% dq_qcld = dq; clear dq;
% 
% 
% qcid = zeros(1,2810489);
% dq  = zeros(1,2810489);
% i = 1;
% 
% qcid_vs_pdfwidth
% dq_qcid = dq; clear dq;


% cd /discover/nobackup/drholdaw/tmp.22293/sens.20130117.000000/
% cf   = zeros(1,5000000);
% cf_T = zeros(1,5000000);
% cf_q = zeros(1,5000000);
% i = 1;
% 
% cf_vsJacTandq_3
% cf(i:end) = [];
% cf_T(i:end) = [];
% cf_q(i:end) = [];
% 
% cd /discover/nobackup/drholdaw/ExperimentData/Journal_Articles/bacmeister_development
% save cf_vs_JacTandq_dt.mat
% clear


% cd /discover/nobackup/drholdaw/tmp.22293/sens.20130117.000000/
% ql   = zeros(1,5000000);
% ql_T = zeros(1,5000000);
% ql_q = zeros(1,5000000);
% i = 1;
% 
% ql_vsJacTandq_3
% ql(i:end) = [];
% ql_T(i:end) = [];
% ql_q(i:end) = [];
% 
% cd /discover/nobackup/drholdaw/ExperimentData/Journal_Articles/bacmeister_development
% save ql_vs_JacTandq_dt.mat
% clear


cd /discover/nobackup/drholdaw/tmp.22293/sens.20130117.000000/
qi   = zeros(1,5000000);
qi_T = zeros(1,5000000);
qi_q = zeros(1,5000000);
i = 1;

qi_vsJacTandq_3
qi(i:end) = [];
qi_T(i:end) = [];
qi_q(i:end) = [];

cd /discover/nobackup/drholdaw/ExperimentData/Journal_Articles/bacmeister_development
save qi_vs_JacTandq_dt.mat
clear


cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/