function [score_fuzzy, score_classic, score_GMM] = RULA_compare_fun(X, modelpath)
%{
把modelpath加到输入变量中
Date: 2026.3.10

原本是E:\GraduateStudy\UA\OERL Lab\Experiment\Project2 GMM\dataprocessing\Sep4
文件夹中的RULA_compare_fun，使用了以下的修正因子。
EMG_mean = gm2.mu(:,2)'./100;
S0 = [2 1 2];
Sm = EMG_mean - mean(EMG_mean);
alpha = 0.2;
S1 = S0 + alpha.*Sm; %和 CMs 顺序一致（按 Ex 升序）

本次修改时删去了修正因子，直接用[2 1 2]的得分
Date: Jan8
%}
% basepath = 'E:\GraduateStudy\UA\OERL Lab\Experiment\Project2 GMM\dataprocessing\Sep4\data';
% filename = fullfile(basepath,'GMM_model','GMMorig_D2.mat');
% gm2 = load(filename,'gm_orig');
gm2 = load(modelpath,'gm_orig');
gm2 = gm2.gm_orig;

[CMs, ~] = gmm_to_clouds_FCG3(gm2, 1, ... %gm_orig
    'NPerComp', 2500, 'ProportionalToWeight', true, ...
    'RngSeed', 42, 'ShowPlot', false);



score_fuzzy = Fuzzy_RULA3(X, 5);
% disp("A_fuzzy= " + num2str(score_fuzzy.A_soft))
% disp("lower_arm= " + num2str(score_fuzzy.joint_soft(5))) 
score_classic = Classic_RULA(X);

% disp("A_fuzzy= " + num2str(score_classic.A))
% disp("lower_arm= " + num2str(score_classic.joint_scores(5)))  


S1 = [2 1 2];

score_GMM = GMM_RULA_main2(X, CMs, 'LowerArmRiskVals', S1);
% disp("A_fuzzy="+ num2str(score_GMM.A_fuzzy));  
% disp("lower_arm="+ num2str(score_GMM.LowerArm_cont));

