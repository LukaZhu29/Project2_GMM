%{
Jan7_RULA_compare_simulate.m
Data: Jan7
Objective: compare different RULA
Functions:
%}
% 比较3种方法的RULA得分图

clear;clc;
basepath = 'E:\JialeZhu\Paper2 Submission\Reviewer comment reply\code';
% modelpath = fullfile(basepath,'GMM_model','GMM_Oct_1.mat');  
modelpath = fullfile(basepath,'GMM_model','GMM_Oct_2_d4.mat');

%% single joint angle
% single_angle = 30;
% X1 = zeros(1,6);
% X1(:,5) = single_angle;
% [~, ~, score_GMM1] = RULA_compare_fun(X1)
%% 计算不同REBA score
% num = size(Rawdata1,1);
test_length = 120;
X = zeros(test_length,6);
test = 1:test_length;
X(:,5) = test';
[score_fuzzy, score_classic, score_GMM] = RULA_compare_fun(X,modelpath);

% 导出REBA score
elbow_GMM = score_GMM.LowerArm_cont(1:test_length);
elbow_fuzzy = score_fuzzy.joint_soft(1:test_length,5);
elbow_classic = score_classic.joint_scores(1:test_length,5);
REBA_score = [test' elbow_GMM elbow_fuzzy elbow_classic];

% outDir  = 'E:\GraduateStudy\UA\OERL Lab\Experiment\Project2 GMM\dataprocessing\Jan6_2025\GMM_model';    % 目标文件夹
% dataFile = fullfile(outDir, 'REBA_score.xlsx');
% xlswrite(dataFile,REBA_score);


%% 分开画 以及 合并画
len = 1:test_length;

% % 论文里用的是这个顺序
% figure('Position',[100,100,1500,300])
% subplot(1,3,3)
% plot(score_GMM.LowerArm_cont(len),'LineWidth',3.5,'Color',[0 161 241]./255);
% xlabel('elbow angle');ylabel('risk score')
% xlim([0 120]);ylim([0.8 2.2])
% xticks([0 20 40 60 80 100 120])
% yticks([1 1.2 1.4 1.6 1.8 2])
% subplot(1,3,2)
% plot(score_fuzzy.joint_soft(len,5),'LineWidth',3.5,'Color',[246,83,20]./255);
% xlabel('elbow angle');ylabel('risk score')
% xlim([0 120]);ylim([0.8 2.2])
% xticks([0 60 100 180])
% yticks([1 1.2 1.4 1.6 1.8 2])
% subplot(1,3,1)
% plot(score_classic.joint_scores(len,5),'LineWidth',3.5,'Color',[124 187 0]./255);
% xlabel('elbow angle');ylabel('risk score')
% xlim([0 120]);ylim([0.8 2.2])
% xticks([0 60 100 180])
% yticks([1 1.2 1.4 1.6 1.8 2])

% 分开画
% y_lim = [0.8 2.5];
figure('Position',[100,100,1500,300])
subplot(1,3,1)
plot(score_GMM.LowerArm_cont(len),'LineWidth',2,'Color',[0 161 241]./255);
xlabel('elbow angle');ylabel('risk score')
xlim([0 120])
xticks([0 20 40 60 80 100 120])
subplot(1,3,2)
plot(score_fuzzy.joint_soft(len,5),'LineWidth',2,'Color',[246,83,20]./255);
xlabel('elbow angle');ylabel('risk score')
xlim([0 180])
xticks([0 60 100 180])
subplot(1,3,3)
plot(score_classic.joint_scores(len,5),'LineWidth',2,'Color',[124 187 0]./255);
xlabel('elbow angle');ylabel('risk score')
xlim([0 180])
xticks([0 60 100 180])


%%

% 分成3张图
len = 1:test_length;
% y_lim = [0.8 2.5];
figure('Position',[100,100,400,300])
plot(score_GMM.LowerArm_cont(len),'LineWidth',2,'Color',[0 161 241]./255);
xlabel('elbow angle');ylabel('risk score')
xlim([0 120]); xticks([0 20 40 60 80 100 120])
ylim([0.9 2.1]); yticks([1 1.5 2])
figure('Position',[100,100,400,300])
plot(score_fuzzy.joint_soft(len,5),'LineWidth',2,'Color',[246,83,20]./255);
xlabel('elbow angle');ylabel('risk score')
xlim([0 180]); xticks([0 60 100 180])
ylim([0.9 2.1]); yticks([1 1.5 2])
figure('Position',[100,100,400,300])
plot(score_classic.joint_scores(len,5),'LineWidth',2,'Color',[124 187 0]./255);
xlabel('elbow angle');ylabel('risk score')
xlim([0 180]); xticks([0 60 100 180])
ylim([0.9 2.1]); yticks([1 1.5 2])
%% 比较peak torque
len = 1:test_length;
GMM_score = score_GMM.LowerArm_cont(len);
fuzzy_score = score_fuzzy.joint_soft(len,5);
classic_score = score_classic.joint_scores(len,5);
NGMM_score = GMM_score./max(GMM_score);
Nfuzzy_score = fuzzy_score./max(fuzzy_score);
Nclassic_score = classic_score./max(classic_score);

peakTorqueE = [40.2 50.5 52.8 50.7 45.4;
    22.1 27.6 29.2 27.7 25.1]; % Extension
peakTorqueF = [56.4 62.9 63.0 60.3 56.4;
    26.2 30.4 32.0 32.4 30.4]; % Flexion
peakTorque0 = [peakTorqueE; peakTorqueF];
peakTorque = mean(peakTorque0,1)';
NpeakTorque = peakTorque./max(peakTorque);
AnglePoint = [15 40 65 90 110];

% 比较peak torque
Tfit = interp1(AnglePoint, NpeakTorque, len, 'pchip', 'extrap'); % Tfit: Torque fit
figure
plot(NGMM_score,'LineWidth',2,'Color',[0 161 241]./255); hold on;
plot(Nfuzzy_score,'LineWidth',2,'Color',[246,83,20]./255); hold on;
plot(Nclassic_score,'LineWidth',2,'Color',[124 187 0]./255); hold on;
plot(Tfit,'k','LineWidth',2); 
legend("Normalized SEF-REBA score", "Normalized Fuzzy REBA score", "Normalized REBA score", "Normalized peak torque")
ylim([0.2 1.5]);
xlabel("Elbow angle (degree)");
ylabel("Normalized value")

corrResult = zeros(3,2);
% [rho,p] = corr(NGMM_score, Tfit', 'Type', 'Spearman')
[corrResult(1,1) corrResult(1,2)] = corr(NGMM_score, Tfit', 'Type', 'Spearman');
[corrResult(2,1) corrResult(2,2)] = corr(Nfuzzy_score, Tfit', 'Type', 'Spearman');
[corrResult(3,1) corrResult(3,2)] = corr(Nclassic_score, Tfit', 'Type', 'Spearman');
corrResult