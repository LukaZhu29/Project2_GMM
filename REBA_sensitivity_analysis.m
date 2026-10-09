%{
比较BI REBA, fuzzy REBA, REBA的sensitivity
Date: 2025.1.8
GPT 链接：https://chatgpt.com/share/6960598d-779c-8009-bc0f-9814ad247002
%}

clear;clc;
basepath = 'E:\GraduateStudy\UA\OERL Lab\Experiment\Project2 GMM\dataprocessing\Jan6_2025';
filename1 = fullfile(basepath,'GMM_model','REBA_score.xlsx ');

X = xlsread(filename1);
%%
% X is n×4: [theta, gmm, fuzzy, reba], theta = (0:1:120)'
res = rebaSensitivityAnalysis(X);

% 或者用分位数法自动设置epsilon + 保存图片
% opts = struct('epsilon_mode','quantile','q_ref',95,'q_frac',0.01, ...
%               'use_semilogy',true,'save_figs',true,'save_dir','./figs', ...
%               'fig_prefix','elbow');
% res = rebaSensitivityAnalysis(X, opts);
%%
opts = struct();
opts.use_semilogy = false;          % 主图用线性 + ycap（更好看）
opts.zoom_inset = true;             % 开启局部放大 inset
opts.zoom_xwins = [50 70; 90 110];  % 放大 60° 与 100°附近
opts.zoom_use_log_inset = true;     % inset 用对数轴显示 REBA spike
opts.epsilon_mode = 'fixed';
opts.epsilon = 0.01;


res = rebaSensitivityAnalysis2(X, opts);
%% 局部放大，并更改放大图的位置
% 这个很不错
opts = struct();
opts.use_semilogy = false;          % 主图用线性 + ycap（更好看）
opts.zoom_inset = true;             % 开启局部放大 inset
opts.zoom_xwins = [50 70; 90 110];  % 放大 60° 与 100°附近
opts.zoom_use_log_inset = true;     % inset 用对数轴显示 REBA spike
opts.epsilon_mode = 'fixed';
opts.epsilon = 0.01;
opts.inset_pos = [0.14 0.58 0.30 0.33];
opts.inset_transparent = true;

res = rebaSensitivityAnalysis3(X, opts);

%% 把放大图放在右边 
% 最后用的是这个
opts = struct();
opts.use_semilogy        = false;   % 主图线性
% opts.zoom_use_log_panel = true;    % 右侧 panel 用 log
opts.zoom_xwins         = [55 65; 95 105];
opts.zoom_ymax_main     = 0.05;    % 让 GMM / fuzzy 清晰
opts.epsilon            = 0.01;

res = rebaSensitivityAnalysis4(X, opts);
