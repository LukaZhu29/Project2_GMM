% 输入数据是GMMrawdataNS
% 'CovarianceType','full'
% 关于收敛类型的讨论： https://chatgpt.com/share/695f0599-db24-8009-bf98-35732692aa99
%% ===================== GMM 主脚本：2D/3D/4D 自适应 =====================
clear; clc;
%[angle Brachioradialis Triceps Biceps] 
% 关闭 TeX 解释器，避免下划线变下标
set(groot,'defaultTextInterpreter','none');
set(groot,'defaultAxesTickLabelInterpreter','none');
set(groot,'defaultLegendInterpreter','none');

% basepath = 'E:\GraduateStudy\UA\OERL Lab\writing\Project2 Writing\Submission\Reviewer comment reply\code\data\GMM_preparedata';
basepath = 'E:\JialeZhu\Paper2 Submission\Reviewer comment reply\code\data\GMM_preparedata';
weight = ["02", "04"];
bpm = ["40", "60", "80"];
task = ["_Isokin_L_", "_Isokin_R_"];
RawData = [];
count = 0;
for i = 1:length(weight)
    for j = 1:length(bpm)
        for k = 1:length(task)
            f1 = fullfile(basepath,'Combined_GMMrawdataNS'+ task(k)+ weight(i) +'kg_'+ bpm(j) +'bpm.xlsx');
            rawData = xlsread(f1);
            RawData = [RawData;rawData];
            count = count+1;
        end
    end
end


RawData(:,2:4) = RawData(:,2:4) * 100;   % EMG -> %MVC 

%% ---------- 选择建模维度 ----------
target_d = 4;        % <<< 改成 2 / 3 / 4
switch target_d
    case 2
        % 2D：Angle + 平均EMG
        data = [RawData(:,1), mean(RawData(:,2:4),2)];
        varNames = {'Angle (deg)','MVC avg (%)'};
    case 3
        % 3D：Angle + 两块肌肉（示例：EMG1、EMG2；需要可改）
        data = [RawData(:,1), mean(RawData(:,[2 4]),2), RawData(:,3)]; % Brachioradialis Biceps; Triceps
        varNames = {'Angle (deg)','EMG1 (%)','EMG2 (%)'};
    case 4
        % 4D：Angle + 三块肌肉
        data = RawData(:,1:4);
        varNames = {'Angle (deg)','EMG1 (%)','EMG2 (%)','EMG3 (%)'};
    otherwise
        error('仅支持 target_d = 2 / 3 / 4');
end
[N, d] = size(data);
% % ---------- ±nσ 过滤  ----------
nSigma = 2;   % 想用 ±2σ 就写 2，±1σ 就写 1
[data, mask_keep] = sigmaFilterND(data, nSigma);
[N, d] = size(data);   % 重新更新 N, d
% 如果后面还要用到 RawData 对应的行，可以：RawData = RawData(mask_keep,:);

% ---------- z-score（仅用于训练） ----------
mu_d = mean(data,1);
sd_d = std(data,[],1); sd_d(sd_d < 1e-12) = 1e-12;
data_z = (data - mu_d) ./ sd_d;

%% ---------- 拟合 GMM（在 z 空间） ----------
K = 3;   % <<< 分量数，可改
gm_z = fitgmdist(data_z, K, 'Start','plus','Replicates',20, ...
    'RegularizationValue',1e-6, ...
    'CovarianceType','full',...
    'Options',statset('MaxIter',1000,'TolFun',1e-6,'TolX',1e-6,'Display','final'));%'Display'设置：'fianl/off'对应打印/不打印 repetition信息
assert(gm_z.Converged, 'EM 未收敛，请提高 MaxIter/Replicates 或正则化');

% ---------- 反变换回原单位（统一 full 协方差） ----------
gm_orig = unscaleGMM(gm_z, mu_d, sd_d);

% fprintf('GMM（原单位）：d=%d, K=%d\n', d, K);
% for k = 1:K
%     fprintf('  Comp %d | w=%.3f | mu=[%s]\n', k, gm_orig.ComponentProportion(k), ...
%         strjoin(arrayfun(@(x) sprintf('%.2f',x), gm_orig.mu(k,:), 'uni',0), ', '));
% end

% ---------- 聚类 ----------
idx = cluster(gm_orig, data);
%% ---------- 基于 Angle(第1维) 的 GMM 分量构造 fuzzy membership + 边界 ----------
angle_dim = 1;

% 在函数内部完成：提取 mu/sigma、按 mu 排序、求相邻等概率交点边界、构建 membership
% fuzzyAngleModel = gaussianEqualPosteriorBoundaries(gm_orig, angle_dim);

% 可选：打印边界
% fprintf('\n=== Angle-dim boundaries (equal-posterior, no weights) ===\n');
% for t = 1:numel(fuzzyAngleModel.boundary)
%     fprintf('Boundary %d (Comp %d vs %d): %.4f deg\n', ...
%         t, t, t+1, fuzzyAngleModel.boundary(t));
% end

% 可选：画 membership + 边界（建议开，便于检查）
% plotMembershipAndBoundaries(fuzzyAngleModel, data(:,angle_dim), varNames{angle_dim});

%% 保存模型 先不保存
% outDir  = 'E:\GraduateStudy\UA\OERL Lab\writing\Project2 Writing\Submission\Reviewer comment reply\code\GMM_model';    % 目标文件夹
outDir  = 'E:\JialeZhu\Paper2 Submission\Reviewer comment reply\code\GMM_model';
if ~exist(outDir,'dir'); mkdir(outDir); end

modelFile = fullfile(outDir, 'GMM_Oct_2_d4.mat');
save(modelFile, 'gm_orig'); 

% save(modelFile, 'gm_orig', 'fuzzyAngleModel');

% dataFile = fullfile(outDir, 'GMMorig_data_2.xlsx');
% GMM_data = [data idx];
% xlswrite(dataFile,GMM_data);

%% ---------- d=3/4：选择"可视化维度" = [1(Angle), vis_dim2] ---------- 先不可视化
if d >= 3
    C = corr(data, 'rows','pairwise');
    [~, mxi] = max(abs(C(1,2:end)));
    vis_dim2 = mxi + 1;            % 与角度相关性最强的那一维
    fprintf('可视化维度选择：始终包含 Dim1=%s；另选 Dim%d=%s（与角度相关性最强，|r|=%.3f）\n', ...
        varNames{1}, vis_dim2, varNames{vis_dim2}, abs(C(1,vis_dim2)));
else
    vis_dim2 = 2;                   % 2D 时第二维就是唯一的EMG
end

% ---------- 导出每个分量各维的一元 mean/std/weight ----------
print_component_marginals(gm_orig, varNames, basepath);

%% ===================== 可视化 =====================

%% 载入保存模型
% basepath = 'E:\GraduateStudy\UA\OERL Lab\Experiment\Project2 GMM\dataprocessing\Jan6_2025\GMM_model';
% filename = fullfile(basepath,'GMM_Oct_2_d4.mat');
% gm = load(filename,'gm_orig');
% gm_orig = gm.gm_orig;

% —— 单变量分布（仅画 Angle 和 选择的 vis_dim2；满宽度）——
figure('Position',[100,100,900,650]);
tiledlayout(2, 1, 'TileSpacing','compact','Padding','compact');

nexttile;
plot_univariate_dist(data(:,1), gm_orig, 1, K, varNames{1});

nexttile;
plot_univariate_dist(data(:,vis_dim2), gm_orig, vis_dim2, K, varNames{vis_dim2});

% —— 二维等高线：始终画 (Angle, vis_dim2) —— 
figure('Position',[100,100,800,600]);
plot_gmm_contour2d(gm_orig, data, [1, vis_dim2], 'mean', varNames, idx)

% —— 可选：d=3 时画 3D 散点与 3σ 椭球 —— 
if d == 3
    figure('Position',[960,100,820,650]);
    colors = lines(K);
    scatter3(data(:,1), data(:,2), data(:,3), 6, idx, '.'); hold on; grid on; box on;
    xlabel(varNames{1}); ylabel(varNames{2}); zlabel(varNames{3});
    mu = gm_orig.mu; plot3(mu(:,1),mu(:,2),mu(:,3),'kp','MarkerSize',14,'MarkerFaceColor','y','LineWidth',1.2);
    title('GMM Scatter (3D)');
    draw_ellipsoid = true;
    if draw_ellipsoid
        for k = 1:K, draw_gaussian_ellipsoid(mu(k,:).', gm_orig.Sigma(:,:,k), 3, colors(k,:)); end
        legend({'Data','Means','3\sigma Ellipsoids'});
    end
end
%% CM plot
% basepath = 'E:\JialeZhu\Paper2 Submission\Reviewer comment reply\code';
% filename = fullfile(basepath,'GMM_model','GMM_Oct_1.mat');
% gm2 = load(filename,'gm_orig');
% gm2 = gm2.gm_orig;

[CMs, ~] = gmm_to_clouds_FCG3(gm_orig, 1, ... %gm_orig
    'NPerComp', 2500, 'ProportionalToWeight', true, ...
    'RngSeed', 42, 'ShowPlot', false);

plot_cloud_scatter(CMs, ...
    'XLabel', 'Elbow_flexion', ...
    'YLabel', 'Membership', ...
    'Title',  '', ...
    'DotSize', 12, ...
    'Alpha',   0.65, ...
    'LegendPrefix', 'CM_', ...
    'Colors', [0.85,0.33,0.10; 0.00,0.45,0.74; 0.93,0.69,0.13],...
    'YLim', [0 1.02]);


%% ===================== 辅助函数 =====================

function gm_orig = unscaleGMM(gm_z, mu_d, sd_d)
% 把 z 空间 GMM 映回原单位；统一返回 full 协方差
    d = numel(mu_d); K = size(gm_z.mu,1); D = diag(sd_d);
    mu_orig = gm_z.mu .* sd_d + mu_d;
    Sigma_orig = zeros(d,d,K);
    sz = size(gm_z.Sigma);
    if numel(sz) >= 3 && sz(2) == 1
        % diagonal: d x 1 x K
        for k = 1:K
            Sigma_orig(:,:,k) = D * diag(gm_z.Sigma(:,1,k)) * D;
        end
    else
        % full: d x d x K
        for k = 1:K
            Sigma_orig(:,:,k) = D * gm_z.Sigma(:,:,k) * D;
        end
    end
    gm_orig = gmdistribution(mu_orig, Sigma_orig, gm_z.ComponentProportion);
end


function plot_univariate_dist(x, gm, dim, K, labelStr)
% 单变量直方图 + 混合密度（原单位）
    h = histogram(x, 'BinMethod','fd', 'Normalization','pdf', ...
        'FaceColor',[0.7 0.7 0.9], 'EdgeColor','none'); hold on; grid on; box on;
%     h.NumBins = 25;
    xr = linspace(min(x), max(x), 1000)'; f = zeros(size(xr));
    for k = 1:K
        mu_k = gm.mu(k,dim);
        s_k  = sqrt(gm.Sigma(dim,dim,k));
        w_k  = gm.ComponentProportion(k);
        f = f + w_k * normpdf(xr, mu_k, s_k);
    end
    plot(xr, f, 'k-', 'LineWidth', 2);
    cols = lines(K);
    for k = 1:K
        mu_k = gm.mu(k,dim);
        s_k  = sqrt(gm.Sigma(dim,dim,k));
        w_k  = gm.ComponentProportion(k);
        plot(xr, w_k*normpdf(xr, mu_k, s_k), '--', 'Color', cols(k,:), 'LineWidth',1.2);
    end
    xlabel(labelStr, 'Interpreter','none'); ylabel('PDF');
    title(['Univariate GMM on ', char(labelStr)], 'Interpreter','none');
end

function plot_gmm_contour2d(gm, data, dims, fixedMode, varNames, idx)
% 在任意维 GMM 上画二维等高线 + 自定义图例（Cluster i / Mean i）
% - dims: 要画的两维（例如 [1 2]）
% - fixedMode: 'mean' 或长度为 d 的向量（≥3D 时用于固定其余维度）
% - idx: 外部传入的聚类标签（建议传入）；若缺省则在函数内计算

    if nargin < 4 || isempty(fixedMode), fixedMode = 'mean'; end
    if nargin < 5 || isempty(varNames)
        varNames = arrayfun(@(j) sprintf('Dim %d',j), 1:size(gm.mu,2), 'uni', 0);
    end
    d = size(gm.mu,2);
    K = size(gm.mu,1);
    i = dims(1); j = dims(2);

    % 网格
    x1 = linspace(min(data(:,i)), max(data(:,i)), 150);
    x2 = linspace(min(data(:,j)), max(data(:,j)), 150);
    [X1, X2] = meshgrid(x1, x2);

    % 计算二维切片上的 PDF
    if d == 2
        P = pdf(gm, [X1(:), X2(:)]);
    else
        if ischar(fixedMode) && strcmpi(fixedMode,'mean')
            xfix = mean(data,1);
        elseif isnumeric(fixedMode) && numel(fixedMode)==d
            xfix = fixedMode;
        else
            error('fixedMode 必须为 ''mean'' 或长度为 d 的数值向量');
        end
        Xfull = repmat(xfix, numel(X1), 1);
        Xfull(:,i) = X1(:); Xfull(:,j) = X2(:);
        P = pdf(gm, Xfull);
    end

    % ---- 绘图（自定义图例）----
    hold on; grid on; box on;
    colors = lines(K);

    % 1) 按簇绘制散点（避免 gscatter 自动图例文字）
    if nargin < 6 || isempty(idx)
        idx = cluster(gm, data);  % 若未传入则内部计算
    end
    hScat = gobjects(K,1);
    for k = 1:K
        mk = (idx == k);
        hScat(k) = scatter(data(mk,i), data(mk,j), 6, ...
                           'Marker', '.', 'MarkerEdgeColor', colors(k,:), ...
                           'DisplayName', sprintf('Cluster %d', k));
    end

    % 2) 轮廓线（不进图例）
%     [~, hContour] = contour(X1, X2, reshape(P, size(X1)), 10, 'LineWidth', 1.2);
%     set(get(get(hContour,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');

    % 3) 分量中心（星形）
    mu = gm.mu;
    hMean = gobjects(K,1);
    for k = 1:K
        hMean(k) = plot(mu(k,i), mu(k,j), 'p', ...
            'MarkerSize', 11, ...
            'MarkerFaceColor', colors(k,:), ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 1.0, ...
            'DisplayName', sprintf('Mean %d', k));
    end

    % 4) 图例：先簇、后中心
    legLabels = [ arrayfun(@(k) sprintf('Cluster %d',k), 1:K, 'uni', false), ...
                  arrayfun(@(k) sprintf('Mean %d',k),    1:K, 'uni', false) ];
    legend([hScat(:); hMean(:)], legLabels, 'Location','best');

    xlabel(varNames{i}); ylabel(varNames{j});
    ylim([0 100]);
    title(sprintf('GMM Contours (%s vs %s)', varNames{i}, varNames{j}));
end


function print_component_marginals(gm, varNames, outdir)
% 打印并导出每个分量在每一维的一元 N(mean, std)
    K = size(gm.mu,1); d = size(gm.mu,2);
    if nargin < 2 || isempty(varNames)
        varNames = arrayfun(@(j) sprintf('Dim %d',j), 1:d, 'uni',0);
    end
    mu_comp  = gm.mu;
    std_comp = zeros(K,d);
    for k = 1:K, std_comp(k,:) = sqrt(diag(gm.Sigma(:,:,k)))'; end
    w = gm.ComponentProportion(:);

    fprintf('\n=== Component-wise marginal N(mean, std) ===\n');
    for k = 1:K
        fprintf('Comp %d | weight = %.3f\n', k, w(k));
        for j = 1:d
            fprintf('  %-12s: mean=%8.3f, std=%8.3f\n', varNames{j}, mu_comp(k,j), std_comp(k,j));
        end
    end

    varNames_mu  = matlab.lang.makeValidName(strcat(varNames, ' mu'));
    varNames_std = matlab.lang.makeValidName(strcat(varNames, ' std'));
    T = table((1:K)', w, 'VariableNames', {'Component','Weight'});
    T = [T, array2table(mu_comp,  'VariableNames', varNames_mu), ...
            array2table(std_comp, 'VariableNames', varNames_std)];
    if nargin >= 3 && ~isempty(outdir)
        try
            writetable(T, fullfile(outdir, 'GMM_component_marginals.xlsx'));
        catch
            writetable(T, fullfile(outdir, 'GMM_component_marginals.csv'));
        end
    end
end

function draw_gaussian_ellipsoid(mu, Sigma, nsig, colorRGB)
% 画 3σ 高斯椭球（仅 3D）
    if numel(mu) ~= 3 || ~all(size(Sigma) == [3 3]), return; end
    [V,D] = eig((Sigma+Sigma')/2);
    radii = nsig * sqrt(max(diag(D),0));
    [x,y,z] = sphere(30);
    E = (V * diag(radii)) * [x(:)'; y(:)'; z(:)'];
    X = reshape(E(1,:)+mu(1), size(x));
    Y = reshape(E(2,:)+mu(2), size(y));
    Z = reshape(E(3,:)+mu(3), size(z));
    s = surf(X,Y,Z);
    set(s,'FaceAlpha',0.08,'EdgeAlpha',0.1,'FaceColor',colorRGB,'EdgeColor',colorRGB);
end

function [data_filt, mask_keep] = sigmaFilterND(data, nSigma)
%SIGMAFILTERND  按 ±nσ 准则对多维数据做简单异常值过滤
%
% INPUT
%   data    : N x D 数据矩阵（N 个样本，D 个特征）
%   nSigma  : 标准差倍数（例如 1 表示 ±1σ, 2 表示 ±2σ）
%
% OUTPUT
%   data_filt : 过滤后的数据（只保留每一维都落在 [μ- nσ σ, μ+ nσ σ] 范围内的样本）
%   mask_keep : N x 1 逻辑向量，true 表示该行被保留

    if nargin < 2
        error('sigmaFilterND 需要两个输入：data 和 nSigma');
    end
    if ~isnumeric(data) || ~ismatrix(data)
        error('data 必须是数值型矩阵');
    end
    if ~isscalar(nSigma) || ~isnumeric(nSigma) || nSigma <= 0
        error('nSigma 必须是正的标量（例如 1、2、3）');
    end

    % 逐列计算均值与标准差
    mu   = mean(data, 1);          % 1 x D
    sigma = std(data, 0, 1);       % 1 x D (无偏 or 有偏都问题不大)
    sigma(sigma < 1e-12) = 1e-12;  % 防止除以 0

    lower = mu - nSigma .* sigma;  % 1 x D
    upper = mu + nSigma .* sigma;  % 1 x D

    % 对于每一行，检查所有维度是否都在 [lower, upper] 之内
    mask_keep = all(data >= lower & data <= upper, 2);   % N x 1 逻辑数组

    data_filt = data(mask_keep, :);
end



