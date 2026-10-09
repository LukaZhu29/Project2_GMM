% [CMs, orderIdx] = gmm_to_clouds_FCG(gm_orig, 1, ...
%     'NPerComp', 4000, ...            % 每个分量大约云滴数
%     'ProportionalToWeight', true, ...% 按GMM权重分配云滴总量
%     'RngSeed', 42, ...               % 固定随机种子可复现
%     'ShowPlot', true);               % 直接画出结果

basepath = 'E:\JialeZhu\Paper2 Submission\Reviewer comment reply\code\GMM_model';
filename = fullfile(basepath,'GMM_Oct_2_d4.mat');
gm2 = load(filename,'gm_orig');
gm2 = gm2.gm_orig;
%%
[CMs, ~] = gmm_to_clouds_FCG3(gm2, 1, ... %gm_orig
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
%% 
% 可视化GMM
f1 = fullfile(basepath,'GMMorig_data_2.xlsx');
GMM_data = xlsread(f1);
data = GMM_data(:,1:2);
idx = GMM_data(:,3);
K = 3;
varNames = {'Angle (deg)','MVC avg (%)'};
vis_dim2 = 2;

% —— 单变量分布（仅画 Angle 和 选择的 vis_dim2；满宽度）——
figure('Position',[100,100,900,650]);
tiledlayout(2, 1, 'TileSpacing','compact','Padding','compact');

nexttile;
plot_univariate_dist(data(:,1), gm2, 1, K, varNames{1});

nexttile;
plot_univariate_dist(data(:,vis_dim2), gm2, vis_dim2, K, varNames{vis_dim2});

% —— 二维等高线：始终画 (Angle, vis_dim2) —— 
figure('Position',[100,100,800,600]);
plot_gmm_contour2d(gm2, data, [1, vis_dim2], 'mean', varNames, idx)

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
% ylim([0 55])
%% ===================== 辅助函数 =====================
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