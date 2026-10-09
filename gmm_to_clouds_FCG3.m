function [CMs, orderIdx] = gmm_to_clouds_FCG3(gm, dim, varargin)
%GMM_TO_CLOUDS_FCG  将 GMM 的某一维(如第1维角度)各高斯分量映射为云模型并用FCG生成云滴
%
% INPUTS
%   gm        : fitgmdist 对象（K个分量，D维）
%   dim       : 目标维度索引（角度维，通常为1）
%   参数 (Name-Value):
%     'NPerComp'            - 每个分量生成的云滴数，默认 3000
%     'ProportionalToWeight'- 是否按 GMM 权重分配云滴数，默认 true
%     'AlphaClip'           - 将 alpha 裁剪到 [0,1]，默认 true
%     'RngSeed'             - 随机数种子（固定可复现），默认 []
%     'ShowPlot'            - 是否绘图，默认 true
%     'GridX'               - 会员度曲线评估网格（向量），默认自动
%
% OUTPUTS
%   CMs: 结构体数组(1×K)，每个元素含字段：
%        .k, .weight, .mu, .sigma, .alpha, .Ex, .En, .He, .CD
%        .xCloud (Nx1), .muCloud (Nx1)
%        .xGrid, .muGrid_mean (基于En的解析高斯)，.muGrid_mc (Monte Carlo均值)
%   orderIdx: 分量按目标维均值升序后的索引（K×1）
%
% 公式：
%   alpha_左 = (μ_k - μ_{k-1}) / (3*(σ_{k-1}+σ_k))
%   alpha_右 = (μ_{k+1} - μ_k) / (3*(σ_k+σ_{k+1}))
%   alpha_k  = min([1, alpha_左(若有), alpha_右(若有)]) 裁剪到 [0,1]
%   Ex_k = μ_k
%   En_k = (1 + alpha_k) * σ_k / 2
%   He_k = (1 - alpha_k) * σ_k / 6
%   CD_k = 3*He_k/En_k = (1 - alpha_k)/(1 + alpha_k)
%
% 示例用法：
%   % 已有 gm (fitgmdist)：
%   [CMs, idx] = gmm_to_clouds_FCG(gm, 1, 'NPerComp', 4000, 'ShowPlot', true);

% -------------------- 解析参数 --------------------
p = inputParser;
p.addParameter('NPerComp', 3000, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('ProportionalToWeight', true, @(x)islogical(x));
p.addParameter('AlphaClip', true, @(x)islogical(x));
p.addParameter('RngSeed', [], @(x)isempty(x)||isscalar(x));
p.addParameter('ShowPlot', true, @(x)islogical(x));
p.addParameter('GridX', [], @(x)isnumeric(x)&&isvector(x));
p.parse(varargin{:});

NPerComp   = p.Results.NPerComp;
propW      = p.Results.ProportionalToWeight;
alphaClip  = p.Results.AlphaClip;
rngSeed    = p.Results.RngSeed;
doPlot     = p.Results.ShowPlot;
xGridUser  = p.Results.GridX;

if ~isempty(rngSeed), rng(rngSeed); end

% -------------------- 提取均值与方差(目标维)并排序 --------------------
K = gm.NumComponents;
mu_all = gm.mu(:, dim);                 % K×1
sigma_all = zeros(K,1);                 % 目标维的标准差

for k = 1:K
    sigma_all(k) = get_sigma_dim(gm, k, dim);
end

% 按目标维均值升序排序，记录重排索引
[mu_sorted, orderIdx] = sort(mu_all, 'ascend');
sigma_sorted = sigma_all(orderIdx);
w_sorted = gm.ComponentProportion(orderIdx);

% -------------------- 计算每个分量的 alpha --------------------
alpha = ones(K,1);
for k = 1:K
    candidates = [];
    if k > 1
        numL = mu_sorted(k) - mu_sorted(k-1);
        denL = 3*(sigma_sorted(k-1) + sigma_sorted(k));
        if denL > 0, candidates(end+1) = numL/denL; end
    end
    if k < K
        numR = mu_sorted(k+1) - mu_sorted(k);
        denR = 3*(sigma_sorted(k) + sigma_sorted(k+1));
        if denR > 0, candidates(end+1) = numR/denR; end
    end
    if ~isempty(candidates)
        alpha(k) = min([1, candidates]); % 不超过1
    else
        alpha(k) = 1; % 单分量或无邻居时给1（无边界重叠）
    end
end
if alphaClip
    alpha = max(0, min(1, alpha));      % 裁剪到[0,1]
end

% -------------------- 映射到云模型参数 Ex/En/He 与混淆度 CD --------------------
Ex = mu_sorted;
En = (1 + alpha) .* sigma_sorted / 2;
He = (1 - alpha) .* sigma_sorted / 6;
CD = 3*He ./ En; % = (1 - alpha)./(1 + alpha)

% -------------------- 生成云滴（前向云生成 FCG） --------------------
% 云滴数量：等量或按权重分配
if propW
    N_all = NPerComp * K;
    Nk = max(1, round(N_all * w_sorted / sum(w_sorted)));
else
    Nk = repmat(NPerComp, K, 1);
end

% 自动网格用于绘制解析/MC曲线
if isempty(xGridUser)
    xmin = min(Ex - 5*sigma_sorted);
    xmax = max(Ex + 5*sigma_sorted);
    xGrid = linspace(xmin, xmax, 801).';
else
    xGrid = xGridUser(:);
end

CMs = repmat(struct('k',[], 'weight',[], 'mu',[], 'sigma',[], ...
    'alpha',[], 'Ex',[], 'En',[], 'He',[], 'CD',[], ...
    'xCloud',[], 'muCloud',[], 'xGrid',[], 'muGrid_mean',[], 'muGrid_mc',[]), K, 1);

for k = 1:K
    % FCG: 先采 En' ~ N(En, He^2)，再 x ~ N(Ex, (En')^2)，μ = exp(-(x-Ex)^2/(2(En')^2))
%     [xCloud, muCloud] = forward_cloud_generate(Ex(k), En(k), He(k), Nk(k));
    [xCloud, muCloud] = forward_cloud_generate(Ex(k), En(k), He(k), Nk(k), alpha(k), sigma_sorted(k));
    
    % 解析曲线（基于 En 的“理想”高斯）
    muGrid_mean = exp(- (xGrid - Ex(k)).^2 ./ (2*(En(k).^2)));
    
    % Monte Carlo 平均曲线：近似 E_{En'~N(En,He^2)} [ exp(-(x-Ex)^2/(2En'^2)) ]
    % 这里进行修改 2025.1.7
%     % 这里用较小样本近似，提高速度
%     [~, muMC] = forward_cloud_generate(Ex(k), En(k), He(k), min(4000, max(1000, round(Nk(k)/3))));
%     % 为了得到曲线，需要把 xCloud 投到栅格上做局部平均（简化：用 ksdensity 做权重平均）
%     muGrid_mc = grid_mc_curve(xCloud, muCloud, xGrid);

    % Monte Carlo 平均曲线：用独立采样的 (xMC, muMC) 来估计 （修改的地方 2025.1.7）
    [xMC, muMC] = forward_cloud_generate(Ex(k), En(k), He(k), ...
    min(4000, max(1000, round(Nk(k)/3))), alpha(k), sigma_sorted(k));
    muGrid_mc = grid_mc_curve(xMC, muMC, xGrid);


    
    CMs(k).k       = k;
    CMs(k).weight  = w_sorted(k);
    CMs(k).mu      = mu_sorted(k);
    CMs(k).sigma   = sigma_sorted(k);
    CMs(k).alpha   = alpha(k);
    CMs(k).Ex      = Ex(k);
    CMs(k).En      = En(k);
    CMs(k).He      = He(k);
    CMs(k).CD      = CD(k);
    CMs(k).xCloud  = xCloud;
    CMs(k).muCloud = muCloud;
    CMs(k).xGrid   = xGrid;
    CMs(k).muGrid_mean = muGrid_mean;
    CMs(k).muGrid_mc   = muGrid_mc;
end

% -------------------- 可视化 --------------------
if doPlot
    figure('Color','w'); 
    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
    
    % (1) 云滴散点 + 解析/MC 曲线
    nexttile;
    hold on; grid on; box on;
    cmap = lines(K);
    legtxt = cell(K,1);
    for k = 1:K
        scatter(CMs(k).xCloud, CMs(k).muCloud, 6, 'MarkerEdgeAlpha',0.25, ...
            'MarkerEdgeColor', cmap(k,:), 'MarkerFaceAlpha', 0.08, 'MarkerFaceColor', cmap(k,:));
        plot(CMs(k).xGrid, CMs(k).muGrid_mean, '-', 'LineWidth', 1.6, 'Color', cmap(k,:));
        plot(CMs(k).xGrid, CMs(k).muGrid_mc, '--', 'LineWidth', 1.0, 'Color', cmap(k,:));
        legtxt{k} = sprintf('Comp %d: Ex=%.3f, En=%.3f, He=%.3f, CD=%.3f', ...
            k, CMs(k).Ex, CMs(k).En, CMs(k).He, CMs(k).CD);
    end
    xlabel(sprintf('Dimension %d (Angle)', dim));
    ylabel('Membership \mu');
    title('Forward Cloud Generation (Cloud drops + mean/MC curves)');
    legend(legtxt, 'Location','bestoutside');
    
    % (2) 邻接边界示意：μ±3ασ
    nexttile;
    hold on; grid on; box on;
    for k = 1:K
        xL = mu_sorted(k) - 3*alpha(k)*sigma_sorted(k);
        xR = mu_sorted(k) + 3*alpha(k)*sigma_sorted(k);
        plot([xL xR], [k k], '-', 'Color', cmap(k,:), 'LineWidth', 4);
        plot(mu_sorted(k), k, 'o', 'MarkerFaceColor', cmap(k,:), 'MarkerEdgeColor','k');
    end
    xlabel(sprintf('Dimension %d (Angle)', dim));
    yticks(1:K); yticklabels(arrayfun(@(i)sprintf('Comp %d',i), 1:K, 'UniformOutput', false));
    title('Overlap control via \alpha: [\mu-3\alpha\sigma, \mu+3\alpha\sigma]');
end

end % <-- 主函数结束

% ==================== 工具函数 ====================

function s = get_sigma_dim(gm, k, dim)
% 取得第 k 个分量在第 dim 维上的标准差
S = gm.Sigma;
if strcmpi(gm.CovarianceType, 'full')
    if gm.SharedCovariance
        s = sqrt(S(dim, dim));
    else
        s = sqrt(S(dim, dim, k));
    end
else % 'diagonal'
    if gm.SharedCovariance
        Sd = S(:); % Dx1
        s = sqrt(Sd(dim));
    else
        s = sqrt(S(dim, 1, k)); % Dx1xK
    end
end
end

% function [xCloud, muCloud] = forward_cloud_generate(Ex, En, He, N)
% % 前向云生成（FCG）
% Enp = En + He .* randn(N,1);
% Enp(Enp <= eps) = En; % 保证正数（极少数采样可能<=0）
% xCloud = Ex + Enp .* randn(N,1);
% muCloud = exp(- (xCloud - Ex).^2 ./ (2*(Enp.^2)));
% end

function [xCloud, muCloud] = forward_cloud_generate(Ex, En, He, N, alpha, sigma)
% 前向云生成（FCG）+ 区间裁剪，使 En' 落在 [alpha*sigma, sigma]，避免过度分散

% 1) 采样 En'
Enp = En + He .* randn(N,1);

% 2) 将 En' 限制在论文“受控粒度范围”内
lo = max(eps, alpha * sigma);
hi = max(lo, sigma);          % 防御：确保 hi>=lo
Enp = min(max(Enp, lo), hi);

% 3) 生成云滴与隶属度
xCloud = Ex + Enp .* randn(N,1);
muCloud = exp(- (xCloud - Ex).^2 ./ (2*(Enp.^2)));
end


function muGrid = grid_mc_curve(x, mu, xGrid)
% 将离散云滴 (x, mu) 投影到栅格，做局部加权平均近似 MC 曲线
% 这里用简单的箱平均法（速度快、依赖少）
nbin = numel(xGrid);
muGrid = zeros(nbin,1);
% 为避免空箱，构造栅格边界
edges = [-inf; (xGrid(1:end-1)+xGrid(2:end))/2; inf];
[~,~,bin] = histcounts(x, edges);
for i = 1:nbin
    idx = (bin == i);
    if any(idx)
        muGrid(i) = mean(mu(idx));
    else
        % 如果该格没有点，使用最近邻插值
        [~, j] = min(abs(x - xGrid(i)));
        muGrid(i) = mu(j);
    end
end
% 平滑一下（可选）
muGrid = movmean(muGrid, 5);
end
