function plot_cloud_scatter(CMs, varargin)
%PLOT_CLOUD_SCATTER  绘制云模型的云滴散点图（类似示例图）
%
% 必要输入：
%   CMs : 由 gmm_to_clouds_FCG 返回的结构体数组（1×K）
%
% Name-Value 可选参数：
%   'XLabel'   : x 轴标签（默认 'Angle'）
%   'YLabel'   : y 轴标签（默认 'Membership'）
%   'Title'    : 标题（默认 ''）
%   'DotSize'  : 散点大小（默认 14）
%   'Alpha'    : 透明度（0~1，默认 0.6）
%   'Colors'   : K×3 颜色矩阵（默认 lines(K)）
%   'XLim'     : x 轴范围（默认自动）
%   'YLim'     : y 轴范围（默认 [0 1.02]）
%   'Legend'   : 是否显示图例（默认 true）
%   'LegendPrefix' : 图例前缀（默认 'CM_'）
%
% 说明：按 Ex 升序绘制，确保从左到右的云与角度顺序一致。

p = inputParser;
p.addParameter('XLabel', 'Angle', @ischar);
p.addParameter('YLabel', 'Membership', @ischar);
p.addParameter('Title',  '', @ischar);
p.addParameter('DotSize', 14, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('Alpha', 0.6, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('Colors', [], @(x)isnumeric(x)&&size(x,2)==3);
p.addParameter('XLim', [], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('YLim', [0 1.02], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('Legend', true, @(x)islogical(x));
p.addParameter('LegendPrefix', 'CM_', @ischar);
p.parse(varargin{:});

xl       = p.Results.XLabel;
yl       = p.Results.YLabel;
ttl      = p.Results.Title;
msz      = p.Results.DotSize;
alphaVal = p.Results.Alpha;
C        = p.Results.Colors;
xlimUser = p.Results.XLim;
ylimUser = p.Results.YLim;
showLeg  = p.Results.Legend;
pre      = p.Results.LegendPrefix;

% 1) 按 Ex 升序排列
K  = numel(CMs);
Ex = arrayfun(@(s)s.Ex, CMs);
[~, idx] = sort(Ex, 'ascend');
CMs = CMs(idx);

if isempty(C), C = lines(K); end

% 2) 绘图
figure('Color','w'); hold on; box on; grid on;
set(gca, 'GridLineStyle',':', 'GridAlpha',0.3, 'LineWidth',1.0);

legtxt = cell(K,1);
for k = 1:K
    % 仅散点：x 为角度，y 为 membership
    scatter(CMs(k).xCloud, CMs(k).muCloud, msz, ...
        'MarkerFaceColor', C(k,:), 'MarkerEdgeColor', C(k,:), ...
        'MarkerFaceAlpha', alphaVal, 'MarkerEdgeAlpha', alphaVal);
    legtxt{k} = sprintf('%s%d', pre, k);
end

% 轴与外观
xlabel(xl, 'Interpreter','none');
ylabel(yl);
title(ttl, 'Interpreter','none');
set(gca, 'FontName','Times New Roman', 'FontSize',12);

% 轴范围：y 固定到 [0,1.02]，x 可以给定或自动
if ~isempty(xlimUser)
    xlim(xlimUser);
else
    xmin = min(cellfun(@(x)min(x), {CMs.xCloud}));
    xmax = max(cellfun(@(x)max(x), {CMs.xCloud}));
    pad  = 0.05*(xmax - xmin + eps);
    xlim([xmin - pad, xmax + pad]);
end
xlim([0 130]);
ylim(ylimUser);

% 图例
if showLeg
    lg = legend(legtxt, 'Location','northeast');
    lg.Box = 'on';
end
end
