function res = rebaSensitivityAnalysis4(X, opts)
%{
rebaSensitivityAnalysis2： 增添了局部放大功能 
rebaSensitivityAnalysis3： 修复了放大图挡原图的问题
%}
% rebaSensitivityAnalysis
% Input:
%   X    : n×4 matrix [theta, gmmScore, fuzzyScore, rebaScore]
%          theta expected to be 0:1:120 (degree)
%   opts : (optional) struct with fields:
%          .epsilon_mode  = 'fixed' (default) | 'quantile'
%          .epsilon       = 0.01 (default, score/deg) when epsilon_mode='fixed'
%          .q_ref         = 95   (default) percentile used when epsilon_mode='quantile'
%          .q_frac        = 0.01 (default) fraction multiplied to percentile when epsilon_mode='quantile'
%          .use_semilogy  = true (default) use log-y for derivative plot
%          .save_figs     = false (default)
%          .save_dir      = pwd (default)
%          .fig_prefix    = 'reba_sens' (default)
%
% Output:
%   res struct with fields:
%     .theta, .scores
%     .theta_mid, .dAbs
%     .epsilon
%     .activeMask
%     .metricsTable  (Coverage, ActiveMean, AUC)
%
% Notes:
% - Derivative is computed as forward difference magnitude:
%     S(i) = |R(i+1)-R(i)| / dtheta
%   on midpoint theta_mid = theta(1:end-1) + dtheta/2
% - Coverage is computed over derivative grid length (n-1).

    if nargin < 2 || isempty(opts), opts = struct(); end
    opts = fillOpts(opts);

    % ---------- Step 1: validate/sort ----------
    assert(size(X,2) == 4, 'X must be n×4: [theta, gmm, fuzzy, reba].');
    X = sortrows(X, 1);
    theta = X(:,1);
    Rgmm  = X(:,2);
    Rfuz  = X(:,3);
    Rreba = X(:,4);

    % Check uniform grid
    dtheta = median(diff(theta));
    if any(abs(diff(theta) - dtheta) > 1e-9)
        warning('Angle grid is not perfectly uniform; using median dtheta=%.6g.', dtheta);
    end
    if ~(abs(dtheta-1) < 1e-9)
        warning('Expected dtheta=1 deg; current dtheta=%.6g. Metrics will use this dtheta.', dtheta);
    end

    % ---------- Step 2: numeric derivative (abs) ----------
    theta_mid = theta(1:end-1) + dtheta/2;
    Sgmm  = abs(diff(Rgmm))  ./ dtheta;
    Sfuz  = abs(diff(Rfuz))  ./ dtheta;
    Sreba = abs(diff(Rreba)) ./ dtheta;

    % ---------- Step 3: epsilon selection ----------
    if strcmpi(opts.epsilon_mode, 'fixed')
        epsilon = opts.epsilon;
    elseif strcmpi(opts.epsilon_mode, 'quantile')
        Sall = [Sgmm(:); Sfuz(:); Sreba(:)];
        qv = prctile(Sall, opts.q_ref);
        epsilon = opts.q_frac * qv;
    else
        error('Unknown epsilon_mode: %s (use fixed|quantile).', opts.epsilon_mode);
    end

    active_gmm  = (Sgmm  > epsilon);
    active_fuz  = (Sfuz  > epsilon);
    active_reba = (Sreba > epsilon);

    % ---------- Step 4: metrics ----------
    metrics = struct();
    metrics.GMM  = computeMetrics(Sgmm,  active_gmm,  dtheta);
    metrics.Fuzzy= computeMetrics(Sfuz,  active_fuz,  dtheta);
    metrics.REBA = computeMetrics(Sreba, active_reba, dtheta);

    methods = {'GMM-REBA'; 'Fuzzy REBA'; 'REBA'};
    Coverage    = [metrics.GMM.Coverage; metrics.Fuzzy.Coverage; metrics.REBA.Coverage];
    ActiveMean  = [metrics.GMM.ActiveMean; metrics.Fuzzy.ActiveMean; metrics.REBA.ActiveMean];
    AUC         = [metrics.GMM.AUC; metrics.Fuzzy.AUC; metrics.REBA.AUC];
    metricsTable = table(Coverage, ActiveMean, AUC, 'RowNames', methods);

    % ---------- Step 5: plots ----------
    plotDerivative(theta_mid, Sgmm, Sfuz, Sreba, epsilon, opts);
    plotMetricsBars(metricsTable, opts);

    % ---------- Step 6: pack output ----------
    res = struct();
    res.theta = theta;
    res.scores = struct('GMM',Rgmm,'Fuzzy',Rfuz,'REBA',Rreba);
    res.theta_mid = theta_mid;
    res.dAbs = struct('GMM',Sgmm,'Fuzzy',Sfuz,'REBA',Sreba);
    res.epsilon = epsilon;
    res.activeMask = struct('GMM',active_gmm,'Fuzzy',active_fuz,'REBA',active_reba);
    res.metricsTable = metricsTable;

    % Print a compact summary
    fprintf('\n=== Sensitivity metrics (epsilon = %.6g, mode=%s) ===\n', epsilon, opts.epsilon_mode);
    disp(metricsTable);
end

% ========================= Helpers =========================

function opts = fillOpts(opts)
    if ~isfield(opts,'epsilon_mode'), opts.epsilon_mode = 'fixed'; end
    if ~isfield(opts,'epsilon'),      opts.epsilon = 0.01; end % score/deg
    if ~isfield(opts,'q_ref'),        opts.q_ref = 95; end
    if ~isfield(opts,'q_frac'),       opts.q_frac = 0.01; end
    if ~isfield(opts,'use_semilogy'), opts.use_semilogy = true; end
    if ~isfield(opts,'save_figs'),    opts.save_figs = false; end
    if ~isfield(opts,'save_dir'),     opts.save_dir = pwd; end
    if ~isfield(opts,'fig_prefix'),   opts.fig_prefix = 'reba_sens'; end
    if ~isfield(opts,'zoom_inset'),    opts.zoom_inset = true; end
    if ~isfield(opts,'zoom_xwins'),    opts.zoom_xwins = [55 65; 95 105]; end % two local windows
    if ~isfield(opts,'zoom_ymax_main'),opts.zoom_ymax_main = []; end % auto if empty
    if ~isfield(opts,'zoom_use_log_inset'), opts.zoom_use_log_inset = true; end
    if ~isfield(opts,'inset_pos'),          opts.inset_pos = [0.12 0.58 0.30 0.33]; end
    % 建议默认放左上（通常不挡 55–105°那段主变化），你也可以改成右下：[0.62 0.18 0.30 0.33]
    if ~isfield(opts,'inset_transparent'),  opts.inset_transparent = true; end
    if ~isfield(opts,'zoom_use_log_panel'),opts.zoom_use_log_panel = true; end
end

function m = computeMetrics(S, activeMask, dtheta)
    % Coverage over derivative grid
    m.Coverage = mean(activeMask);

    % Conditional mean sensitivity over active region
    if any(activeMask)
        m.ActiveMean = mean(S(activeMask));
    else
        m.ActiveMean = 0;
    end

    % AUC / total variation
    m.AUC = sum(S) * dtheta;
end

function plotDerivative(theta_mid, Sgmm, Sfuz, Sreba, epsilon, opts)
% plotDerivative (paper-style layout)
% Left: main panel (full range, capped y)
% Right: zoom panel (threshold regions), no overlap

    % ---------- figure ----------
    figure('Name','|dR/d\theta| comparison','Color','w');

    % ---------- determine y-cap for main panel ----------
    if isempty(opts.zoom_ymax_main)
        sCont = [Sgmm(:); Sfuz(:)];
        if all(sCont == 0)
            ycap = max(10*epsilon, 0.05);
        else
            ycap = prctile(sCont, 99);
            ycap = max(ycap, 10*epsilon);
        end
    else
        ycap = opts.zoom_ymax_main;
    end

    % ====================================================
    % Left panel: main axes (full range)
    % ====================================================
    axMain = axes('Position',[0.08 0.15 0.62 0.75]); %#ok<LAXES>
    hold(axMain,'on');

    plot(axMain, theta_mid, Sgmm,  'LineWidth',1.6);
    plot(axMain, theta_mid, Sfuz,  'LineWidth',1.6);
    plot(axMain, theta_mid, Sreba, 'LineWidth',1.0);

    yline(axMain, epsilon, '--', 'LineWidth',1.2);

    xlabel(axMain,'\theta (deg)');
    ylabel(axMain,'|dR/d\theta| (score/deg)');
    title(axMain,'Score–angle sensitivity (full range)');
    grid(axMain,'on');

    xlim(axMain,[min(theta_mid) max(theta_mid)]);
    ylim(axMain,[0 ycap]);

    lgd = legend(axMain, ...
        {'SEF-REBA','Fuzzy REBA','REBA','\epsilon threshold'}, ...
        'Location','northwest');
    lgd.AutoUpdate = 'off';

    % Highlight zoom regions on main panel
    if ~isempty(opts.zoom_xwins)
%         yl = ylim(axMain);
        yl = [0 0.15];
        for k = 1:size(opts.zoom_xwins,1)
            xw = opts.zoom_xwins(k,:);
            patch(axMain, ...
                [xw(1) xw(2) xw(2) xw(1)], ...
                [yl(1) yl(1) yl(2) yl(2)], ...
                'k', 'FaceAlpha',0.04, ...
                'EdgeColor','none', ...
                'HandleVisibility','off');
        end
    end

    % ====================================================
    % Right panel: zoom axes (threshold regions)
    % ====================================================
    axZoom = axes('Position',[0.74 0.15 0.22 0.75]); %#ok<LAXES>
    hold(axZoom,'on');

    if opts.zoom_use_log_panel
        semilogy(axZoom, theta_mid, Sgmm,  'LineWidth',1.4);
        semilogy(axZoom, theta_mid, Sfuz,  'LineWidth',1.4);
        semilogy(axZoom, theta_mid, Sreba, 'LineWidth',1.0);
        yline(axZoom, max(epsilon, eps), '--', 'LineWidth',1.1, ...
              'HandleVisibility','off');
        ylabel(axZoom,'|dR/d\theta| (score/deg)');
        xlabel(axZoom,'\theta (deg)');
    else
        plot(axZoom, theta_mid, Sgmm,  'LineWidth',1.4);
        plot(axZoom, theta_mid, Sfuz,  'LineWidth',1.4);
        plot(axZoom, theta_mid, Sreba, 'LineWidth',1.0);
        yline(axZoom, epsilon, '--', 'LineWidth',1.1, ...
              'HandleVisibility','off');
        ylabel(axZoom,'|dR/d\theta|');
        xlabel(axZoom,'\theta (deg)');
    end

    % x-limits: union of zoom windows
    xmin = min(opts.zoom_xwins(:,1));
    xmax = max(opts.zoom_xwins(:,2));
    xlim(axZoom,[xmin xmax]);

    grid(axZoom,'on');
    title(axZoom,'Zoom near thresholds');
    set(axZoom,'FontSize',10);

    % y-limits for log scale
    if opts.zoom_use_log_panel
        ymin = max(min([Sgmm(:); Sfuz(:); Sreba(:); epsilon]), eps);
        ymax = max([Sgmm(:); Sfuz(:); Sreba(:); epsilon]);
        if ymax <= ymin, ymax = ymin*10; end
        ylim(axZoom,[ymin ymax]);
    end

    % ---------- save ----------
    if isfield(opts,'save_figs') && opts.save_figs
        ensureDir(opts.save_dir);
        f = fullfile(opts.save_dir, sprintf('%s_derivative_panel.png', opts.fig_prefix));
        exportgraphics(gcf, f, 'Resolution', 300);
    end
end


function plotMetricsBars(metricsTable, opts)
    % 3 bar charts: Coverage, ActiveMean, AUC
    names = metricsTable.Properties.RowNames;

    % Coverage
    figure('Name','Coverage','Color','w');
    bar(metricsTable.Coverage);
    set(gca,'XTickLabel',names,'XTickLabelRotation',15);
    ylabel('Coverage P (fraction of angles with |dR/d\theta| > \epsilon)');
    title('Effective responsive angular range (Coverage)');
    grid on;
    ylim([0 1]);
    if opts.save_figs
        ensureDir(opts.save_dir);
        f = fullfile(opts.save_dir, sprintf('%s_coverage.png', opts.fig_prefix));
        exportgraphics(gcf, f, 'Resolution', 300);
    end

    % ActiveMean
    figure('Name','ActiveMean','Color','w');
    bar(metricsTable.ActiveMean);
    set(gca,'XTickLabel',names,'XTickLabelRotation',15);
    ylabel('Mean(|dR/d\theta| \mid |dR/d\theta| > \epsilon)  (score/deg)');
    title('Conditional mean sensitivity over active regions');
    grid on;
    if opts.save_figs
        ensureDir(opts.save_dir);
        f = fullfile(opts.save_dir, sprintf('%s_activemean.png', opts.fig_prefix));
        exportgraphics(gcf, f, 'Resolution', 300);
    end

    % AUC
    figure('Name','AUC','Color','w');
    bar(metricsTable.AUC);
    set(gca,'XTickLabel',names,'XTickLabelRotation',15);
    ylabel('\int |dR/d\theta| d\theta  (score)');
    title('Total score variation induced by angle changes (AUC / total variation)');
    grid on;
    if opts.save_figs
        ensureDir(opts.save_dir);
        f = fullfile(opts.save_dir, sprintf('%s_auc.png', opts.fig_prefix));
        exportgraphics(gcf, f, 'Resolution', 300);
    end
end

function ensureDir(d)
    if ~exist(d,'dir')
        mkdir(d);
    end
end
