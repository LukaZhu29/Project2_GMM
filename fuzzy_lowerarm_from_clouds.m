function [muFA, sFA_cont] = fuzzy_lowerarm_from_clouds(thetaFA, CMs, sFA_singletons)
% 输入：thetaFA (N×1)，CMs 为 3×1 或 1×3 struct，字段含 .Ex .En
%      sFA_singletons(1×3) 为每个 CM 的后件单点分数（如 [2 1 2]）

    % >>> 关键改动：将CMs标准化为列向量 <<<
    CMs = CMs(:);

    N = numel(thetaFA);
    K = numel(CMs);
    assert(K==3, '此实现假定 lower arm 只有 3 个 CM（与 RULA 前臂三档一致）。');

    muFA = zeros(N, K);
    for k = 1:K
        Ex = CMs(k).Ex;
        En = CMs(k).En;
        muFA(:,k) = exp(- (thetaFA(:)-Ex).^2 ./ (2*(En.^2)));  % 解析均值曲线
    end

    w = sum(muFA, 2) + eps;
    sFA_cont = (muFA * sFA_singletons(:)) ./ w;   % CoA：单点后件的加权平均
end
