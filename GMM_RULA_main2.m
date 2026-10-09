function out = GMM_RULA_main2(AnglesDeg, CMs, varargin)
% ...（前略注释同前版本）...

    p = inputParser;
    p.addParameter('LowerArmRiskVals', [2, 1, 2], @(x)isnumeric(x)&&numel(x)==3);
    p.addParameter('UseCoAGrid', false, @(x)islogical(x)&&isscalar(x));
    p.parse(varargin{:});
    sFA = p.Results.LowerArmRiskVals(:).';  % 1x3

    % >>> 关键改动：将CMs标准化为列向量 <<<
    CMs = CMs(:);
    assert(numel(CMs)==3, 'CMs 应该包含 3 个云模型（前臂三档）。');

    AnglesDeg = double(AnglesDeg);
    assert(size(AnglesDeg,2)==6, 'AnglesDeg must be N×6: [Neck Trunk Leg UpperArm LowerArm Wrist].');

    N = size(AnglesDeg,1);
    Neck      = AnglesDeg(:,1);
    Trunk     = AnglesDeg(:,2);
    Leg       = AnglesDeg(:,3);
    UpperArm  = AnglesDeg(:,4);
    LowerArm  = AnglesDeg(:,5);
    Wrist     = AnglesDeg(:,6);

    % 1) 非前臂
    sNeck     = neckScore_basic(Neck);
    sTrunk    = trunkScore_basic(Trunk);
    sLeg      = legScore_basic(Leg);
    sUpper    = upperArmScore_basic(UpperArm);
    sWrist    = wristScore_basic(Wrist);

    % 2) 前臂（模糊）
    [muFA, sFA_cont] = fuzzy_lowerarm_from_clouds(LowerArm, CMs, sFA);


    % 3) Table-A 前臂维模糊插值
    A_fuzzy = rula_tableA_fuzzy_interp(sUpper, sWrist, muFA);

    % 4) 后续表
    A_rounded = max(1, min(9, round(A_fuzzy)));
    [B, C, Grand] = rula_tables_basic(A_rounded, sNeck, sTrunk, sLeg);

    out = struct();
    out.A_fuzzy        = A_fuzzy;
    out.B              = B;
    out.C              = C;
    out.Grand          = Grand;
    out.LowerArm_cont  = sFA_cont;
    out.LowerArm_mu    = muFA;       % N×3
    out.UpperArm       = sUpper;
    out.Wrist          = sWrist;
    out.Neck           = sNeck;
    out.Trunk          = sTrunk;
    out.Leg            = sLeg;
end

function s = neckScore_basic(theta)
% RULA-颈部（无调整）：常用简化
% 0~10 ->1 ; 10~20 ->2 ; >20 ->3 ; (<0 轻微后伸) ->1
    t = abs(theta);
    s = ones(size(theta));
    s(t>10 & t<=20) = 2;
    s(t>20) = 3;
end

function s = trunkScore_basic(theta)
% RULA-躯干（无调整）：简化
% 0~20 ->1 ; 20~60 ->2 ; >60 ->3 ; 轻度后伸 ->1
    t = abs(theta);
    s = ones(size(theta));
    s(t>20 & t<=60) = 2;
    s(t>60) = 3;
end

function s = legScore_basic(thetaLike)
% RULA-腿：简化为是否“平衡/支撑”
% 这里以一个“角度差/稳定度”代理：|θ|<=5 ->1(平衡) ; 否则 2(不平衡)
    s = ones(size(thetaLike));
    s(abs(thetaLike)>5) = 2;
end

function s = upperArmScore_basic(theta)
% RULA-上臂（无调整）简化（常见阈值）
% <=20° ->1 ; 20~45 ->2 ; >45 ->3 ; （轻度后伸亦按 2/3 类似处理）
    s = ones(size(theta));
    s(theta>20 & theta<=45) = 2;
    s(theta>45) = 3;
    s(theta<-20 & theta>=-45) = 2;
    s(theta<-45) = 3;
end

function s = wristScore_basic(theta)
% RULA-腕（仅屈伸，无偏斜/扭转调整）
% 0~15 ->1 ; 15~30 ->2 ; >30 ->3 ; （若要兼容 1~4 档，可把极端>45 设为4）
    t = abs(theta);
    s = ones(size(theta));
    s(t>15 & t<=30) = 2;
    s(t>30) = 3;
end
function val = rula_tableA_lookup(U, L, W)
% Table-A(U,L,W)：U=1..6, L=1..3, W=1..4
% 下面矩阵来源于常见 RULA 实作的无调整项版本（只按位置分），
% 若你已有官方表，请替换为你的版本。
%
% 实现方式：给出 W=1..4 四个 6x3 矩阵；val = T{W}(U,L)

    persistent T;
    if isempty(T)
        % ---- Wrist=1 时的 6x3 (U×L) ----
        T1 = [ 1 2 2;
               2 3 3;
               3 4 4;
               4 5 5;
               5 6 6;
               6 7 7 ];
        % ---- Wrist=2 ----
        T2 = [ 2 2 3;
               3 3 4;
               4 4 5;
               5 5 6;
               6 6 7;
               7 7 8 ];
        % ---- Wrist=3 ----
        T3 = [ 3 3 4;
               4 4 5;
               5 5 6;
               6 6 7;
               7 7 8;
               8 8 9 ];
        % ---- Wrist=4 ----
        T4 = [ 3 4 4;
               4 5 5;
               5 6 6;
               6 7 7;
               7 8 8;
               8 9 9 ];

        T = {T1,T2,T3,T4};
    end

    U = max(1,min(6,round(U)));
    L = max(1,min(3,round(L)));
    W = max(1,min(4,round(W)));

    val = T{W}(U,L);
end

function [B, C, Grand] = rula_tables_basic(A, sNeck, sTrunk, sLeg)
% Table-B：由颈/躯干/腿得到 B（简化版的常用表）
% Table-C：由 A 与 B 得到 C，再作为 Grand（此简化版不加“肌肉/负载/活动”）
% 若你已有权威表，请直接替换。
    N = numel(A);
    B  = zeros(N,1);
    C  = zeros(N,1);
    Grand = zeros(N,1);

    for i = 1:N
        % ---- Table-B（简化 3D 表：N(1..3) × T(1..3) × L(1..2) -> 1..7）----
        Nn = max(1,min(3,round(sNeck(i))));
        Tn = max(1,min(3,round(sTrunk(i))));
        Ln = max(1,min(2,round(sLeg(i))));
        B(i) = rula_tableB_lookup(Nn, Tn, Ln);

        % ---- Table-C（A × B -> 1..9）----
        Ai = max(1,min(9,round(A(i))));
        Bi = max(1,min(7,round(B(i))));
        C(i) = rula_tableC_lookup(Ai, Bi);

        Grand(i) = C(i); % 无任何后续调整时，Grand=C
    end
end

function val = rula_tableB_lookup(Nn, Tn, Ln)
% 一个常见的简化 Table-B（不含调整），3×3×2 -> 1..7
% 这里给出一个广泛使用的近似表；如你有官方矩阵请直接替换。
    persistent B1 B2
    if isempty(B1)
        % L=1（腿平衡）
        B1 = [ 1 2 3;
               2 3 4;
               3 4 5 ];
        % L=2（腿不平衡）
        B2 = [ 2 3 4;
               3 4 5;
               4 5 6 ];
    end
    Nn = max(1,min(3,Nn));
    Tn = max(1,min(3,Tn));
    Ln = max(1,min(2,Ln));
    if Ln==1, val = B1(Nn,Tn); else, val = B2(Nn,Tn); end
end

function val = rula_tableC_lookup(Ai, Bi)
% 一个常见的简化 Table-C（A×B -> 1..9），无调整项
% 给出趋势型矩阵：风险随 A/B 单调上升
    persistent Cmat
    if isempty(Cmat)
        % 行：A=1..9；列：B=1..7
        Cmat = [ 1 2 2 3 3 4 4;
                 2 2 3 3 4 4 5;
                 2 3 3 4 4 5 5;
                 3 3 4 4 5 5 6;
                 3 4 4 5 5 6 6;
                 4 4 5 5 6 6 7;
                 4 5 5 6 6 7 8;
                 5 5 6 6 7 8 8;
                 5 6 6 7 8 8 9 ];
    end
    Ai = max(1,min(9,Ai));
    Bi = max(1,min(7,Bi));
    val = Cmat(Ai, Bi);
end
