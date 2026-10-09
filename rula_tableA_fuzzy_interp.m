function A_fuzzy = rula_tableA_fuzzy_interp(sUpper, sWrist, muFA)
% sUpper, sWrist: N×1（离散整数）
% muFA:          N×3，分别对应 前臂 L=1,2,3 的触发强度
% 返回：A_fuzzy:  N×1，Table-A 在“前臂维度”的模糊插值结果
%
% 做法：A = (μ1*A(U,1,W) + μ2*A(U,2,W) + μ3*A(U,3,W)) / (μ1+μ2+μ3)

    N = numel(sUpper);
    A_fuzzy = zeros(N,1);
    mu_sum = sum(muFA,2) + eps;

    for i = 1:N
        U = clamp_int(sUpper(i), 1, 6);
        W = clamp_int(sWrist(i), 1, 4);
        A1 = rula_tableA_lookup(U, 1, W);
        A2 = rula_tableA_lookup(U, 2, W);
        A3 = rula_tableA_lookup(U, 3, W);
        A_fuzzy(i) = (muFA(i,1)*A1 + muFA(i,2)*A2 + muFA(i,3)*A3) / mu_sum(i);
    end
end

function x = clamp_int(x, a, b)
    x = max(a, min(b, round(x)));
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