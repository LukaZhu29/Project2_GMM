function out = Classic_RULA(AnglesDeg)
% RULA_SCORE_CLASSIC
%  依据 RULA 工作表（去掉所有 adjustment factors），用 A/B/C 三表逐级查表，
%  计算离散的 RULA 分数（非模糊版）。
%
% 输入:
%   AnglesDeg : N×6 矩阵，列依次为
%               [Neck, Trunk, Leg, UpperArm, LowerArm, Wrist] (单位：度)
%               Wrist 为屈伸角（正负均可）
%
% 输出:
%   out.joint_scores : N×6，六个关节的离散分数（按你给的阈值规则）
%   out.A, out.B     : N×1，Section-A / Section-B 的离散分
%   out.Grand        : N×1，最终 RULA 等级（1..7）
%   out.info         : 说明与域定义
%
% 约束 & 简化（与需求一致）：
%   - 忽略所有 adjustment factors（外展/扭转/耸肩/负荷/肌力等）；
%   - Wrist twist 视为中性（=1），不加额外列变化；
%   - LowerArm “>120° = 2”；
%   - Wrist 两档：|角度|≤15 => 1；>15 => 2（映射到官方 1..4 档中的 1/2）。
%
% 参考：与 fuzzy 版本同样使用官方 A/B/C 查表（内置），但本函数不做模糊化，
%       直接“硬分段→查表”得到 A/B/C。
%
NECK=1; TRUNK=2; LEG=3; UARM=4; LARM=5; WRIST=6;

X = AnglesDeg;
N = size(X,1);

% ------------ 0) 官方 A/B/C 表（去掉调整项，固定 twist=1 的切片） ------------
[T_A, T_B, T_C, dom] = get_rula_tables_official();
% 我们仅使用：
%   A: UpperArm(1..6) × LowerArm(1..3) × Wrist(1..4) × Twist(=1) → 1..9
%   B: Neck(1..6) × Trunk(1..6) × Leg(1..2) → 1..9
%   C: A(1..9) × B(1..9) → Grand(1..7)（>8、>7 已在表内扩展处理）

% ------------ 1) 六关节硬阈值打分（完全按你给的规则） ------------
joint = zeros(N,6);

% Neck：负角=4；0–10=1；10–20=2；>20=3
for n=1:N
    x = X(n,NECK);
    if x < 0
        s = 4;
    elseif x <= 10
        s = 1;
    elseif x <= 20
        s = 2;
    else
        s = 3;
    end
    % B 表需要 Neck 1..6，这里只会到 4，不用 5/6
    joint(n,NECK) = s;
end

% Trunk：0=1；1–20=2；20–60=3；>60=4（不考虑侧弯/扭转+1）
for n=1:N
    x = X(n,TRUNK);
    if x == 0
        s = 1;
    elseif x <= 20
        s = 2;
    elseif x <= 60
        s = 3;
    else
        s = 4;
    end
    joint(n,TRUNK) = s;
end

% Leg：0–30=1；30–60=2（>60 也视作 2）
for n=1:N
    x = X(n,LEG);
    if x <= 30
        s = 1;
    else
        s = 2;
    end
    joint(n,LEG) = s;
end

% UpperArm：<-20→2；[-20,20]→1；(20,45]→2；(45,90]→3；≥90→4
% A 表允许到 6；我们按需求只产生 1..4
for n=1:N
    x = X(n,UARM);
    if x < -20
        s = 2;
    elseif x <= 20
        s = 1;
    elseif x <= 45
        s = 2;
    elseif x <= 90
        s = 3;
    else
        s = 4;
    end
    joint(n,UARM) = s;
end

% LowerArm：0–60=2；60–100=1；100–120=2；>120=2
% A 表允许到 3；这里只产生 1 或 2
for n=1:N
    x = X(n,LARM);
    if x <= 60
        s = 2;
    elseif x <= 100
        s = 1;
    elseif x <= 120
        s = 2;
    else
        s = 2;
    end
    joint(n,LARM) = s;
end

% Wrist（屈伸）：|角度|≤15→1；>15→2
% A 表 wrist 原有 1..4 档；此处仅 1/2，3/4 不使用
for n=1:N
    w = abs(X(n,WRIST));
    if w <= 15
        s = 1;
    else
        s = 2;
    end
    joint(n,WRIST) = s;
end

% 关节分数输出
out.joint_scores = joint;
out.info.joint_order = {'Neck','Trunk','Leg','UpperArm','LowerArm','Wrist'};

% ------------ 2) A 表查表： (UpperArm, LowerArm, Wrist, Twist=1) → A ------------
U = clamp_to(joint(:,UARM), 1, 6);
L = clamp_to(joint(:,LARM), 1, 3);
W = clamp_to(joint(:,WRIST),1, 4);
twist = 1; % 固定中性

A = zeros(N,1);
for n=1:N
    A(n,1) = T_A(U(n), L(n), W(n), twist);
end
out.A = A;

% ------------ 3) B 表查表： (Neck, Trunk, Leg) → B ------------
Nneck = clamp_to(joint(:,NECK), 1, 6);
Ntrunk= clamp_to(joint(:,TRUNK),1, 6);
Nleg  = clamp_to(joint(:,LEG),  1, 2);

B = zeros(N,1);
for n=1:N
    B(n,1) = T_B(Nneck(n), Ntrunk(n), Nleg(n));
end
out.B = B;

% ------------ 4) C 表查表： (A, B) → Grand（1..7） ------------
% C 表内部已经把 A>8 合并到 8 行、B>7 合并到 7 列
Aq = clamp_to(A, 1, 9);
Bq = clamp_to(B, 1, 9);

Grand = zeros(N,1);
for n=1:N
    Grand(n,1) = T_C(Aq(n), Bq(n));
end
out.Grand = Grand;

% ------------ 附：信息 ------------
out.info.tables = 'official(A/B/C), twist fixed=1, no adjustment factors';
out.info.dom = dom;

% ==================== 内部函数 ====================
function y = clamp_to(x, lo, hi)
    y = min(max(x, lo), hi);
end

end

% ==================== 附：A/B/C 表（去调整项；内置） ====================
function [T_A, T_B, T_C, dom] = get_rula_tables_official()
    % 离散域
    dom.Neck     = 1:6;
    dom.Trunk    = 1:6;
    dom.Leg      = 1:2;
    dom.UpperArm = 1:6;
    dom.LowerArm = 1:3;
    dom.Wrist    = 1:4;
    dom.A        = 1:9;
    dom.B        = 1:9;
    dom.Grand    = 1:7;

    % --- A 表：U×L×W×Twist(1/2) → A(1..9)
    % 行向量顺序：W1_T1,W1_T2,W2_T1,W2_T2,W3_T1,W3_T2,W4_T1,W4_T2
    A_rows = {};
    % U=1
    A_rows{end+1} = [2 2 2 2 3 3 3 3]; % L=1
    A_rows{end+1} = [2 2 2 2 2 3 3 3]; % L=2
    A_rows{end+1} = [2 3 3 3 3 3 4 4]; % L=3
    % U=2
    A_rows{end+1} = [2 3 3 3 3 4 4 4];
    A_rows{end+1} = [3 3 3 3 3 4 4 4];
    A_rows{end+1} = [3 4 4 4 4 4 5 5];
    % U=3
    A_rows{end+1} = [3 3 4 4 4 4 5 5];
    A_rows{end+1} = [3 4 4 4 4 4 5 5];
    A_rows{end+1} = [4 4 4 4 4 5 5 5];
    % U=4
    A_rows{end+1} = [4 4 4 4 4 5 5 5];
    A_rows{end+1} = [4 4 4 4 4 5 5 5];
    A_rows{end+1} = [4 4 4 5 5 5 6 6];
    % U=5
    A_rows{end+1} = [5 5 5 5 5 6 6 7];
    A_rows{end+1} = [5 6 6 6 6 7 7 7];
    A_rows{end+1} = [6 6 6 7 7 7 7 8];
    % U=6
    A_rows{end+1} = [7 7 7 7 7 8 8 9];
    A_rows{end+1} = [8 8 8 8 8 9 9 9];
    A_rows{end+1} = [9 9 9 9 9 9 9 9];

    U=6; L=3; W=4; Tw=2;
    T_A = zeros(U,L,W,Tw);
    for iu=1:U
        for il=1:L
            row = A_rows{(iu-1)*L + il};
            for iw=1:W
                T_A(iu,il,iw,1) = row(2*iw-1); % Twist=1（中性）
                T_A(iu,il,iw,2) = row(2*iw);   % Twist=2（未使用）
            end
        end
    end

    % --- B 表：Neck×Trunk×Leg → B(1..9)
    % 每个 Trunk 有 Leg=1/2 两列，顺序：(T1_L1, T1_L2, T2_L1, T2_L2, ...)
    B_mat = [...
        1 3 2 3 3 4 5 5 6 6 7 7;   % N=1
        2 3 2 3 4 5 5 5 6 7 7 7;   % N=2
        3 3 3 4 4 5 5 6 6 7 7 7;   % N=3
        5 5 5 6 6 7 7 7 7 7 8 8;   % N=4
        7 7 7 7 7 8 8 8 8 8 8 8;   % N=5
        8 8 8 8 8 8 8 9 9 9 9 9];  % N=6
    Nn=6; Nt=6; Nl=2;
    T_B = zeros(Nn,Nt,Nl);
    for in=1:Nn
        row = B_mat(in,:);
        for it=1:Nt
            T_B(in,it,1) = row(2*it-1);
            T_B(in,it,2) = row(2*it);
        end
    end

    % --- C 表：A×B → Grand(1..7)
    % 行 A=1..7,8+；列 B=1..7+；扩展到 9×9 矩阵，>8 行、>7 列按边界并入
    C_mat = [...
        1 2 3 3 4 5 5;  % A=1
        2 2 3 4 4 5 5;  % A=2
        3 3 3 4 4 5 6;  % A=3
        3 3 3 4 5 6 6;  % A=4
        4 4 4 5 6 7 7;  % A=5
        4 4 5 6 6 7 7;  % A=6
        5 5 6 6 7 7 7;  % A=7
        5 5 6 7 7 7 7]; % A=8+
    T_C = zeros(9,9);
    T_C(1:8,1:7) = C_mat;      % 基本区
    T_C(9,1:7)   = C_mat(8,:); % A=9 → 8+
    T_C(:,8)     = T_C(:,7);   % B=8 → 7+
    T_C(:,9)     = T_C(:,7);   % B=9 → 7+
end
