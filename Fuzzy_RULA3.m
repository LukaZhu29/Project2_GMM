function out = Fuzzy_RULA3(AnglesDeg, I)
% FUZZY_RULA_IFTHEN_COA
%  按论文：RULA 的 A/B/C 查表 → IF–THEN 规则；
%  前件 AND=min，规则聚合=max，输出用 CoA 去模糊。
%
% 输入:
%   AnglesDeg : N×6 [Neck, Trunk, Leg, UpperArm, LowerArm, Wrist] (deg)
%   I         : 平滑阈值(°)，默认 5
%
% 输出:
%   out.* : 关节层 soft/hard，A/B/Grand 的 soft/hard
%
if nargin<2 || isempty(I), I = 5; end

NECK=1; TRUNK=2; LEG=3; UARM=4; LARM=5; WRIST=6;

% ------- 0) 读入/内置 RULA A/B/C 三表 -------
[T_A, T_B, T_C, dom] = get_rula_tables_official();

X = AnglesDeg; 
N = size(X,1);

% ------- 1) 关节层模糊化（显式用 struct 容器，避免 mu 被覆盖） -------
MUJ = struct('N',[],'T',[],'L',[],'U',[],'F',[],'W',[]);

% Neck：负角=4；0–10=1；10–20=2；>20=3（兼容 dom.Neck=1..6，5/6 不用）
MUJ.N = zeros(N, numel(dom.Neck));
for n=1:N
    x = X(n,NECK);
    MUJ.N(n,4) = trap_left (x, -Inf, 0, I);
    MUJ.N(n,1) = trap_mid  (x, 0, 10, I);
    MUJ.N(n,2) = trap_mid  (x, 10, 20, I);
    MUJ.N(n,3) = trap_right(x, 20, Inf, I);
end
[softN, hardN] = defuzz_disc_batch(dom.Neck, MUJ.N, [1 2 3 4], true);

% Trunk：0=1；1–20=2；20–60=3；>60=4（兼容 dom.Trunk=1..6，5/6 不用）
MUJ.T = zeros(N, numel(dom.Trunk));
for n=1:N
    x = X(n,TRUNK);
    MUJ.T(n,1) = max(trap_mid(x,0,1,I), trap_mid(x,0,0+I,I));
    MUJ.T(n,2) = trap_mid (x,1,20,I);
    MUJ.T(n,3) = trap_mid (x,20,60,I);
    MUJ.T(n,4) = trap_right(x,60,Inf,I);
end
[softT, hardT] = defuzz_disc_batch(dom.Trunk, MUJ.T, [1 2 3 4], true);

% Leg：0–30=1；30–60=2；>60→2
MUJ.L = zeros(N, numel(dom.Leg));
for n=1:N
    x = X(n,LEG);
    MUJ.L(n,1) = trap_mid (x,0,30,I);
    MUJ.L(n,2) = max(trap_mid(x,30,60,I), trap_right(x,60,Inf,I));
end
[softL, hardL] = defuzz_disc_batch(dom.Leg, MUJ.L, [1 2], true);

% UpperArm：<-20→2；[-20,20]→1；(20,45]→2；(45,90]→3；≥90→4（dom 1..6，5/6 不用）
MUJ.U = zeros(N, numel(dom.UpperArm));
for n=1:N
    x = X(n,UARM);
    part2a    = trap_left(x,-Inf,-20,I);
    MUJ.U(n,1)= trap_mid(x,-20,20,I);
    part2b    = trap_mid(x,20,45,I);
    MUJ.U(n,2)= max(part2a, part2b);
    MUJ.U(n,3)= trap_mid(x,45,90,I);
    MUJ.U(n,4)= trap_right(x,90,Inf,I);
end
[softU, hardU] = defuzz_disc_batch(dom.UpperArm, MUJ.U, [1 2 3 4], true);

% LowerArm（按你要求仅 1/2 两档；3 档恒 0）
MUJ.F = zeros(N, numel(dom.LowerArm)); % 1..3
for n=1:N
    x = X(n,LARM);
    mu2 = max(trap_mid(x,0,60,I), trap_mid(x,100,120,I));
    mu2 = max(mu2, trap_right(x,120,Inf,I));
    mu1 = trap_mid(x,60,100,I);
    MUJ.F(n,1) = mu1;   % score=1
    MUJ.F(n,2) = mu2;   % score=2
    MUJ.F(n,3) = 0;     % 不用
end
[softF, hardF] = defuzz_disc_batch(dom.LowerArm, MUJ.F, [1 2], false);

% Wrist：把 |角度| 映射到 1..4 档（≤15，15–30，30–45，>45）
MUJ.W = zeros(N, numel(dom.Wrist));
for n=1:N
    w = abs(X(n,WRIST));
    MUJ.W(n,1) = trap_mid (w,0,15,I);
    MUJ.W(n,2) = peak_tri (w,15,30,I);
    MUJ.W(n,3) = peak_tri (w,30,45,I);
    MUJ.W(n,4) = trap_right(w,45,Inf,I);
end
[softW, hardW] = defuzz_disc_batch(dom.Wrist, MUJ.W, [1 2 3 4], true);

out.joint_soft = [softN softT softL softU softF softW];
out.joint_hard = [hardN hardT hardL hardU hardF hardW];
out.info.joint_order = {'Neck','Trunk','Leg','UpperArm','LowerArm','Wrist'};
out.info.I = I;

% ------- 2) A 表推理（twist 固定=1） -------
T_A_use = T_A(:,:,:,1);
[muA_disc, A_soft, A_hard] = infer_ifthen_coa(...
    {MUJ.U, MUJ.F, MUJ.W}, {dom.UpperArm, dom.LowerArm, dom.Wrist}, T_A_use, dom.A);
out.A_soft = A_soft; out.A_hard = A_hard;

% ------- 3) B 表推理 -------
[muB_disc, B_soft, B_hard] = infer_ifthen_coa(...
    {MUJ.N, MUJ.T, MUJ.L}, {dom.Neck, dom.Trunk, dom.Leg}, T_B, dom.B);
out.B_soft = B_soft; out.B_hard = B_hard;

% ------- 4) C 表推理（A×B→Grand）-------
[~, G_soft, ~] = infer_ifthen_coa(...
    {muA_disc, muB_disc}, {dom.A, dom.B}, T_C, dom.Grand);
out.grand_soft = G_soft;
out.grand_hard = ceil(G_soft);

out.info.tables = 'official(A/B/C), twist=1, LowerArm only {1,2} used, Wrist mapped 1..4';

% ================= 内部函数 =================
function [mu_out_disc, soft_out, hard_out] = infer_ifthen_coa(mu_in_list, dom_in_list, tableVals, domOut)
    N = size(mu_in_list{1},1);
    S = numel(domOut);
    mu_out_disc = zeros(N, S);
    D = numel(mu_in_list);

    if D==3
        A = numel(dom_in_list{1}); B = numel(dom_in_list{2}); C = numel(dom_in_list{3});
        for i=1:A, for j=1:B, for k=1:C
            concl = tableVals(i,j,k);
            str = min(mu_in_list{1}(:,i), min(mu_in_list{2}(:,j), mu_in_list{3}(:,k)));
            idx = find(domOut==concl,1);
            mu_out_disc(:,idx) = max(mu_out_disc(:,idx), str);
        end, end, end
    else
        A = numel(dom_in_list{1}); B = numel(dom_in_list{2});
        for i=1:A, for j=1:B
            concl = tableVals(i,j);
            str = min(mu_in_list{1}(:,i), mu_in_list{2}(:,j));
            idx = find(domOut==concl,1);
            mu_out_disc(:,idx) = max(mu_out_disc(:,idx), str);
        end, end
    end

    scores = domOut(:).';
    den = sum(mu_out_disc,2); den(den==0)=eps;
    soft_out = sum(mu_out_disc .* scores, 2) ./ den;
    [~,idxMax] = max(mu_out_disc, [], 2);
    hard_out = domOut(idxMax);
end

function [softVec, hardVec] = defuzz_disc_batch(ALL_SCORES, MU, usedScores, isClip)
    if nargin<4, isClip=true; end
    N=size(MU,1); K=size(MU,2);
    mask = false(1,K); mask(usedScores)=true;
    if isClip, MU(:,~mask)=0; end
    scores = 1:K; scores = scores(:).';
    softVec=zeros(N,1); hardVec=zeros(N,1);
    for t=1:N
        mu = MU(t,:);
        if sum(mu)==0
            softVec(t) = mean(usedScores);
            [~,j]=min(abs(usedScores-softVec(t)));
            hardVec(t)= usedScores(j);
        else
            softVec(t) = sum(mu.*scores)/sum(mu);
            [~,idx]=max(mu); hardVec(t)= idx;
        end
    end
end

% —— 基本隶属度函数 ——
function mu = trap_left(x, L, B, I)
    if isinf(L), L=-1e6; end
    mu = trapmf_scalar(x,[L L B-I B+I]);
end
function mu = trap_right(x, A, R, I)
    if isinf(R), R=1e6; end
    mu = trapmf_scalar(x,[A-I A+I R R]);
end
function mu = trap_mid(x, A, B, I)
    mu = trapmf_scalar(x,[A-I A+I B-I B+I]);
end
function mu = peak_tri(x, L, R, I)
    a=L-I; b=L+I; c=R-I; d=R+I;
    muL = trapmf_scalar(x,[a b b c]);
    muR = trapmf_scalar(x,[b c c d]);
    mu  = min(muL, muR);
end
function mu = trapmf_scalar(x, abcd)
    a=abcd(1); b=abcd(2); c=abcd(3); d=abcd(4);
    if x<=a || x>=d, mu=0; return; end
    if x>=b && x<=c, mu=1; return; end
    if x>a && x<b, mu=(x-a)/(b-a); return; end
    if x>c && x<d, mu=(d-x)/(d-c); return; end
end

% —— 官方 A/B/C 表（去肌力，只保留姿势与 twist） ——
function [T_A, T_B, T_C, dom] = get_rula_tables_official()
    dom.Neck     = 1:6;  dom.Trunk = 1:6;  dom.Leg = 1:2;
    dom.UpperArm = 1:6;  dom.LowerArm = 1:3; dom.Wrist = 1:4;
    dom.A = 1:9; dom.B = 1:9; dom.Grand = 1:7;

    % ---- A 表：U×L×W×Twist ----
    A_rows = {};
    % U=1
    A_rows{end+1} = [2 2 2 2 3 3 3 3];
    A_rows{end+1} = [2 2 2 2 2 3 3 3];
    A_rows{end+1} = [2 3 3 3 3 3 4 4];
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
                T_A(iu,il,iw,1) = row(2*iw-1); % twist=1
                T_A(iu,il,iw,2) = row(2*iw);   % twist=2
            end
        end
    end

    % ---- B 表：Neck×Trunk×Leg ----
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

    % ---- C 表：A×B → Grand ----
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
    T_C(1:8,1:7) = C_mat;
    T_C(9,1:7)   = C_mat(8,:); % A=9 行视作 8+
    T_C(:,8)     = T_C(:,7);   % B>=7 视作 7+
    T_C(:,9)     = T_C(:,7);
end
end
