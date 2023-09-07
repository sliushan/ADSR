clc
clear
close all
%% Remarks
% 1 Generators               
% 2 Main (Primary) Buses
% 3 Secondary Buses
% 4 Load   1 ESS  
Num_G    = 1;
Num_MB   = 2;
Num_SB   = 3;
Num_Load = 4;
Num_ESS  = 1;
W=4; %3 Load
T=15;
%% Case Parameters
P_G_upper            = 1.500;   %发电机出力上界    
P_G_lower            = 0.100;      %发电机出力下界

P_MB_upper           = 1.500; %主要母线容量上界
P_SB_upper           = 1.500;    %次级母线容量上界
P_ESS_upper          = 0.100;     %ESS功率进/出上限
P_ESS_lower          = 0;    %ESS功率进/出下限
epsilon_MB = 0.95;
epsilon_SB = 0.85;
epsilon_P_ESS=0.9;
epsilon_p_d=1;
ESS_cap            = 0.800;   %SoC容量KW*h
SoC_upper            = 0.95;   %SoC上限
SoC_lower            = 0.25;   %SoC下限
SoC_initial          = 0.5;   %SoC初始
Pk_obj=150;
ESS_obj=50;
% Del_up_obj=0.25;
%% uncertainty set
if T==15
    load data_sample.mat;
    %抽取样本，按照先t后i，每个t对应三个i
    PD_mean=P_d;
    VP_positive=set_upbound;
    VP_negative=set_lowerbound;
    VP=r;   
else
    disp('the T is not right')
end
M=size(VP,2);
%% Piecewise linear decision rule
bN = 3; % 分段数
bpN = bN-1; % 断点数
bp = zeros(bN - 1,1); % 等分点向量
kl = 1+bN; %一个ksi的维度
verts = []; % set of vertices
for t = 1:T
    for i = 1:W
        step = (VP_positive((t-1)*W+i) - VP_negative((t-1)*W+i)) / bN; %计算每段区间长度
        bp = linspace(VP_negative((t-1)*W+i)+ step,VP_positive((t-1)*W+i)-step,bN-1)'; %从[a+s,b-s]中间区间取等分点，因为linspace函数只处理闭区间
        verts = [verts, vers_nonb(bN,bp,VP_positive((t-1)*W+i),VP_negative((t-1)*W+i))];
    end
end
vN = length(verts)/(W*T); %结点数
% 将所有ksi做决策处理
for m = 1:M
    ksi = VP(:,m);
    ksi_hat = [];
    ksi_hat2 = [];
    ksi_wan_now = [];
    ksi_test_now=[];
    for t = 1:T
        for i = 1:W
            step = (VP_positive((t-1)*W+i) - VP_negative((t-1)*W+i)) / bN; %计算每段区间长度
            bp = linspace(VP_negative((t-1)*W+i)+ step,VP_positive((t-1)*W+i)-step,bN-1)'; %从[a+s,b-s]中间区间取等分点，因为linspace函数只处理闭区间
            bp = [VP_negative((t-1)*W+i);bp;VP_positive((t-1)*W+i)];
            ksi_hat = L_hat(ksi((t-1)*W+i),bN,bp);
            ksi_hat2 = [ksi_hat2; L_hat(ksi((t-1)*W+i),bN,bp)];
            ksi_wan_now = [ksi_wan_now;ksi((t-1)*W+i);ksi_hat];
            ksi_test_now=[ksi_test_now;ksi_hat];
        end
    end
    ksi_test{m}=ksi_test_now;
    ksi_wan{m} = [ ones(T,1);ones(T,1);ksi_wan_now]; %按ksi0;ksi_hat0;ksi_tine0;ksi;ksi_hat;ksi_tine排列
    ksi_wan2{m} =  ksi_hat2; %按ksi_hat;ksi_tine排列
end
CC=1;
if CC==1
    % 构造wasserstein模糊集
    eta = 0.95; %设置置信度
    rho = sdpvar(1); %模糊集半径
    sum_C = 0;
    for i = 1:M
        mid = VP(:,i);
        mid = rho * norm(mid,1)^2;
        mid = exp(mid);
        sum_C = sum_C + mid;
    end
    sum_C = sum_C / M;
    obj_C = 1 + log(sum_C);
    Constraints_C =  rho >= 0;
    Objective_C = 2 *  ( ((1 / (2 * rho)) * obj_C) ^(1/2) );
    options_C = sdpsettings('verbose',0,'debug',1,'savesolveroutput',1);%, 'fmincon.TolX',1e-4
    sol_C = optimize(Constraints_C,Objective_C,options_C);
    C = sol_C.solveroutput.fmin;
    mid_eps = log( 1 / (1-eta));
    mid_eps = mid_eps / M;
    eps = C * sqrt(mid_eps);
else
    eps =1.3;
end
%% 计算ksi长度总长度
leng_k = 0;
for k = 1:T
    leng_k = leng_k+k;
end
%end of 计算ksi长度

%% variables
% 连续变量定义 每个时段变量数不同，在用到时定义
Pk = sdpvar(Num_G ,T+W*bN*leng_k); %generator
Pm = sdpvar(Num_MB,T+W*bN*leng_k); % mian bus
Pj = sdpvar(Num_SB,T+W*bN*leng_k); %second bus
Pch = sdpvar(Num_ESS,T+W*bN*leng_k); %charge
Pdis= sdpvar(Num_ESS,T+W*bN*leng_k); %discharge
Soc= sdpvar(Num_ESS,T+W*bN*leng_k); %storage
PD = sdpvar(W,T+W*bN*leng_k); %Load
%% 构造约束
% 等式约束
cons=[];
%约束1
for t = 1:T
    for n=1:Num_G
        Zac = zeros(1,1+t*W*bN); %连续型变量决策规则
        dac = 0; %常数项
        Zac = Zac + [Pk(n,t),Pk(n,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)];
        Zac = Zac -(1/epsilon_MB)* sum([Pm(1:end,t),Pm(1:end,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)],1);
        cons = [cons, Zac(1,1) == dac];
        cons = [cons, Zac(2:end)' == 0];
    end
end

%约束2-3
for t = 1:T
    % Num_MB=1
    Zac = zeros(1,1+t*W*bN); %连续型变量决策规则
    dac = -Reg(1,t); %常数项
    Zac = Zac + [Pm(1,t),Pm(1,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)];
    Zac = Zac -(1/epsilon_SB)*sum([Pj(1,t),Pj(1,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)],1)...
        -[PD(1,t),PD(1,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)];
    cons = [cons, Zac(1,1) == dac];
    cons = [cons, Zac(2:end)' == 0];
    % Num_MB=2
    Zac = zeros(1,1+t*W*bN); %连续型变量决策规则
    dac = 0; %常数项
    Zac = Zac + [Pm(2,t),Pm(2,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)];
    Zac = Zac -(1/epsilon_SB)*sum([Pj(2:3,t),Pj(2:3,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)],1);
    cons = [cons, Zac(1,1) == dac];
    cons = [cons, Zac(2:end)' == 0];
end

%约束4-5
for t = 1:T
    % Num_SB=1
    Zac = zeros(1,1+t*W*bN); %连续型变量决策规则
    dac = 0; %常数项
    Zac = Zac + [Pj(1,t),Pj(1,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)];
    Zac = Zac -[PD(2,t),PD(2,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)]...
        -[Pch(1,t),Pch(1,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)]...
         +[Pdis(1,t),Pdis(1,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)];
    cons = [cons, Zac(1,1) == dac];
    cons = [cons, Zac(2:end)'==0];
      % Num_SB=2
    Zac = zeros(1,1+t*W*bN); %连续型变量决策规则
    dac = 0; %常数项
    Zac = Zac + [Pj(2,t),Pj(2,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)];
    Zac = Zac -[PD(3,t),PD(3,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)];
    cons = [cons, Zac(1,1) == dac];
    cons = [cons, Zac(2:end)'==0];
      % Num_SB=3
    Zac = zeros(1,1+t*W*bN); %连续型变量决策规则
    dac = 0; %常数项
    Zac = Zac + [Pj(3,t),Pj(3,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)];
    Zac = Zac -[PD(4,t),PD(4,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)];
    cons = [cons, Zac(1,1) == dac];
    cons = [cons, Zac(2:end)'==0];
end
%约束11
for t = 1:T
    if t==1
    Zac = zeros(1,1+t*W*bN); %连续型变量决策规则
    dac = SoC_initial ; %常数项
    Zac = Zac + [Soc(1,t),Soc(1,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)];
    Zac = Zac -(epsilon_P_ESS/ESS_cap)*[Pch(1,t),Pch(1,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)]...
         +(1/(epsilon_P_ESS*ESS_cap))*[Pdis(1,t),Pdis(1,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)];
    cons = [cons, Zac(1,1) == dac];
    cons = [cons, Zac(2:end)'==0];
    else
    Zac = zeros(1,1+t*W*bN); %连续型变量决策规则
    dac = 0 ; %常数项
    Zac = Zac + [Soc(1,t),Soc(1,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)]-[Soc(1,t-1),Soc(1,T+1+Sumt(t-2)*W*bN:T+Sumt(t-1)*W*bN),zeros(1,W*bN)];
    Zac = Zac -(epsilon_P_ESS/ESS_cap)*[Pch(1,t),Pch(1,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)]...
         +(1/(epsilon_P_ESS*ESS_cap))*[Pdis(1,t),Pdis(1,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)];
    cons = [cons, Zac(1,1) == dac];
    cons = [cons, Zac(2:end)'==0];      
    end
end
% 不等式约束

% 不等式约束统一构造
% 构造通用的U矩阵，h向量 公用常量
Uc_l = []; %迭代生成顶点矩阵
Uc_c = []; %lamda系数矩阵
Uc = [];
Uconstant_0 = [];
Uconstant_0 = blkdiag(1, Uconstant_0);%在最开始为Uconstant矩阵补上系数1 for ksi0 常数项
Uconstant_0 = blkdiag(1, Uconstant_0);%在最开始为Uconstant矩阵补上系数1 for ksi_hat0 常数项
Ulamda_0 = [];
Ulamda_0 = blkdiag(-1, Ulamda_0); %在最开始为Ulamda矩阵补上系数-1 for ksi0 常数项
Ulamda_0 = blkdiag(-1, Ulamda_0); %在最开始为Ulamda矩阵补上系数-1 for ksi_hat0 常数项
for t = 1:T
    Uconstant = []; %顶点系数lamda组合为1约束的系数矩阵
    for w = 1:W
        Uconstant = blkdiag(Uconstant,ones(1,vN));
    end
    %change
    Ulamda = blkdiag(-verts(:,1+(t-1)*4*vN:(4*t-3)*vN),-verts(:,1+(4*t-3)*vN:(4*t-2)*vN),-verts(:,1+(4*t-2)*vN:(4*t-1)*vN),-verts(:,1+(4*t-1)*vN:4*vN*t)); %顶点凸组合表示ksi约束的系数矩阵
    Uc_l = blkdiag(Ulamda_0,Uc_l); % U矩阵上部分，由凸包顶点构成
    Uc_l = blkdiag(Uc_l,Ulamda); % U矩阵上部分，由凸包顶点构成
    Uc_c = blkdiag(Uconstant_0,Uc_c); %U矩阵下部分，由1组成实现顶点系数求和为1
    Uc_c = blkdiag(Uc_c,Uconstant); %U矩阵下部分，由1组成实现顶点系数求和为1
    Uc = [Uc_l;Uc_c]; %凸包等式约束U矩阵
    Wc_l = diag(ones(t*W*(1+bN)+2*t,1)); %W矩阵上部分，单位阵
    Wc_c = zeros(t*W+2*t,t*W*(1+bN)+2*t); %W矩阵下部分，0矩阵
    Wc = [Wc_l;Wc_c];
    hc = [zeros(t*W*(1+bN)+2*t,1);ones(t*W+2*t,1)]; %构造常数项向量 = 两倍的Z向量长度
    for n = 1:Num_G %遍历所有Num_G
        %约束6-1 upbound
        Lamda_Pkup = [sdpvar(t*W*(1+bN)+2*t,1); sdpvar(t*W+2*t,1)]; %构造等式对偶乘子lamda向量 = h向量维度
        Lamda2_Pkup = sdpvar((t*W)*vN+2*t,1); %构造不等式对偶乘子lamda向量 = verts维度
        ZPkup = -Zhat_nonb(Pk,t,T,W,n,P_G_upper ,bN);
        cons = [cons, Wc'* Lamda_Pkup == ZPkup'];
        cons = [cons, Uc'* Lamda_Pkup + Lamda2_Pkup == 0];
        cons = [cons, hc'* Lamda_Pkup >= 0];
        cons = [cons, Lamda2_Pkup >= 0];
        %约束6-2 lowerbound
        Lamda_Pklow = [sdpvar(t*W*(1+bN) + 2*t,1); sdpvar(t*W + 2*t,1)];
        Lamda2_Pklow = sdpvar((t*W)*vN + 2*t,1);
        ZPklow = Zhat_nonb(Pk,t,T,W,n,P_G_lower,bN);
        cons = [cons, Wc'* Lamda_Pklow == ZPklow'];
        cons = [cons, Uc'* Lamda_Pklow + Lamda2_Pklow == 0];
        cons = [cons, hc'* Lamda_Pklow >= 0];
        cons = [cons, Lamda2_Pklow >= 0];      
    end
   
    for n = 1:Num_MB %遍历所有Num_MB
        %约束7-1
        Lamda_PmUp = [sdpvar(t*W*(1+bN) + 2*t,1); sdpvar(t*W + 2*t,1)]; 
        Lamda2_PmUp = sdpvar((t*W)*vN + 2*t,1);
        ZPmUp = -Zhat_nonb(Pm,t,T,W,n,P_MB_upper,bN);
        cons = [cons, Wc'* Lamda_PmUp == ZPmUp'];
        cons = [cons, Uc'* Lamda_PmUp + Lamda2_PmUp == 0];
        cons = [cons, hc'* Lamda_PmUp >= 0];
        cons = [cons, Lamda2_PmUp >= 0];
        %约束7-2 lowerbound
        Lamda_Pmlow = [sdpvar(t*W*(1+bN) + 2*t,1); sdpvar(t*W + 2*t,1)];
        Lamda2_Pmlow = sdpvar((t*W)*vN + 2*t,1);
        ZPmlow = Zhat_nonb(Pm,t,T,W,n,0,bN);
        cons = [cons, Wc'* Lamda_Pmlow == ZPmlow'];
        cons = [cons, Uc'* Lamda_Pmlow + Lamda2_Pmlow == 0];
        cons = [cons, hc'* Lamda_Pmlow >= 0];
        cons = [cons, Lamda2_Pmlow >= 0];
    end
    
    
    for n = 1:Num_SB %遍历所有Num_SB
        %约束8-1
        Lamda_PjUp = [sdpvar(t*W*(1+bN) + 2*t,1); sdpvar(t*W + 2*t,1)]; 
        Lamda2_PjUp = sdpvar((t*W)*vN + 2*t,1);
        ZPjUp = -Zhat_nonb(Pj,t,T,W,n,P_SB_upper,bN);
        cons = [cons, Wc'* Lamda_PjUp == ZPjUp'];
        cons = [cons, Uc'* Lamda_PjUp + Lamda2_PjUp == 0];
        cons = [cons, hc'* Lamda_PjUp >= 0];
        cons = [cons, Lamda2_PjUp >= 0];
        %约束8-2 lowerbound
        Lamda_Pjlow = [sdpvar(t*W*(1+bN) + 2*t,1); sdpvar(t*W + 2*t,1)];
        Lamda2_Pjlow = sdpvar((t*W)*vN + 2*t,1);
        ZPjlow = Zhat_nonb(Pj,t,T,W,n,0,bN);
        cons = [cons, Wc'* Lamda_Pjlow == ZPjlow'];
        cons = [cons, Uc'* Lamda_Pjlow + Lamda2_Pjlow == 0];
        cons = [cons, hc'* Lamda_Pjlow >= 0];
        cons = [cons, Lamda2_Pjlow >= 0];
    end
   
    for n = 1:Num_ESS %遍历所有Num_ESS
        %约束9-1
        Lamda_PchUp = [sdpvar(t*W*(1+bN) + 2*t,1); sdpvar(t*W + 2*t,1)]; 
        Lamda2_PchUp = sdpvar((t*W)*vN + 2*t,1);
        ZPchUp = -Zhat_nonb(Pch,t,T,W,n,P_ESS_upper,bN);
        cons = [cons, Wc'* Lamda_PchUp == ZPchUp'];
        cons = [cons, Uc'* Lamda_PchUp + Lamda2_PchUp == 0];
        cons = [cons, hc'* Lamda_PchUp >= 0];
        cons = [cons, Lamda2_PchUp >= 0];
    end
    for n = 1:Num_ESS %遍历所有Num_ESS        
        %约束9-2
        Lamda_Pchlow = [sdpvar(t*W*(1+bN) + 2*t,1); sdpvar(t*W + 2*t,1)]; 
        Lamda2_Pchlow = sdpvar((t*W)*vN + 2*t,1);
        ZPchlow = Zhat_nonb(Pch,t,T,W,n,P_ESS_lower,bN);
        cons = [cons, Wc'* Lamda_Pchlow == ZPchlow'];
        cons = [cons, Uc'* Lamda_Pchlow + Lamda2_Pchlow == 0];
        cons = [cons, hc'* Lamda_Pchlow >= 0];
        cons = [cons, Lamda2_Pchlow >= 0];        
    end
    
    for n = 1:Num_ESS %遍历所有Num_ESS
        %约束10-1
        Lamda_PdisUp = [sdpvar(t*W*(1+bN) + 2*t,1); sdpvar(t*W + 2*t,1)]; 
        Lamda2_PdisUp = sdpvar((t*W)*vN + 2*t,1);
        ZPdisUp = -Zhat_nonb(Pdis,t,T,W,n,P_ESS_upper,bN);
        cons = [cons, Wc'* Lamda_PdisUp == ZPdisUp'];
        cons = [cons, Uc'* Lamda_PdisUp + Lamda2_PdisUp == 0];
        cons = [cons, hc'* Lamda_PdisUp >= 0];
        cons = [cons, Lamda2_PdisUp >= 0];
    end
    for n = 1:Num_ESS %遍历所有Num_ESS        
        %约束10-2
        Lamda_Pdislow = [sdpvar(t*W*(1+bN) + 2*t,1); sdpvar(t*W + 2*t,1)]; 
        Lamda2_Pdislow = sdpvar((t*W)*vN + 2*t,1);
        ZPdislow = Zhat_nonb(Pdis,t,T,W,n,P_ESS_lower,bN);
        cons = [cons, Wc'* Lamda_Pdislow == ZPdislow'];
        cons = [cons, Uc'* Lamda_Pdislow + Lamda2_Pdislow == 0];
        cons = [cons, hc'* Lamda_Pdislow >= 0];
        cons = [cons, Lamda2_Pdislow >= 0];  
    end
    %约束11
    for n = 1:Num_ESS %遍历所有Num_ESS        
        %约束12-1
        Lamda_SOCupper = [sdpvar(t*W*(1+bN) + 2*t,1); sdpvar(t*W + 2*t,1)]; 
        Lamda2_SOCupper = sdpvar((t*W)*vN + 2*t,1);
        ZSOCupper = -Zhat_nonb(Soc,t,T,W,n,SoC_upper,bN);
        cons = [cons, Wc'* Lamda_SOCupper == ZSOCupper'];
        cons = [cons, Uc'* Lamda_SOCupper + Lamda2_SOCupper == 0];
        cons = [cons, hc'* Lamda_SOCupper >= 0];
        cons = [cons, Lamda2_SOCupper >= 0];  
    end        
    for n = 1:Num_ESS %遍历所有Num_ESS        
        %约束12-2
        Lamda_SOClow = [sdpvar(t*W*(1+bN) + 2*t,1); sdpvar(t*W + 2*t,1)]; 
        Lamda2_SOClow = sdpvar((t*W)*vN + 2*t,1);
        ZSOClow = Zhat_nonb(Soc,t,T,W,n,SoC_lower,bN);
        cons = [cons, Wc'* Lamda_SOClow == ZSOClow'];
        cons = [cons, Uc'* Lamda_SOClow + Lamda2_SOClow == 0];
        cons = [cons, hc'* Lamda_SOClow >= 0];
        cons = [cons, Lamda2_SOClow >= 0];  
    end    
    
   
    %约束13
    for n=1:Num_Load
        Lamda_D = [sdpvar(t*W*(1+bN) + 2*t,1); sdpvar(t*W + 2*t,1)];
        Lamda2_D = sdpvar((t*W)*vN + 2*t,1);
        AD = [PD(n,t)-PD_mean((t-1)*W+n,1), PD(n,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)];
        ZD = [];
        for tt = 1:t
            if tt == t
                ZD = [ZD,0,AD(1,1)];
            else
                ZD = [ZD,0,0];
            end
        end
        for tt = 1:t
            for w = 1:W
                if (tt==t && w==n)
                    ZD = [ZD,-1,AD(1,2+(tt-1)*W*bN+(w-1)*bN:1+(tt-1)*W*bN+w*bN)];
                else
                    ZD = [ZD,0,AD(1,2+(tt-1)*W*bN+(w-1)*bN:1+(tt-1)*W*bN+w*bN)];
                end
            end
        end
        cons = [cons, Wc'* Lamda_D == ZD'];
        cons = [cons, Uc'* Lamda_D + Lamda2_D == 0];
        cons = [cons, hc'* Lamda_D >= 0];
        cons = [cons, Lamda2_D >= 0];
    end
 
end

% 构造目标函数约束

%成本函数
Aobj = zeros(1,T*W*bN); %ksi_hat系数矩阵
Aobj0 = zeros(1,T); %ksi_hat0系数矩阵
Uobj_l = []; %目标函数凸包约束等式约束顶点矩阵
Uobj_c = []; %目标函数凸包约束等式约束lamda求和为1系数矩阵
for t = 1:T
    for n = 1:Num_G
        Aobj0 = Aobj0+ Pk_obj*[zeros(1,t-1),Pk(n,t),zeros(1,T-t)]+ESS_obj*[zeros(1,t-1),Pch(n,t)+Pdis(n,t),zeros(1,T-t)];
        Aobj = Aobj + Pk_obj*[Pk(n,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN),zeros(1,(T-t)*W*bN)]...
            + ESS_obj*[Pch(n,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)+Pdis(n,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN), zeros(1,(T-t)*W*bN)];
    end
    Ulamda = blkdiag(verts(:,1+(t-1)*4*vN:(4*t-3)*vN),verts(:,1+(4*t-3)*vN:(4*t-2)*vN),verts(:,1+(4*t-2)*vN:(4*t-1)*vN),verts(:,1+(4*t-1)*vN:4*vN*t)); %顶点凸组合表示ksi约束的系数矩阵
    Uconstant = []; %顶点系数lamda组合为1约束的系数矩阵
    for w = 1:W
        Uconstant = blkdiag(Uconstant,ones(1,vN));
    end
    Uobj_l = blkdiag(Ulamda_0,Uobj_l);
    Uobj_l = blkdiag(Uobj_l,-Ulamda);
    Uobj_c = blkdiag(Uconstant_0,Uobj_c);
    Uobj_c = blkdiag(Uobj_c,Uconstant);
end
U = [Uobj_l;Uobj_c];
W_l = diag(ones(T*W*kl+2*t,1));
W_c = zeros(T*W+2*t,T*W*kl+2*t);
W_obj = [W_l;W_c];
Zobj = [];
Zobj_0 = [];
for t = 1:T
    Zobj_0 = [Zobj_0,0,Aobj0(1,t)];
    for w = 1:W
        %change
        Zobj = [Zobj,0,Aobj(1,(4*(t-1)+w-1)*bN+1:(4*(t-1)+w)*bN)];
    end
end
Zobj = [Zobj_0,Zobj]; %ksi0,ksi_hat0,ksi_tine0,ksi,ksi_hat,ksi_tine

hobj = [zeros(1,T*W*kl+2*t),ones(1,T*W+2*t)]'; %凸包常数项
v = sdpvar(1,M,'full'); %构造第一次lagrangian等式约束乘子v
beta = sdpvar(1,1,'full'); %构造第一次lagrangian不等式约束乘子beta
y1 = sdpvar((T*W*(kl+1))+2*2*t,M,'full');%构造第二次等式lagrangian乘子y1，与h同维，每个样本都有一个这样的乘子向量
constraints = [];
constraints = [constraints, beta >= 0]; %不等式乘子要大于等于0
for m = 1:M %开始构造目标函数转换带出的约束
    y_j = y1(:,m); 
    mid_obj = Zobj' - W_obj'*y_j;
    constraints = [constraints, (mid_obj' * ksi_wan{m} + hobj' * y_j) <= M*v(1,m)];
    constraints = [constraints, abs(mid_obj) <= beta];
    constraints = [constraints, U' * y_j >= 0];
end
% end of 构造目标函数约束


constraints = [constraints, cons]; %添加物理约束

objective = sum(v,2) + eps*beta; %设置目标函数
options = sdpsettings('verbose',2,'solver','gurobi','debug',1,'savesolveroutput',1,'savesolverinput',1);
sol = optimize(constraints,objective,options);

% Analyze error flags
if sol.problem == 0
    Obj = value(objective);
    result_Pk=value(Pk);
    result_Pm=value(Pm);
    result_Pj=value(Pj);
    result_Pch=value(Pch);
    result_Pdis=value(Pdis);
    result_Soc=value(Soc);  
     result_PD=value(PD);    
    
    disp(sol.solvertime);
else
    disp('Oh shit!, something was wrong!');
    sol.info
    yalmiperror(sol.problem)
end
save('Inspection_data.mat','T','result_Soc','result_Pk','result_Pm',...
     'result_Pj','result_Pch','result_Pdis','result_PD');



