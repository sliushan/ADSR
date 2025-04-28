clc
clear
close all
%% Remarks
% 1 Main generators              
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
P_G_upper            = 1.5;   %发电机出力上界    
P_G_lower            = 0.1;      %发电机出力下界

P_MB_upper           = 1.5; %主要母线容量上界
P_SB_upper           = 1.5;    %次级母线容量上界
P_ESS_upper          = 0.1;     %ESS功率进/出上限
P_ESS_lower          = 0;    %ESS功率进/出下限
epsilon_MB = 0.95;
epsilon_SB = 0.85;
epsilon_P_ESS=0.9;
epsilon_p_d=1;
ESS_cap            = 0.8;   %SoC容量
SoC_upper            = 0.95;   %SoC上限
SoC_lower            = 0.25;   %SoC下限
SoC_initial          = 0.5;   %SoC初始
Pk_obj=150;
ESS_obj=50;
%% uncertainty set
if T==15
    load data_sample.mat;
    %抽取样本，按照先t后i，每个t对应三个i
    PD_mean=P_d;
    VP_positive=set_upbound;
    VP_negative=set_lowerbound;
    VP=r;
    % ksi 支撑集
    ksi_positive = [1;VP_positive];
    ksi_negative = [1;VP_negative];   
else
    disp('the T is not right')
end
M=size(VP,2);
%% decision rule
bN = 1; 
CC=0;
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
    eps =1.5;
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
% 约束1
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
    cons = [cons, Zac(2:end)'==0];
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


Uc_l = []; %上界
Uc_c = []; %下界

for t = 1:T
    Ulamda = eye(W); %顶点凸组合表示ksi约束的系数矩阵
    Uc_l = blkdiag(Ulamda,Uc_l); % U矩阵上部分，由凸包顶点构成
    Uc_c = blkdiag(-Ulamda,Uc_c); %U矩阵下部分，由1组成实现顶点系数求和为1
    Uc=[Uc_l;Uc_c];
    d_l=VP_positive(1:W*t,1);
    d_c=-VP_negative(1:W*t,1);
    d=[d_l;d_c];
    for n = 1:Num_G %遍历所有Num_G
        %约束6-1 upbound
        Lamda_Pkup = sdpvar(1,2*W*t,'full'); %构造等式对偶乘子lamda向量 = h向量维度
        [ZPkup1, ZPkup2]=Zhat_nonb_withoutlift(Pk,t,T,W,n,bN);
        cons=[cons,ZPkup1+Lamda_Pkup*d<=P_G_upper ];
        cons = [cons, ZPkup2-Lamda_Pkup*Uc == 0];
        cons = [cons, Lamda_Pkup>= 0];
    end
    for n = 1:Num_G %遍历所有Num_G
        %约束6-2 lowerbound
        Lamda_Pklow =sdpvar(1,2*W*t,'full');
        [ZPklow1,ZPklow2] = Zhat_nonb_withoutlift(Pk,t,T,W,n,bN);
        cons = [cons, ZPklow1-Lamda_Pklow*d>=P_G_lower];
        cons = [cons, ZPklow2+Lamda_Pklow*Uc == 0];
        cons = [cons, Lamda_Pklow >= 0];
    end

    for n = 1:Num_MB %遍历所有Num_MB
        %约束7
        Lamda_PmUp = sdpvar(1,2*W*t,'full'); 
        [ZPmUp1,ZPmUp2] = Zhat_nonb_withoutlift(Pm,t,T,W,n,bN);
        cons = [cons, ZPmUp1+Lamda_PmUp*d<=P_MB_upper];
        cons = [cons, ZPmUp2-Lamda_PmUp*Uc == 0];
        cons = [cons, Lamda_PmUp>= 0];
    end
    for n = 1:Num_SB %遍历所有Num_SB
        %约束8
        Lamda_PjUp = sdpvar(1,2*W*t,'full'); 
        [ZPjUp1,ZPjUp2] = Zhat_nonb_withoutlift(Pj,t,T,W,n,bN);
        cons = [cons,ZPjUp1+Lamda_PjUp *d<=P_SB_upper];
        cons = [cons,ZPjUp2-Lamda_PjUp*Uc == 0];
        cons = [cons,  Lamda_PjUp>= 0];
    end
   
    for n = 1:Num_ESS %遍历所有Num_ESS
        %约束9-1
        Lamda_PchUp = sdpvar(1,2*W*t,'full'); 
        [ZPchUp1,ZPchUp2] = Zhat_nonb_withoutlift(Pch,t,T,W,n,bN);
        cons = [cons, ZPchUp1+Lamda_PchUp *d<=P_ESS_upper];
        cons = [cons, ZPchUp2-Lamda_PchUp*Uc == 0];
        cons = [cons,Lamda_PchUp>= 0];
    end
    for n = 1:Num_ESS %遍历所有Num_ESS        
        %约束9-2
        Lamda_Pchlow = sdpvar(1,2*W*t,'full'); 
        [ZPchlow1,ZPchlow2] = Zhat_nonb_withoutlift(Pch,t,T,W,n,bN);
        cons = [cons,ZPchlow1-Lamda_Pchlow *d>=P_ESS_lower];     
        cons = [cons,  ZPchlow2+Lamda_Pchlow*Uc == 0];
        cons = [cons,Lamda_Pchlow >= 0];     
    end
         
    for n = 1:Num_ESS %遍历所有Num_ESS
        %约束10-1
        Lamda_PdisUp = sdpvar(1,2*W*t,'full'); 
        [ZPdisUp1, ZPdisUp2] = Zhat_nonb_withoutlift(Pdis,t,T,W,n,bN);
        cons = [cons, ZPdisUp1+Lamda_PdisUp  *d<=P_ESS_upper];
        cons = [cons,  ZPdisUp2-Lamda_PdisUp *Uc == 0];
        cons = [cons, Lamda_PdisUp >= 0];
    end
    for n = 1:Num_ESS %遍历所有Num_ESS        
        %约束10-2
        Lamda_Pdislow = sdpvar(1,2*W*t,'full'); 
        [ZPdislow1,ZPdislow2] = Zhat_nonb_withoutlift(Pdis,t,T,W,n,bN);
        cons = [cons,ZPdislow1-Lamda_Pdislow *d>=P_ESS_lower];
        cons = [cons,  ZPdislow2+Lamda_Pdislow*Uc == 0];
        cons = [cons,Lamda_Pdislow >= 0];  
    end
    for n = 1:Num_ESS %遍历所有Num_ESS        
        %约束12-1
        Lamda_SOCupper =  sdpvar(1,2*W*t,'full'); 
        [ZSOCupper1,ZSOCupper2] = Zhat_nonb_withoutlift(Soc,t,T,W,n,bN);
        cons = [cons, ZSOCupper1+Lamda_SOCupper  *d<=SoC_upper,];
        cons = [cons,  ZSOCupper2-Lamda_SOCupper *Uc == 0];
        cons = [cons, Lamda_SOCupper >= 0];
    end        
    for n = 1:Num_ESS %遍历所有Num_ESS        
        %约束12-2
        Lamda_SOClow =  sdpvar(1,2*W*t,'full'); 
        [ZSOClow1, ZSOClow2] = Zhat_nonb_withoutlift(Soc,t,T,W,n,bN);
        cons = [cons,ZSOClow1-Lamda_SOClow *d>=SoC_lower];
        cons = [cons,  ZSOClow2+Lamda_SOClow*Uc == 0];
        cons = [cons,Lamda_SOClow>= 0];  
    end    

   
    %约束13
    for n=1:Num_Load
        Lamda_D = sdpvar(1,2*W*t,'full');
        AD1 = PD(n,t)-PD_mean((t-1)*W+n,1);
        AD2=PD(n,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN);
        ZD = zeros(1,W*t);
        if n==1
        ZD(1,4*t-3)=-1;  
        elseif n==2
        ZD(1,4*t-2)=-1;                  
        elseif n==3
        ZD(1,4*t-1)=-1;   
        else
         ZD(1,4*t)=-1;           
        end
        AD2=AD2+ZD;
        cons = [cons, AD1-Lamda_D*d>=0];       
        cons = [cons, AD2+Lamda_D*Uc == 0];
        cons = [cons, Lamda_D>= 0];
    end
end

%成本函数
Aobj0 = 0; %ksi_hat系数矩阵
Aobj = zeros(1,T*W); %ksi_hat0系数矩阵

for t = 1:T
    for n = 1:Num_G
        Aobj0 = Aobj0+ Pk_obj*Pk(n,t);
        Aobj = Aobj + Pk_obj*[Pk(n,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN),zeros(1,(T-t)*W*bN)]; 
        if n==1
            Aobj0 = Aobj0+ESS_obj*Pch(n,t)+ESS_obj*Pdis(n,t);
            Aobj = Aobj +ESS_obj*[Pch(n,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN) + Pdis(n,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN), zeros(1,(T-t)*W*bN)]; %将所有Pch、Pdis和Pm决策系数对应加起来
        end
    end
end



v = sdpvar(1,M,'full'); %构造第一次lagrangian等式约束乘子v
beta = sdpvar(1,1,'full'); %构造第一次lagrangian不等式约束乘子beta
y1 = sdpvar(T*W*2,M,'full');%构造第二次等式lagrangian乘子y1，与h同维，每个样本都有一个这样的乘子向量
constraints = [];
constraints = [constraints, beta >= 0]; %不等式乘子要大于等于0
for m = 1:M %开始构造目标函数转换带出的约束
    y_j = y1(:,m); 
    mid_obj = Aobj - y_j'*Uc;
    constraints = [constraints, y_j'*d+ mid_obj *VP(:,m)+Aobj0 <= M*v(1,m)];
    constraints = [constraints,  y_j >= 0];
    constraints = [constraints, abs(mid_obj) <= beta];
end
% end of 构造目标函数约束
constraints = [constraints, cons]; %添加物理约束
objective= sum(v,2) + eps*beta; %设置目标函数
options = sdpsettings('verbose',2,'solver','gurobi','debug',1,'savesolveroutput',1,'savesolverinput',1);
sol = optimize(constraints,objective,options);

% Analyze error flags
if sol.problem == 0
    % Extract and display value
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
save('Inspection_data_withoutlifting.mat','result_Soc','result_Pk','result_Pm',...
     'result_Pj','result_Pch','result_PD','result_Pdis');


