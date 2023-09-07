%本例主要为目标函数分段线性化，将工作区间分为三类,引入0-1变量进行分离
%
%
clc
clear
close all
Num_G    = 1;
Num_MB   = 2;
Num_SB   = 3;
Num_Load = 4;
Num_ESS  = 1;
W=4; %4 Load
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
ESS_cap            = 0.8;   %SoC容量KW*h
SoC_upper            = 0.95;   %SoC上限
SoC_lower            = 0.25;   %SoC下限
SoC_initial          = 0.5;   %SoC初始
Pk_obj=150;
ESS_obj=50;
%设备可用性
if T==15
    load data_sample.mat; 

    %抽取样本，按照先t后i，每个t对应三个i
    PD_mean=P_d;
    PD=PD_mean;
    VP_positive=set_upbound;
    PD=PD_mean+VP_positive;%不确定变量 u 的取值应为所描述的波动区间的边界
   
else
    disp('the T is not right')
end

%% variables
% 连续变量定义 每个时段变量数不同，在用到时定义
Pk = sdpvar(Num_G ,T,'full'); %generator
Pm = sdpvar(Num_MB,T,'full'); % mian bus
Pj = sdpvar(Num_SB,T,'full'); %second bus
Pch = sdpvar(Num_ESS,T,'full'); %charge
Pdis= sdpvar(Num_ESS,T,'full'); %discharge
Soc= sdpvar(Num_ESS,T,'full'); %storage
%% cons
cons=[];
% end
for t=1:T
    %约束1
    cons=[cons,Pk(1,t)-(1/epsilon_MB)* sum(Pm(1:end,t),1)==0];
    
    %约束2-3
    cons=[cons,Reg(1,t)+Pm(1,t)-(1/epsilon_SB)*sum(Pj(1,t),1)-PD(4*t-3,1)==0];%HVAC1
    cons=[cons,Pm(2,t)-(1/epsilon_SB)*sum(Pj(2:3,t),1)==0];
    
    %约束4-5
    cons=[cons,Pj(1,t)-(PD(4*t-2,1)+Pch(1,t)-Pdis(1,t))==0];
    cons=[cons,Pj(2,t)-PD(4*t-1,1)==0];
    cons=[cons,Pj(3,t)-PD(4*t,1)==0];
    
    %约束6-10
    
    cons=[cons,Pk(1,t)<=P_G_upper];
    cons=[cons,P_G_lower<=Pk(1,t)];
    
    for n=1:Num_MB
        cons=[cons,Pm(n,t)<=P_MB_upper];
        cons=[cons,0<=Pm(n,t)];        
    end
    
    for n=1:Num_SB
        cons=[cons,Pj(n,t)<=P_SB_upper];
        cons=[cons,0<=Pj(n,t)];
    end 
    cons=[cons,Pch(1,t)<=P_ESS_upper];
    cons=[cons,Pdis(1,t)<=P_ESS_upper];
    cons=[cons,0<=Pch(1,t)];
    cons=[cons,0<=Pdis(1,t)];
    
    
    %约束11
    
    if t==1
        cons=[cons,Soc(1,t)-SoC_initial-(epsilon_P_ESS/ESS_cap)*Pch(1,t)+(1/(epsilon_P_ESS*ESS_cap))*Pdis(1,t)==0];
    else
        cons=[cons,Soc(1,t)-Soc(1,t-1)-(epsilon_P_ESS/ESS_cap)*Pch(1,t)+(1/(epsilon_P_ESS*ESS_cap))*Pdis(1,t)==0];
    end
    
    %约束12
    
    cons=[cons,Soc(1,t)<=SoC_upper];
    cons=[cons,SoC_lower<=Soc(1,t)];
end
objective=Pk_obj*sum(Pk,'all')+ESS_obj*sum(Pch+Pdis,'all');
options = sdpsettings('verbose',2,'solver','gurobi','debug',1,'savesolveroutput',1,'savesolverinput',1);
sol = optimize(cons,objective,options);

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
     min(result_Pk)
     max(result_Pk)
    disp(sol.solvertime);
else
    disp('Oh shit!, something was wrong!');
    sol.info
    yalmiperror(sol.problem)
end
save('Deter.mat','T','result_Soc','result_Pk','result_Pm',...
     'result_Pj','result_Pch',...
     'result_Pdis');


