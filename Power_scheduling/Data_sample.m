clc
clear
clean all
load pwr600.mat;
T=15;
EMY=15;%紧急状况下非必要性负荷全部脱落
Reg_T=zeros(1,15);
% Reg_T=[0,0,0,89.551879680000,231.306766800000,254.864661680000,135.240601480000,58.511278240000,0,0,0,49,85,0,0];
Reg=Reg_T(:,15-T+1:end);
pwr(601,:)=[];%删除多余的列
pwr = pwr/10;
%defined the RSD 0.1
RSD=0.1;
sigma=pwr*RSD;
rng('default') % For reproducibility
c = rng;
rng(c);
pwr_var=zeros(length(pwr),1);%波动负荷样本
%% 正态分布
for i=1:length(pwr)
        pwr_var(i)=normrnd(pwr(i),sigma(i));
end
SampleingRate = floor(size(pwr_var,1)/T);

%负荷1
for t = 1:T
    Ws_1(t,:) = pwr_var((t-1)*SampleingRate+1:t*SampleingRate)';
end
%负荷2
for t = 1:T
    if t<=EMY
        Ws_2(t,:) = pwr_var((t-1)*SampleingRate+1:t*SampleingRate)';
    else
        Ws_2(t,:) = pwr_var((t-1)*SampleingRate+1:t*SampleingRate)'/2;
    end
end
%负荷3
for t = 1:T
    if t<=EMY
        Ws_3(t,:) = pwr_var((t-1)*SampleingRate+1:t*SampleingRate)';
    else 
        Ws_3(t,:) = pwr_var((t-1)*SampleingRate+1:t*SampleingRate)'/2;
    end
end
%负荷4
for t = 1:T
    if t<=EMY
        Ws_4(t,:) = pwr_var((t-1)*SampleingRate+1:t*SampleingRate)';
    else 
        Ws_4(t,:) = pwr_var((t-1)*SampleingRate+1:t*SampleingRate)'/2;
    end
end
Ws_hat(:,:,1)=Ws_1;%HVAC load
Ws_hat(:,:,2)=Ws_2;%115AC load
Ws_hat(:,:,3)=Ws_3;%28DC load
Ws_hat(:,:,4)=Ws_4;%HVDC load
W=4;
n=40;%样本量
for i=1:W
    for t=1:T
        s = RandStream('mt19937ar');%For reproducibility
        R(t,:,i)=randsample(s,Ws_hat(t,:,i),n);
        SAA_P_d(t,:,i)=R(t,:,i);
        P_d(t,:,i)=mean(R(t,:,i));
        R(t,:,i)=R(t,:,i)-mean(R(t,:,i));
        set_upbound(t,:,i)=ceil(max(R(t,:,i),[],2));
        set_lowerbound(t,:,i)=floor(min(R(t,:,i),[],2));
    end
end
P_d=reshape(squeeze(P_d)',[],1)/1000;
set_upbound=reshape(squeeze(set_upbound)',[],1)/1000;
set_lowerbound=reshape(squeeze(set_lowerbound)',[],1)/1000;
for t=1:T
    r(4*t-3,:)=R(t,:,1)/1000;
    r(4*t-2,:)=R(t,:,2)/1000;
    r(4*t-1,:)=R(t,:,3)/1000;
    r(4*t,:)=R(t,:,4)/1000;    
end
SAA_P_d=SAA_P_d/1000;
Reg=Reg/1000;
save('data_sample.mat','P_d','set_upbound','set_lowerbound','r','SAA_P_d','Reg'); 

