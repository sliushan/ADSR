 %%  Zhat_nonb  by yy 2022.4.10
 %对变量进行排序
function [Z2,Z3] = Zhat_nonb_withoutlift(Y,t,T,W,n,bN)
Z2 =Y(n,t);
Z3 = Y(n,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN);
end