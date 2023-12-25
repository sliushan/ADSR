clc
clear 
close all
A=["S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S54","G1","MB1","MB2","SB1","SB2","SB3","L1","L2","L3","L4","ESS1"];
Num_load=5;% load number
path=zeros(Num_load,size(A,2));
SN_G=13; %Serial number of generator
% Adiacency Matrix
W = [0	inf	inf	inf	inf	inf	inf	inf	inf	inf	inf	inf	1	1	inf	inf	inf	inf	inf	inf	inf	inf	inf
inf	0	inf	inf	inf	inf	inf	inf	inf	inf	inf	inf	1	inf	1	inf	inf	inf	inf	inf	inf	inf	inf
inf	inf	0	inf	inf	inf	inf	inf	inf	inf	inf	inf	inf	1	inf	inf	inf	inf	1	inf	inf	inf	inf
inf	inf	inf	0	inf	inf	inf	inf	inf	inf	inf	inf	inf	1	inf	1	inf	inf	inf	inf	inf	inf	inf
inf	inf	inf	inf	0	inf	inf	inf	inf	inf	inf	inf	inf	inf	1	inf	1	inf	inf	inf	inf	inf	inf
inf	inf	inf	inf	inf	0	inf	inf	inf	inf	inf	inf	inf	1.5	inf	inf	inf	1.5	inf	inf	inf	inf	inf
inf	inf	inf	inf	inf	inf	0	inf	inf	inf	inf	inf	inf	inf	1	inf	inf	1	inf	inf	inf	inf	inf
inf	inf	inf	inf	inf	inf	inf	0	inf	inf	inf	inf	inf	inf	inf	1	inf	inf	inf	1	inf	inf	inf
inf	inf	inf	inf	inf	inf	inf	inf	0	inf	inf	inf	inf	inf	inf	1.5	1.5	inf	inf	inf	inf	inf	inf
inf	inf	inf	inf	inf	inf	inf	inf	inf	0	inf	inf	inf	inf	inf	inf	1	inf	inf	inf	1	inf	inf
inf	inf	inf	inf	inf	inf	inf	inf	inf	inf	0	inf	inf	inf	inf	inf	inf	1	inf	inf	inf	1	inf
inf	inf	inf	inf	inf	inf	inf	inf	inf	inf	inf	0	inf	inf	inf	1	inf	inf	inf	inf	inf	inf	1
1	1	inf	inf	inf	inf	inf	inf	inf	inf	inf	inf	0	inf	inf	inf	inf	inf	inf	inf	inf	inf	inf
1	inf	1	1	inf	1.5	inf	inf	inf	inf	inf	inf	inf	0	inf	inf	inf	inf	inf	inf	inf	inf	inf
inf	1	inf	inf	1	inf	1	inf	inf	inf	inf	inf	inf	inf	0	inf	inf	inf	inf	inf	inf	inf	inf
inf	inf	inf	1	inf	inf	inf	1	1.5	inf	inf	1	inf	inf	inf	0	inf	inf	inf	inf	inf	inf	inf
inf	inf	inf	inf	1	inf	inf	inf	1.5	1	inf	inf	inf	inf	inf	inf	0	inf	inf	inf	inf	inf	inf
inf	inf	inf	inf	inf	1.5	1	inf	inf	inf	1	inf	inf	inf	inf	inf	inf	0	inf	inf	inf	inf	inf
inf	inf	1	inf	inf	inf	inf	inf	inf	inf	inf	inf	inf	inf	inf	inf	inf	inf	0	inf	inf	inf	inf
inf	inf	inf	inf	inf	inf	inf	1	inf	inf	inf	inf	inf	inf	inf	inf	inf	inf	inf	0	inf	inf	inf
inf	inf	inf	inf	inf	inf	inf	inf	inf	1	inf	inf	inf	inf	inf	inf	inf	inf	inf	inf	0	inf	inf
inf	inf	inf	inf	inf	inf	inf	inf	inf	inf	1	inf	inf	inf	inf	inf	inf	inf	inf	inf	inf	0	inf
inf	inf	inf	inf	inf	inf	inf	inf	inf	inf	inf	1	inf	inf	inf	inf	inf	inf	inf	inf	inf	inf	0
];
tic
[D, P] = floydWarshallWithPath_1(W);
toc
for i=1:Num_load
    SN_L=size(A,2)-Num_load+i;%Serial number of load
    if D(SN_G,SN_L) == inf
        path(i,:) = 0;
    else
        OPF_N = [SN_L];%OPF path node
        while OPF_N(1) ~= SN_G
            OPF_N = [P(SN_G,OPF_N(1)),OPF_N];
            path(i,1:size(OPF_N,2))=OPF_N;
        end
    end
end
%graph representation
Gr=zeros(size(A,2));
for i=1:Num_load
    j=1;
    while path(i,j+1)~=0
       Gr(path(i,j), path(i,j+1))=1;
    j=j+1;
    end
end
G = digraph(Gr,A);
plot(G,'Layout','layered','Sources',SN_G);% SN_G is the sources node
%print the OPF path
for i=1:Num_load
    OPF=path(i,:);
    OPF(OPF==0)=[];
    disp(A(1,OPF));
end