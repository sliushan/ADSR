clc
clear 
close all

% Adiacency Matrix
% W = [0,1,1;
%     1,0,inf;
%     1,inf,0;  
% ];
% W = [0,1.5,1.5;
%     1.5,0,inf;
%     1.5,inf,0;  
% ];
W = [0,inf,1,inf,1;
    inf,0,inf,1.5,1.5;
    1,inf,0,inf,inf;
     inf,1.5,inf,0,inf;   
    1,1.5,inf,inf,0;  
];
tic
[D, P] = floydWarshallWithPath_1(W);
toc
