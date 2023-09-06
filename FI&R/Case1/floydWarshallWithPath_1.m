function [D, P] = floydWarshallWithPath_1(W)
% W是邻接矩阵，表示有向图或无向图
% i和j是起点和终点的编号
% D是距离矩阵，P是前驱矩阵，path是从i到j的最短路径
n = size(W,1);
D = W; % 初始化距离矩阵为邻接矩阵
P = repmat((1:n)',1,n); % 初始化前驱矩阵为每个点到自身的路径
for k = 1:n
    for i = 1:n
        for j = 1:n
            if D(i,k) + D(k,j) < D(i,j)
                % 如果经过k点可以缩短i到j的距离，更新距离和前驱矩阵
                D(i,j) = D(i,k) + D(k,j);
                P(i,j) = P(k,j);
            end
        end
    end
end
end
