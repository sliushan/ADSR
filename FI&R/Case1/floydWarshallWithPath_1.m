function [D, P] = floydWarshallWithPath_1(W)
% W���ڽӾ��󣬱�ʾ����ͼ������ͼ
% i��j�������յ�ı��
% D�Ǿ������P��ǰ������path�Ǵ�i��j�����·��
n = size(W,1);
D = W; % ��ʼ���������Ϊ�ڽӾ���
P = repmat((1:n)',1,n); % ��ʼ��ǰ������Ϊÿ���㵽�����·��
for k = 1:n
    for i = 1:n
        for j = 1:n
            if D(i,k) + D(k,j) < D(i,j)
                % �������k���������i��j�ľ��룬���¾����ǰ������
                D(i,j) = D(i,k) + D(k,j);
                P(i,j) = P(k,j);
            end
        end
    end
end
end
