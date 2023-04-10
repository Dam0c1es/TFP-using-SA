function BestZ1 = CalZ1(X,Y)
%CALZ1 计算上层函数
%   此处显示详细说明
global Biao_1;
global Dt;
global S;
global N;
global D;
f=zeros(size(N));
for i=1:length(N)
    for j=1:length(N)
        m=0;
        for s=1:length(N)
            if ~isempty(find(S{s}{j}==i))
                n=find(S{s}{j}==i);
                m=m+f(s,j)*X{s}{j}(n);
            end
        end
        f(i,j)=N(i,j)+m;
    end
end
%计算车站处理车数Fk
for k=1:length(N)
    m=0;
    for i=1:length(N)
        for j=1:length(N)
            if ~isempty(find(S{i}{j}==k, 1))
                n=find(S{i}{j}==k);
                m=m+X{i}{j}(n)*f(i,j);
            end
        end
    end
    F(k)=m;
end
%计算吸引车流量矩阵D
for i=1:length(N)
    for j=1:length(N)
        m=0;
        for l=1:length(N)
            if ~isempty(find(S{i}{l}==j, 1))
                n=find(S{i}{l}==j);
                m=m+f(i,l)*X{i}{l}(n);
            end
        end
        D(i,j)=f(i,j)*Y(i,j);
    end
end

%惩罚函数项
M1=F-Biao_1(:,4);
M2=round(sum(D,2)/200)-Biao_1(:,5);
BestZ1=func1(Biao_1,N,Y,Dt,D,F,M1,M2);


end

