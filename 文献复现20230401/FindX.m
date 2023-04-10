function X = FindX(Y)
global N;
global S;
%FINDX 在Y确定的情况下选择X
%   此处显示详细说明
for i=1:length(N)
    for j=1:length(N)
        if i~=j
            X{i}{j}=zeros(1,length(S{i}{j}));%初始化x^k_i,j
        else
            X{i}{j}=[];
        end
    end
end

for i=1:length(N)
    for j=1:length(N)
        if Y(i,j)==0 && ~isempty(S{i}{j})
            p=randi([1,length(S{i}{j})],1,1);
            X{i}{j}(p)=1;
        end
    end
end
end

