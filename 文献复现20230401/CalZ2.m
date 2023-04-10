function f = CalZ2(x)
%CALZ2 计算Z2
%   此处显示详细说明
global N;
global S;
global d;
f=0;
for i=1:length(N)
    for j=1:length(N)
        if i~=j
            if ~isempty(x{i}{j})
                for k=1:length(x{i}{j})
                    f=f+d(i,S{i}{j}(k))*x{i}{j}(k);%目标函数初始化完成
                end
            end
        end
    end
end
end

