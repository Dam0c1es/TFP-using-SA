function BestX = FindBestX(Y)
%FINDBESTX 找到最优的X
global N;
global S;
T1=1000;
K1=0.98;
L1=100;
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
PreX=X;
BestX=FindX(Y);
count=0;
pp=1;%记录取值次数
while T1>0.01
    count=count+1;
    for i=1:length(L1)
        NextX=FindX(Y);
        Z2_1=CalZ2(NextX);
        Z2_2=CalZ2(BestX);
        if Z2_1<Z2_2
            preBestX=BestX;
            BestX=NextX;
        end
        if CalZ2(NextX)<CalZ2(PreX)
            PreX=NextX;
            pp=pp+1;
        else
            changer=-1*(CalZ2(NextX)-CalZ2(PreX))/T1;
            p=exp(changer);
            if p>rand
                PreX=NextX;
                pp=pp+1;
            end
        end
        trace(pp)=CalZ2(BestX);
    end
    T1=K1*T1;
end


end

