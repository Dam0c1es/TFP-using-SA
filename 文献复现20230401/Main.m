clc,clear
close all;
global N;
global S;
global d;
global Dt;
global Biao_1;
global D;
count1=1;
%启发式算法找第2层最优解
%%%%%%%%%%%%%%%%%%导入数据%%%%%%%%%%%%%%%%%%%
Biao_1=xlsread("表1 沿途7个支点站的主要参数.xlsx");%第二列为c_k,第3列为tao_k,第4列为R_k,第五列为H_k,第6列T_k
Biao_2=xlsread("表2 支点站之间车流矩阵.xlsx");%OD矩阵
Biao_3=xlsread("表3 相邻2个支点站之间管内列车的参数.xlsx");%相关参数，分上下行，分别为区段旅行时间、摘挂旅行时间、区段车流临界值
load Zhanming.mat;
Biao_1(:,2)=[];%删除Nan
N=Biao_2';%OD矩阵
%%%%%%%%%%%%%%%%%%输入模拟退火初始参数%%%%%%%%%
T=100;%初始温度
L=100;%Markov链长度
K=0.9;%退火系数
Yz=10^(-5);%容差
P=0;%记录取值次数
%%%%%%%%%%%%%%%%%%%初始化模型与决策变量%%%%%%%
%%集合S_ij:i——j中的车站，不包含i,j
for i=1:length(N)
    for j=i:length(N)
        m=find(Biao_1(:,1)==i);%点i的所在行号
        n=find(Biao_1(:,1)==j);%点j所在行号
        l=[];
        for ii=m:n
            l=[l;Biao_1(ii,1)];
        end
        l(l==m)=[];
        l(l==n)=[];
        S{i}{j}=l;
        S{j}{i}=fliplr(l);
    end
end
%%
% 确定越阶函数Epsilon
%首先确定临界矩阵Dt,i-j的临界值就是i-j区段的最小值
%E是Dij^T矩阵
for i=1:length(N)
    for j=1+i:length(N)
        m=find(Biao_3(:,1)==i);%点i的所在行号
        n=find(Biao_3(:,2)==j);%点j所在行号
        ll1=[];
        ll2=[];
        for ii=m:n
            ll1=[ll1;Biao_3(ii,5)];%下行方向
            ll2=[ll2,Biao_3(ii,8)];%上行方向
        end
        if ~isempty(ll1) || ~isempty(ll2)
            Dt(i,j)=min(ll1);%临界矩阵
            Dt(j,i)=min(ll2);
        else
            Dt(i,j)=0;
            Dt(j,i)=0;
        end
    end
end

%%
 %确定一组初始解y_ij
for i=1:length(N)
    for j=1:length(N)
        if i~=j
            Y(i,j)=1-Epsilon(i,j,N,Dt);%一组初始解
        else
            Y(i,j)=0;
        end
    end
end
%%确定中间变量dik
for i=1:length(N)
    for k=1:length(N)
        if Y(i,k)==1
            d(i,k)=Biao_1(i,k)*Y(i,k);
        else
            d(i,k)=400;
        end
    end
end
%%
%%初始化初始化x^k_i,j
for i=1:length(N)
    for j=1:length(N)
        if i~=j
            x{i}{j}=zeros(1,length(S{i}{j}));%初始化x^k_i,j
        else
            x{i}{j}=[];
        end
    end
end

for i=1:length(N)
    for j=1:length(N)
        if Y(i,j)==0 && ~isempty(S{i}{j})
            p=randi([1,length(S{i}{j})],1,1);
            x{i}{j}(p)=1;
        end
    end
end
T1=10000;
K1=0.98;
L1=100;
%初始化目标函数
Z2=0;
for i=1:length(N)
    for j=1:length(N)
        if i~=j
            if ~isempty(x{i}{j})
                for k=1:length(x{i}{j})
                    Z2=Z2+d(i,S{i}{j}(k))*x{i}{j}(k);%目标函数初始化完成
                end
            end
        end
    end
end
%%%%%注：点k在S{i}{j}(k)搜索就行
%模拟退火求解下层最优
PreX=x;
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
plot(trace);
title("Z2最小值：",num2str(CalZ2(BestX)));


%%
%计算车流矩阵f
f=zeros(size(N));
for i=1:length(N)
    for j=1:length(N)
        m=0;
        for s=1:length(N)
            if ~isempty(find(S{s}{j}==i))
                n=find(S{s}{j}==i);
                m=m+f(s,j)*x{s}{j}(n);
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
                m=m+x{i}{j}(n)*f(i,j);
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
                m=m+f(i,l)*x{i}{l}(n);
            end
        end
        D(i,j)=f(i,j)*Y(i,j);
    end
end

%惩罚函数项
M1=F-Biao_1(:,4);
M2=round(sum(D,2)/200)-Biao_1(:,5);
BestZ1=func1(Biao_1,N,Y,Dt,D,F,M1,M2);

%%
%模拟退火过程
OmegaA=Y;%储存所有可能去向yij
% OmegaB=zeros(size(N));%储存已经确定去向的yij
% for i=1:length(N)
%     for j=1:length(N)
%         if Y(i,j)==1
%             OmegaB(i,j)=1;
%         end
%     end
% end%储存集合初始化完成
BestY=Y;%当前Y设置为最好Y
PreY=Y;
PreX=BestX;

%开始模拟退火过程
m=0;
n=0;
pp1=1;%取值次数
while T>15
    for ii=1:L%Markov随机搜索
        NextY=PreY;
        while m==n
            m=randi([1,length(N)],1,1);
            n=randi([1,length(N)],1,1);
        end
        if NextY(m,n)==1
            NextY(m,n)=0;
        else
            NextY(m,n)=1;
        end
        NextbestX=FindBestX(NextY);%在Next情况下求得最优X
        %%%%%%%%计算当前情况下能量函数%%%%%%%%%
        Cal_1=CalZ1(NextbestX,NextY);
        Cal_2=CalZ1(BestX,BestY);
        if  Cal_1<Cal_2
            preBestX=BestX;
            preBestY=BestY;
            BestX=NextbestX;
            BestY=NextY;
        end

        if CalZ1(NextbestX,NextY)<CalZ1(PreX,PreY)
            PreX=NextbestX;
            PreY=NextY;
            pp1=pp1+1;
        else
            changer=-1*(CalZ1(NextbestX,NextY)-CalZ1(PreX,PreY));
            p=exp(changer);
            if p>rand
                PreX=NextbestX;
                PreY=NextY;
                pp1=pp1+1;
            end
        end
        trace2(pp1)=CalZ1(BestX,BestY);
    end
    T=K*T;
    trace3(count1)=CalZ1(BestX,BestY);
    figure(2);
    title(num2str(count1),num2str(CalZ1(BestX,BestY)));
    plot(trace3);
    pause(0.05);
    hold on;
    count1=count1+1;
end
figure(3);
plot(trace2);
title("上层最优为Z1：",num2str(CalZ1(BestX,BestY)))
            





            

        














