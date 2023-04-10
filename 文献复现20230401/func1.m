function f = func1(Biao_1,N,Y,E,D,F,M1,M2)
%FUNC1 用于计算上层的函数
% Cij^difference  此处定为4h
%Y 输入区段决策变量
%E 越阶函数计算矩阵
%D i——j的吸引车流量
%F 车站k的改编量
%mij 平均编成辆数55
global Dt;
f=0;
for i=1:length(N)
    for j=1:length(N)
        f=f+Biao_1(i,2)*55*Y(i,j);%第一项
    end
end
for k=1:length(N)
    f=f+F(k)*Biao_1(k,3);%第二项
end
%第三项
Mu=400*max(0,max(max(M1)))+200*max(0,max(M2));
f=f+Mu;
%第4项
m=0;
for i=1:length(N)
    for j=1:length(N)
        m=Epsilon(i,j,D,Dt)*4*D(i,j);
    end
end
f=f+m;

end

