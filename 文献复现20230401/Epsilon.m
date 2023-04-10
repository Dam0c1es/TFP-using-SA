function f = Epsilon(i,j,D,Dt)
%EPSILON 计算越阶函数
% D_i,j      i——j吸引的车流量
%Dt_i,j      i——j临界车流量
m=Dt(i,j)-D(i,j);
if m>0
    f=1;
else
    f=0;
end
f;
end

