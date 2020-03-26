function m=cgnr_radon_hyper(d,dt,x,p)
%时空域高分辨率Radon变换
%输入参数：
%d:二维数据  dt:采样间隔  x:偏移距  p:扫描斜率或者曲率
%输出参数：
%m:Radon域结果
[nt,nx]=size(d);
np=length(p);
M=ones(nt*np,1);  %初始的加权矩阵
size(M)
d1=reshape(d,nt*nx,1);               %将二维数据拉成一维
L=radon_transform1();                %重载乘法，重载算子L
%采用CGLS算法求解大型的矩阵方程
for j=1:5                           %外部迭代次数
m0=zeros(size(M));
s0=d1-L*m0;
r0=L'*s0;
r0=M.*r0;
p0=r0;
pn=M.*p0;
q0=L*pn;
k=1;
while k<=50                     %内部迭代次数
    a0=norm(r0)^2/norm(q0)^2;
    m1=m0+a0*p0;
    s1=s0-a0*q0;
    r1=L'*s1;
    r1=M.*r1;
    b0=norm(r1)^2/norm(r0)^2;
    p1=r1+b0*p0;
    pn=M.*p1;
    q1=L*pn;
    m0=m1;
    s0=s1;
    r0=r1;
    p0=p1;
    q0=q1;
    k=k+1;
end
m0=m0.*M;
%var_m=0.0001*max(abs(m0).^2);
var_m=0.0001;
M=sqrt(var_m+abs(m0).^2);                       %加权矩阵
end
m=reshape(m0,nt,np);