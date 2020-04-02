function res2 = mtimes(A,x1)
dt=0.002;                             %数据采样间隔                        
nt=1001;                              %采样点数  
dx=0.005;
x=0*dx:dx:59*dx;                          %偏移距
pmax=20;                            %最大斜率或曲率
p=0:pmax/100:pmax;                    %扫描的斜率或者曲率
np=length(p);
nx=length(x);
if A.adjoint == 0                      %A*x
    x2=reshape(x1,nt,np);
    ress=invfwd_tx_sstackn_linear(x2,dt,p,x);  %时空域逆Radon变换
    res2=reshape(ress,nt*nx,1);
else                                   %At*x
    x2=reshape(x1,nt,nx);
    ress=fwd_tx_sstackn_linear(x2,dt,p,x);     %时空域正Radon变换
    res2=reshape(ress,nt*np,1);
end