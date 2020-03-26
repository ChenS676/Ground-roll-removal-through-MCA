function b=invfwd_tx_sstackn_hyper(m,dt,p,x)  
%对原始数据d做时空域逆Radon变换
%输入参数：
%d:二维数据  dt:采样间隔  x:偏移距  p:扫描斜率或者曲率
%输出参数：
%m:Radon域结果
nt=size(m,1);
np=length(p);
nx=length(x);
tz=0*dt:dt:(nt-1)*dt;                       %时间点
data=zeros(nt,nx);
for i=1:np
    for j=1:nt
        for k=1:nx
           %t=tz(j)+p(i)*x(k);                   %叠加路径为直线
           t=sqrt(tz(j)^2+p(i)^2*x(k)^2);      %叠加路径为双曲线
           %t=tz(j)+p(i).^2*x(k).^2;            %叠加路径为抛物线
          tc=t/dt;
          itc=floor(tc);
          it=itc+1;
          fra=tc-itc;
          if it>0&&it<nt
                 data(it,k)=data(it,k)+(1.-fra)*m(j,i);
                 data(it+1,k)=data(it+1,k)+fra*m(j,i);
          end
        end
    end
end
b=data;             