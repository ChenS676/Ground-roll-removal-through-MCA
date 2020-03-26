%%%%%去谐波单道数据处理程序%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; 
clc;
% Add some utility directories to the path
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt  = 0.002; %采样时间间隔，由输入文件决定
N   = 1001;  %每道数据的长度，由输入文件决定
nt=N;
nx=1;  %实际截取道数
dx=0.005;
M=120;       %文件的总道数，由输入文件决定
t=(0:N-1).*dt;
t=t'; 
n=nt*nx;
nsg=zeros(nt,nx);
qumianbo=zeros(nt,nx);
mianbo=zeros(nt,nx);
yuanshishuju=zeros(nt,nx);
% for i=1:M
%     trace_head = fread(fidin, [240,1], '*uchar');%%%%%读240字节道头
% %    
%      nsg(:,i)=fread(fidin,[N,1],'float');
%  end
% fclose(fidin);
% 
% fnq  = 1/dt/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The value of Q is small
pmax1=5;
p1=0:pmax1/100:pmax1;
x1=0*dx:dx:59*dx;
dictWave1 = struct('p1', p1,'x1',x1);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% The value of Q is large
pmax2=20;
p2=0:pmax2/100:pmax2;
x2=0*dx:dx:59*dx;
dictWave2 = struct('p2', p2,'x2',x2);
% Bodywave/groundroll separation
% Call the BCR to sovle the optimization problem   
CQ1 =0.85;%Q1值的小波变换阈值权重
CQ2 =1.5;%Q2值的小波变换阈值权重
Cweight     = struct('CQ1', CQ1, 'CQ2', CQ2);
 thdtype     = 'hard';%阈值类型
itermax 	= 30;
expdecrease	= 1;      %阈值按照线性减小
lambdastop=5e-1;
sigma      = 5e-3;
display		= 0;

fid = fopen('E:\Radon变换字典MCA应用程序\MCAsuppressgroundroll\HarmNoisRemv_去谐波_最终程序\HarmNoisRemv\实验原数据\newgroundmodel.dat','r');%输入文件
[ddp,count]=fread(fid,[1001,120],'float');
% fidout1 = fopen('E:\MCA实际数据分离结果\典型模拟数据newgroundmodel的分离结果\去面波后1.sgy','w');%输出文件
% fidout2 = fopen('E:\MCA实际数据分离结果\典型模拟数据newgroundmodel的分离结果\去掉面波1.sgy','w');%输出文件
% fidout3 = fopen('E:\MCA实际数据分离结果\典型模拟数据newgroundmodel的分离结果\原信号1.sgy','w');%输出文件
%  noprocessing=zeros(1,M);
     
nsg=ddp(:,65);
   parts = MCA_Bcr(nsg,dictWave1,dictWave2,Cweight,thdtype,itermax,expdecrease,lambdastop,sigma,display,dt);
   nsg1=reshape(nsg,n,1);
   parts(:,1)=nsg1-parts(:,2);
   mianbo=reshape(parts(:,2),nt,nx);
   qumianbo=reshape(parts(:,1),nt,nx);
   yuanshishuju=nsg;
%    fwrite(fidout1,parts(:,1),'float');
%    fwrite(fidout2,parts(:,2),'float');
%    fwrite(fidout3,nsg(:,i),'float');
fclose(fid);
% fclose(fidout1);
% fclose(fidout2);
% fclose(fidout3);
% noprocessing
%  这里取阶子带置零
yy=(1:1001)*dt;
toc
figure;plot(yy,yuanshishuju);
figure;plot(yy,mianbo);
figure;plot(yy,qumianbo);
