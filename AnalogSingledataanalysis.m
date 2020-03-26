close all; clear all; 
clc;
load mianbo1.mat -mat
nsg = mianbo;
i=5;
signal=nsg(:,i);
N=length(signal);
dt=0.002;
figure;
fs=1/dt;  %采样频率和数据点数
n=0:N-1;t=n/fs;   %时间序列
y=fft(signal,N);    %对信号进行快速Fourier变换
mag=abs(y);     %求得Fourier变换后的振幅
mmag=max(mag);
mmag=mag./mmag;%%%%%归一化
f=n*fs/N;    %频率序列
plot(f(1:N/2),mmag(1:N/2));   %绘出随频率变化的振幅
xlabel('Frequency/Hz');ylabel('Amplitude');

[S,F,T,P]=spectrogram(signal,256,200,256,1/dt);%绘制原始信号时频图
P=abs(S);
figure
pcolor(T,F,P),shading interp;colorbar;
xlabel('Time (Seconds)'); ylabel('Frequency/Hz');
title('原始信号');set(gcf, 'Renderer', 'ZBuffer');
%%%%%%%%%%%%%%%%滤波
fre1=25/(1/dt/2);
[B,A]=butter(8,fre1,'low');
mianbonew=filtfilt(B,A,signal);
yy=(1:1001)*dt;
figure;plot(yy,signal);
figure;plot(yy,mianbonew);
load qumianbo1.mat -mat
nsg1 = qumianbo;
i=5;
signal1=nsg1(:,i);
N1=length(signal1);
dt=0.002;
figure;
fs=1/dt;  %采样频率和数据点数
n=0:N1-1;t=n/fs;   %时间序列
y1=fft(signal1,N1);    %对信号进行快速Fourier变换
mag1=abs(y1);     %求得Fourier变换后的振幅
mmag1=max(mag1);
mmag1=mag1./mmag1;%%%%%归一化
f1=n*fs/N1;    %频率序列
plot(f1(1:N1/2),mmag1(1:N1/2));   %绘出随频率变化的振幅
xlabel('Frequency/Hz');ylabel('Amplitude');

[S1,F1,T1,P1]=spectrogram(signal1,256,200,256,1/dt);%绘制原始信号时频图
P1=abs(S1);
figure
pcolor(T1,F1,P1),shading interp;colorbar;
xlabel('Time (Seconds)'); ylabel('Frequency/Hz');
title('原始信号');set(gcf, 'Renderer', 'ZBuffer');
%%%%%%%%%%%%%%%%滤波
fre1=200/(1/dt/2);
[B,A]=butter(8,fre1,'low');
qumianbonew=filtfilt(B,A,signal1);
figure;plot(yy,signal1);
figure;plot(yy,qumianbonew);
