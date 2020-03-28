close all; clear all; 
clc;
%%
load mianbo1.mat -mat
nsg = mianbo;
signal=nsg(:,5);
N=length(signal);
dt=0.002;
fs=1/dt;  %采样频率和数据点数
n=0:N-1;   t=n/fs;   %时间序列

figure;
y=fft(signal,N);    
mag=abs(y);     mmag=mag./max(mag);    f=n*fs/N;    % f=(0:N-1)*fs/dt*N;  
plot(f(1:N/2),mmag(1:N/2));   xlabel('Frequency/Hz');ylabel('Amplitude');
%%
[S,F,T,~]=spectrogram(signal,256,200,256,1/dt);
P=abs(S);
figure
pcolor(T,F,P),shading interp;colorbar;
xlabel('Time (Seconds)'); ylabel('Frequency/Hz');
title('Raw data');set(gcf, 'Renderer', 'ZBuffer');
%% 滤波
fre1=25/(1/dt/2);
[B,A]=butter(8,fre1,'low');
mianbonew=filtfilt(B,A,signal);
yy=(1:1001)*dt;
figure;plot(yy,signal);
figure;plot(yy,mianbonew);


%%
load qumianbo1.mat -mat
nsg1 = qumianbo;
signal1=nsg1(:,5);
N1=length(signal1);
dt=0.002;

figure;
fs=1/dt;  
n=0:N1-1;t=n/fs;   
y1=fft(signal1,N1);   
mmag1=abs(y1)./max(mag1);
f1=n*fs/N1;  
plot(f1(1:N1/2),mmag1(1:N1/2));   
xlabel('Frequency/Hz');ylabel('Amplitude');


%% STFT 
[S1,F1,T1,P1]=spectrogram(signal1,256,200,256,1/dt);%绘制原始信号时频图
P1=abs(S1);
figure
pcolor(T1,F1,P1),shading interp;colorbar;
xlabel('Time (Seconds)'); ylabel('Frequency/Hz');
title('原始信号');set(gcf, 'Renderer', 'ZBuffer');
%% Filter
fre1=200/(1/dt/2);
[B,A]=butter(8,fre1,'low');
qumianbonew=filtfilt(B,A,signal1);
figure;plot(yy,signal1);
figure;plot(yy,qumianbonew);
