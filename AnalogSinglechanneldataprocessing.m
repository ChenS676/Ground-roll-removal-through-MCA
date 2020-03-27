close all; clear all; clc;
%% load data frequency analysis and time-frequency analysis
load mianbo1.mat -mat
gr = mianbo;

%% Fourier analysis of arbitrary single shot 1d signal
signal=gr(:,5); % signal
N=length(signal); 
dt=0.002;     % sampling rate
figure;
fs=1/dt;      % sampling frequency
n=0:N-1;t=n/fs;      % time_coordinate
y=fft(signal,N);   
mag=abs(y);  
mmag=mag./max(mag);
f=n*fs/N;
plot(f(1:N/2),mmag(1:N/2));   
xlabel('Frequency/Hz');ylabel('Amplitude');

%% Short Time Fourier Transformation of single signal
[S,F,T,~]=spectrogram(signal,256,200,256,1/dt);
P=abs(S);
figure
pcolor(T,F,P),shading interp;colorbar;
xlabel('Time (Seconds)'); ylabel('Frequency/Hz');
title('Raw data');set(gcf, 'Renderer', 'ZBuffer');


%% Filter in single shot signal
f_C=25/(1/dt/2);
[B,A]=butter(8,f_C,'low');
mianbonew=filtfilt(B,A,signal);
yy=(1:1001)*dt;
figure;plot(yy,signal);
figure;plot(yy,mianbonew);


%% Fourier Analysis of single shot of body wave 
load qumianbo1.mat -mat
bw = qumianbo;
i=5;
signal1=bw(:,i);
N1=length(signal1);
dt=0.002;
figure;
fs=1/dt;  
n=0:N1-1;  t=n/fs;  
y1=fft(signal1,N1);   
mag1=abs(y1);   
mmag1=mag1./max(mag1);
f1=n*fs/N1;    %频率序列
plot(f1(1:N1/2),mmag1(1:N1/2));   %绘出随频率变化的振幅
xlabel('Frequency/Hz');ylabel('Amplitude');

%% STFT Analysis of single shot of body wave
[S1,F1,T1,P1]=spectrogram(signal1,256,200,256,1/dt);
P1=abs(S1);
figure
pcolor(T1,F1,P1),shading interp;colorbar;
xlabel('Time (Seconds)'); ylabel('Frequency/Hz');
title('原始信号');set(gcf, 'Renderer', 'ZBuffer');


%% Filter of signal
f_C=200/(1/dt/2);
[B,A]=butter(8,f_C,'low');
qumianbonew=filtfilt(B,A,signal1);
figure;plot(yy,signal1);
figure;plot(yy,qumianbonew);
