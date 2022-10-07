function [xFreq,yOutput]=F_FFT_20220719(yInput,v_daq)

%Beat signal
% v_daq=1*10^9; % 采样率1G
% t=t_x*1/v_daq; % 形成时间轴矩阵
% f1=15;
% f2=20;
% x=sin(2*pi*f1*t)+sin(2*pi*f2*t); % 原始信号

% 快速傅里叶变换的幅值
y0=abs(fft(yInput)); 
N=length(yInput);
f0=(0:N-1)*v_daq/N; % 将时间横坐标转换为频率

% 调整0频位置
f1=(0:N-1)*v_daq/N-v_daq/2;
y1=abs(fftshift(fft(yInput)));

% 用FFT对信号做频谱分析，只需考察0~Nyquist频率范围内的幅频特性
f2=(N/2:N-1)*v_daq/N-v_daq/2; % 频率范围0~25Hz
y2=2*y1(N/2:N-1)/N; % 幅值修正得到真实幅值

% yOutput=y2;
% xFreq=f2;

figure(10);
y2_1=(y2-min(y2))/(max(y2)-min(y2));
plot(f2,y2_1);
yOutput=y2_1;
xFreq=f2;


end