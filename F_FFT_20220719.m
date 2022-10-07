function [xFreq,yOutput]=F_FFT_20220719(yInput,v_daq)

%Beat signal
% v_daq=1*10^9; % ������1G
% t=t_x*1/v_daq; % �γ�ʱ�������
% f1=15;
% f2=20;
% x=sin(2*pi*f1*t)+sin(2*pi*f2*t); % ԭʼ�ź�

% ���ٸ���Ҷ�任�ķ�ֵ
y0=abs(fft(yInput)); 
N=length(yInput);
f0=(0:N-1)*v_daq/N; % ��ʱ�������ת��ΪƵ��

% ����0Ƶλ��
f1=(0:N-1)*v_daq/N-v_daq/2;
y1=abs(fftshift(fft(yInput)));

% ��FFT���ź���Ƶ�׷�����ֻ�迼��0~NyquistƵ�ʷ�Χ�ڵķ�Ƶ����
f2=(N/2:N-1)*v_daq/N-v_daq/2; % Ƶ�ʷ�Χ0~25Hz
y2=2*y1(N/2:N-1)/N; % ��ֵ�����õ���ʵ��ֵ

% yOutput=y2;
% xFreq=f2;

figure(10);
y2_1=(y2-min(y2))/(max(y2)-min(y2));
plot(f2,y2_1);
yOutput=y2_1;
xFreq=f2;


end