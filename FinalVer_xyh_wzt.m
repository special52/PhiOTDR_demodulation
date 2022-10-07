close all;
clear;%%综合
tic;
%% 读取数据
load('0.mat');%25000（点数）*2500（traces数目）
original_data=s';%转置
data_row = size(original_data,1);%读取列数
data_column = size(original_data,2);%读取行数
data=zeros(data_row,data_column);
load('noise.mat');
noise=s';
B=mean(noise,1);
original_data=double(original_data);
for i=1:data_row
    data(i,:)=original_data(i,:)-B;%减底噪 这一步的意义在于将traces图移植0坐标
end
clear s B original_data noise
%% Reshape in traces
row_start = 1;
row_end = 500;%提取traces数目
column_start = 1; %传感光纤起点 534
column_end = data_column; %传感光纤终点 21062
tempData = data(row_start:row_end,column_start:column_end);%形成背散拍频信号强度矩阵
[trace_N, pos_N]=size(tempData);%
t_x = 1:1:pos_N;
clear data 
%% Beat Signal
figure;
v_daq=1*10^9;%数据采集卡的采集频率1Ga/s
t_daq=1/v_daq;
index_mmf=1.448;%光纤折射率
v_optfiber=3*10^5/index_mmf;%光在光纤中的传播速度
for i=1:10  %10是画的traces条数
    plot(0.5*t_daq*v_optfiber*t_x,tempData(i,:),'Color',[rand(),rand(),rand()]);%两个列向量，t_x是点数，tempData(i,:)是traces的强度
    hold on;
end
% 画出所有traces 没必要！增加了运行时间
% for i=1:trace_N
%     plot(0.5*t_daq*v_optfiber*t_x,tempData(i,:),'Color',[rand(),rand(),rand()]);%两个列向量，t_x是点数，tempData(i,:)是traces的强度
%     hold on;
% end
xlabel('Distance(km)');
ylabel('Amplitude(a.u.)');
%%构建滤波器
%% 滤波器
fs=1e9; %数据采集卡的采样频率
N=data_column;%一条traces的点数
f=(0:N-1)*fs/(N-1);
fy=f-fs/2;
pf1=195e6;
pf2=205e6;
pyyfs1=zeros(1,fix(N*(1/2-pf2/fs)));
pyyfs2=ones(1,fix(N*(1/2-pf1/fs))-fix(N*(1/2-pf2/fs)));
pyyfs3=zeros(1,fix(N*(1/2+pf1/fs))-fix(N*(1/2-pf1/fs))+1);
pyyfs4=ones(1,fix(N*(1/2+pf2/fs))-fix(N*(1/2+pf1/fs)));
pyyfs5=zeros(1,N-fix(N*(1/2+pf2/fs))-1);
pyyfsm=[pyyfs1 pyyfs2 pyyfs3 pyyfs4 pyyfs5];
% figure;plot(fy,pyyfsm);title('filter shape');构建滤波器
clear pyyfs1 pyyfs2 pyyfs3 pyyfs4 pyyfs5
%% filter data
nd=row_end;
data_filter_freq=zeros(nd,N);
data_filter_time=zeros(nd,N);

parfor i=1:nd
data_filter_freq(i,:)=fftshift(fft(tempData(i,:),N)).*pyyfsm;
data_filter_time(i,:)=ifft(ifftshift(data_filter_freq(i,:)),N,'symmetric');
%data_filter_time(i,:)就是后面的tempData
end
figure;
plot(fy,abs(fftshift(fft(tempData(100,:),N))));hold on; 
plot(fy,abs(data_filter_freq(100,:)),'r');
title('frequency shape before and after filtering');
figure;
plot(tempData(100,:));hold on;
plot((data_filter_time(100,:)),'r');title('time domain before and after filtering');
tempData=data_filter_time;
clear data_filter_time;
%% FFT
% 快速傅里叶变换的幅值
y0=abs(fft(tempData(1,:))); 
% N=length(tempData(1,:));跟上面的N大小一样
f0=(0:N-1)*v_daq/N; % 将时间横坐标转换为频率
% 调整0频位置
f1=(0:N-1)*v_daq/N-v_daq/2;
y1=abs(fftshift(fft(tempData(1,:))));
% 用FFT对信号做频谱分析，只需考察0~Nyquist频率范围内的幅频特性
f2=(N/2:N-1)*v_daq/N-v_daq/2; % 频率范围0~25Hz
y2=2*y1(N/2:N-1)/N; % 幅值修正得到真实幅值
figure;
% subplot(2,2,1);plot(t,tempData(1,:));title('Original signal');xlabel('Time');ylabel('Amplitude');%画出traces图
subplot(2,2,2);plot(f0,y0);title('y0=abs(fft(x))');xlabel('Frequency');ylabel('Amplitude');
subplot(2,2,3);plot(f1,y1);title('y1=abs(fftshift(fft(x)))');xlabel('Frequency');ylabel('Amplitude');
subplot(2,2,4);plot(f2,y2);title('y2=2*y1(N/2:N-1)/N');xlabel('Frequency');ylabel('Amplitude');

%% Local signal - Make 2*pi*delt(f)*t
t=t_x*1/v_daq; % 形成时间轴矩阵
frequency_shift=2*10^8;%这玩意是个频移！！！！！不是光速
x = 2*pi*frequency_shift*t_daq*t_x; %phi=fs/frequency_shift=1e9/
refi=cos(x);
refq=sin(x);
%% IQ demodulation
% Filter defination滤波器
b = fir1(200,0.0002);%滤波器的具体参数 待会细看
% figure(10);
freqz(b,1);
ri = zeros(trace_N,pos_N);
rq = zeros(trace_N,pos_N);
% Combination of Beat signal and Local signal震动信号和本地信号的组合
for i=1:trace_N
    tempi = tempData(i,:).*refi;
    tempq = tempData(i,:).*refq;
    ri(i,:) = filter(b,1,tempi);
    rq(i,:) = filter(b,1,tempq);
end
% Phase orthogonal demodulation相位正交解调
aresult = zeros(trace_N, pos_N);
presult = zeros(trace_N, pos_N);
new_phase = zeros(trace_N, pos_N);
phase = presult*0;
for j=1:trace_N
    for i=1:column_end-column_start+1
        aresult(j,i) = sqrt(ri(j,i)^2+rq(j,i)^2);
        presult(j,i) = atan2(rq(j,i),ri(j,i));
    end
    phase(j,:) = unwrap(presult(j,:));%第一次解卷绕（以行为组，列与列解）
end
%% 强度！
% Amplitude Demodulation强度解调
figure;
for i=1:trace_N
    plot(0.5*t_daq*v_optfiber*t_x, aresult(i,:),'Color',[rand(),rand(),rand()]);
    hold on;
end
xlabel('Distance(km)');ylabel('Amplitude(a.u.)');title('Amplitude viration along fiber with time')
% Amplitude positioning强度定位
varonce= std(aresult);% S = std(A)是一个列为随机变量且行为观测值的矩阵，则 S 是一个包含与每列对应的标准差的行向量。
figure;plot(0.5*t_daq*v_optfiber*t_x,varonce/max(varonce));
figure;plot(0.5*t_daq*v_optfiber*t_x,varonce-mean(varonce));%将曲线下移至x轴
xlabel('Distance(km)');ylabel('Amplitude rms');title('Amplitude positioning');
% pulse64=(varonce-mean(varonce))';%取出脉宽为64ns的强度定位图
figure; surf((aresult),'EdgeColor','None');  %画出三维度强度图
[max_var,max_index] = max(varonce);
pzt_distance=max_index*0.5*t_daq*v_optfiber;
%% 尝试另外一个方法的强度解调 通过做差 来得到差值最大的区域从而确定扰动位置 （误差step)
step=2;
Amp_diff=zeros(nd-step,N);
for i=1:(nd-step)
    ip=i+step;
    Amp_diff(i,:)=aresult(i,:)-aresult(i+step,:);
    plot(0.5*t_daq*v_optfiber*t_x, Amp_diff(i,:),'Color',[rand(),rand(),rand()]);
    hold on
end
figure;surf((Amp_diff),'EdgeColor','None');
figure;plot(Amp_diff');title('Amp_diff');
clear Amp_diff
%% 相位
% phase demodulation
column_start =  1; %传感光纤起点 534
column_end = data_column; %传感光纤终点 21062
figure;plot(presult');%相位信息
figure;plot(phase');%解卷绕01相位信息
%消除相位干扰

% for i = 1:column_end-column_start+1% 去掉头尾白噪声处
%         new_phase(:,i) = unwrap(phase(:,i));%第二次解卷绕，以列为组，行与行解
% end
%  figure;plot(new_phase');
phase_diff=zeros(nd,N);
for i=80:25000
    phase_diff(:,i)=phase(:,i)-phase(:,i-79);
    
end
figure;plot(phase_diff');
% figure;plot(phase_point');
% figure;surf(unwrap(phase_diff),'EdgeColor','None');
phase_point=(phase_diff(:,11150));
% figure;plot(phase_point);
F_FFT_20220719(phase_point,25000);
toc;



%%