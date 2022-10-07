close all;
clear;
%% Read in data (in traces)
tic;
s=[];
List =dir('D:\毕设准备\2022-09-16 16_19_58\*.mat');
k =length(List)-2;
for i=1:k
    file_name{i}=List(i).name;
    temp=importdata(file_name{i});
    temp=temp';
    s=[s;temp];
end
path_pri = ('Sig.mat');
save(path_pri,'s','-v7.3')

%%

% [a b]=size(s);
% p=127*ones([a b]);
% s=double(s);
% p=double(p);
% data(i,:)=s-mean(s,1);
data_row = size(s,1);
data_column = size(s,2);
data=zeros(data_row,data_column);
s1 = load('noise.mat');
n0 = struct2cell(s1);
n1 = cell2mat(n0);
n1=n1';
% s1=load('noise.mat');
B=mean(n1,1);
s=double(s);
for i=1:data_row
    data(i,:)=s(i,:)-B;%减底噪
end
clear s;
clear s1;
clear B;
clear n1;


%% Reshape in traces

row_start = 1;
row_end = 10000;
column_start = 1;
column_end = data_column;

tempData = data(row_start:row_end,column_start:column_end);%形成背散拍频信号强度矩阵
[trace_N, pos_N]=size(tempData);
t_x = 1:1:pos_N;
clear data;

%% Beat Signal

figure(1);
v_daq=1*10^9;
t_daq=1/v_daq;
index_mmf=1.5;
v_optfiber=3*10^5/index_mmf;

for i=1:trace_N
    plot(0.5*t_daq*v_optfiber*t_x,tempData(i,:),'Color',[rand(),rand(),rand()]);%横坐标对应时间间隔为采集卡采样点间隔的1/2
    hold on;
end
xlabel('Distance(km)');
ylabel('Amplitude(a.u.)');

%% FFT

%Beat signal
% v_daq=0.5*10^9; % 采样率1G
t=t_x*1/v_daq; % 形成时间轴矩阵
% f1=15;
% f2=20;
% x=sin(2*pi*f1*t)+sin(2*pi*f2*t); % 原始信号

% 快速傅里叶变换的幅值
y0=abs(fft(tempData(1,:))); 
N=length(tempData(1,:));
f0=(0:N-1)*v_daq/N; % 将时间横坐标转换为频率

% 调整0频位置
f1=(0:N-1)*v_daq/N-v_daq/2;
y1=abs(fftshift(fft(tempData(1,:))));

% 用FFT对信号做频谱分析，只需考察0~Nyquist频率范围内的幅频特性
f2=(N/2:N-1)*v_daq/N-v_daq/2; % 频率范围0~25Hz
y2=2*y1(N/2:N-1)/N; % 幅值修正得到真实幅值

figure(2);
subplot(2,2,1);
plot(t,tempData(1,:));
title('Original signal');
xlabel('Time');
ylabel('Amplitude');

subplot(2,2,2);
plot(f0,y0);
title('y0=abs(fft(x))');
xlabel('Frequency');
ylabel('Amplitude');

subplot(2,2,3);
plot(f1,y1);
title('y1=abs(fftshift(fft(x)))');
xlabel('Frequency');
ylabel('Amplitude');

subplot(2,2,4);
plot(f2,y2);
title('y2=2*y1(N/2:N-1)/N');
xlabel('Frequency');
ylabel('Amplitude');

%% Local signal - Make 2*pi*delt(f)*t
% v_daq=1*10^9; % 采样率1G
t=t_x*1/v_daq; % 形成时间轴矩阵
v_shift=2*10^8;
x = 2*pi*v_shift*t_x*t_daq;
yi=cos(x);
yq=sin(x);
refi = yi;
refq = yq;
% %调整坐标
y_0=abs(fft(refi)); 
N=length(refi);
f_0=(0:N-1)*v_daq/N; % 将时间横坐标转换为频率

% 调整0频位置
f_1=(0:N-1)*v_daq/N-v_daq/2;
y_1=abs(fftshift(fft(refi)));

% x = 0:2*pi/5:2*pi/5*4;
% yi = cos(x);
% yq = sin(x);
% 
% refi = zeros(1,pos_N);
% refq = zeros(1,pos_N);
% 
% for i=1:5:column_end - column_start - 3
%     refi(i:i+4) = yi;
%     refq(i:i+4) = yq;
% end
% figure(3);
% subplot(3,1,1);
% plot(t,refi);
% title('Original signal');
% xlabel('Time');
% ylabel('Amplitude');
% 
% subplot(3,1,2);
% plot(f_0,y_0);
% title('y0=abs(fft(refi))');
% xlabel('Frequency');
% ylabel('Amplitude');
% 
% subplot(3,1,3);
% plot(f_1,y_1);
% title('y1=abs(fftshift(fft(refi)))');
% xlabel('Frequency');
% ylabel('Amplitude');
% % plot(abs(fftshift(fft(refi))));
% % title('frequency of local signal');




%% IQ demodulation

% Filter defination
% b = fir1(10,0.002);
b = fir1(100,0.0002);
figure(4);
freqz(b,1);

ri = zeros(trace_N,pos_N);
rq = zeros(trace_N,pos_N);

% refi = data_0hz(1,501:25500);
% refq = data_0hz(1,504:25503);

% Combination of Beat signal and Local signal
for i=1:trace_N
    tempi = tempData(i,:).*refi;
    tempq = tempData(i,:).*refq;
    ri(i,:) = filter(b,1,tempi);
    rq(i,:) = filter(b,1,tempq);
end

% Phase orthogonal demodulation
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

% for i = 1:trace_N
%     for j = 1:column_end-column_start+1
%         phase(:,j) = unwrap(phase(:,j));
%     end
% end
clear rq;
clear ri;
clear aresult;
clear presult

%%
% Amplitude Demodulation
figure(5);
for i=1:trace_N
    plot(0.5*t_daq*v_optfiber*t_x, aresult(i,:),'Color',[rand(),rand(),rand()]);
    hold on;
end
xlabel('Distance(km)');
ylabel('Amplitude(a.u.)');

title('Amplitude viration along fiber with time')

%%
% phase demodulation
% 
    for i = 1:column_end-column_start+1%头尾白噪声处
        new_phase(:,i) = unwrap(phase(:,i));%第二次解卷绕，以列为组，行与行解
    end

figure(6);
new_phase1 = new_phase*0;
for i = 1:trace_N % trace_N
    new_phase1(i,:) =new_phase(i,:)-new_phase(1,:); % 
    plot(0.5*t_daq*v_optfiber*t_x,new_phase1(i,:),'Color',[rand(),rand(),rand()]);
    hold on;
end
xlabel('Distance(km)');
ylabel('Phase(rad)');

title('phase viration along fiber with time')
clear new_phase;
%%
% phase varition along fiber with time

b1 = fir1(50,0.0001);

figure(7); 
new_phase2 = new_phase1*0;
for i = 2:trace_N
    new_phase2(i,:) = filter(b1,1,new_phase1(i,:));
    plot(0.5*t_daq*v_optfiber*t_x,new_phase2(i,:),'Color',[rand(),rand(),rand()]);
    hold on;
end
xlabel('Distance(km)');
ylabel('Phase difference(rad)');
clear new_phase1;
clear new_phase2;

%% Phase demodulation

% % %调整坐标
% y__0=abs(fft(new_phase(1,:))); 
% N=length(new_phase(1,:));
% f__0=(0:N-1)*v_daq/N; % 将时间横坐标转换为频率
% 
% % 调整0频位置
% f__1=(0:N-1)*v_daq/N-v_daq/2;
% y__1=abs(fftshift(fft(new_phase(1,:))));
% 
% 
% figure(8);
% subplot(3,1,1);
% plot(t,new_phase(1,:));
% title('Original signal');
% xlabel('Time');
% ylabel('Amplitude');
% 
% subplot(3,1,2);
% plot(f__0,y__0);
% title('y0=abs(fft(new_phase(1,:)))');
% xlabel('Frequency');
% ylabel('Amplitude');
% 
% subplot(3,1,3);
% plot(f__1,y__1);
% title('y1=abs(fftshift(fft(new_phase(1,:))))');
% xlabel('Frequency');
% ylabel('Amplitude');

% plot(abs(fftshift(fft(new_phase(1,:)))));
% title('frequency of demodulated phase')
% 
% 
% 
% fprintf('..........End calc............\n');

%% test
div=100;
[traceNum,pointsNum]=size(phase);
phaseTemp=zeros(traceNum,pointsNum-div);
for i=1:pointsNum-div
    phaseTemp(:,i)=phase(:,i+div)-phase(:,i);
end
phaseTempU = phaseTemp*0;
for i=1:traceNum
   phaseTempU(i,:)=unwrap(phaseTemp(i,:)); 
end
phaseTempUnwrapDiff=diff(phaseTempU);
figure (8);
subplot(2,1,1);
imagesc(phaseTempU);
subplot(2,1,2);
imagesc(phaseTempUnwrapDiff);



%%
yVibration=phaseTempU(:,11258:8:11386);%相位提取位置扩展范围（多提取几个光纤位置的相位）

figure (9);
plot(yVibration);
xlabel('Traces');
ylabel('Amplitude');
hold on;
p=(1:1:2000);
q=cos(10000*p);
plot(p,q);
title('Phase Extraction')
clear phaseTempUnwrapDiff;
clear phaseTemp;
clear phaseTempU;
%%
% for j=1:size(yVibration,2)
%     F_FFT_20220719(yVibration(:,j),25000);%需要画出多条频率解调曲线，待解决
%     xlabel('Frequency(Hz)');
%     ylabel('Amplitude(a.u.)');
%     title('Frequency Demodulation')
%     cd('F:\OTDR\2022.9.16 Laser-40mW L-2500m F-25kHz P=64ns G=5dBm\MMF\Frequency\10kHz\2022-09-16 10_45_37') %把当前工作目录切换到指定文件夹
%     pictureName = sprintf( '%d.png',j);
%     saveas(gcf,pictureName);
% %         cd('D:\光纤传感\OTDR\To_zihan\3_Mattlab code\3016_OTDR_Simulation_Theory_0411') %切回原工作目录
%    
%     pause('on')
%     pause(0.1)
%     pause('off')
% end
pointnumber = size(yVibration,2);
K = (1:1:pointnumber)';
for j=1:pointnumber
    F_FFT_20220719(yVibration(:,j),25000);%需要画出多条频率解调曲线，待解决
%     legend(num2str(K))
    hold on;
end


%%sss
% F_FFT_20220719(yVibration,25000);


%% 
% figure (11);
% plot(phaseTempUnwrapDiff(1:50,:)');
toc;