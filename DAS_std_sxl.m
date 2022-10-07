%% read data
clear
load('50ns_10khz_pzt200hz_3v.mat')
%% reshape data
nd=1000; % num of traces
N=6250;%length of traces
L2=62500;
data_A=zeros(nd,N);
for i=1:nd
data_A(i,:) = A((1+31000+(i-1)*L2):((N+(i-1)*L2)+31000));
end
figure;plot(data_A(100:120,:)');title('raw data after seperation')
%%
fs=625e6; %sampling frequency
f=(0:N-1)*fs/(N-1);
fy=f-fs/2;
pf1=70e6;
pf2=90e6;
pyyfs1=zeros(1,fix(N*(1/2-pf2/fs)));
pyyfs2=ones(1,fix(N*(1/2-pf1/fs))-fix(N*(1/2-pf2/fs)));
pyyfs3=zeros(1,fix(N*(1/2+pf1/fs))-fix(N*(1/2-pf1/fs))+1);
pyyfs4=ones(1,fix(N*(1/2+pf2/fs))-fix(N*(1/2+pf1/fs)));
pyyfs5=zeros(1,N-fix(N*(1/2+pf2/fs))-1);
pyyfsm=[pyyfs1 pyyfs2 pyyfs3 pyyfs4 pyyfs5];
figure;plot(fy,pyyfsm);title('filter shape');
%% filter data
data_A_filter_freq=zeros(nd,N);
data_A_filter_time=zeros(nd,N);

parfor i=1:nd
data_A_filter_freq(i,:)=fftshift(fft(data_A(i,:),N)).*pyyfsm;
data_A_filter_time(i,:)=ifft(ifftshift(data_A_filter_freq(i,:)),N,'symmetric');

end

figure;plot(fy,abs(fftshift(fft(data_A(100,:),N)))); hold on; plot(fy,abs(data_A_filter_freq(100,:)),'r');title('frequency shape before and after filtering');
figure;plot(data_A(100,:));hold on;plot((data_A_filter_time(100,:)),'r');title('time domain before and after filtering');
%%
k=1:N;
phi=7.8125; %phi=fs/frequency_shift
M_sin=sin(2*pi/phi*(k-1)); 
M_cos=cos(2*pi/phi*(k-1));
data_sin=zeros(nd,N);
data_cos=zeros(nd,N);
parfor i=1:nd
data_sin(i,:)=data_A_filter_time(i,:).*M_sin;
data_cos(i,:)=data_A_filter_time(i,:).*M_cos;
end
%%
f1 = 0e6;  f2=10e6;  N2=N;
yyfs11=zeros(1,fix(N2*(1/2-f2/fs))-1+1);
yyfs22=ones(1,fix(N2*(1/2-f1/fs))-fix(N2*(1/2-f2/fs))-1+1);
yyfs33=zeros(1,fix(N2*(1/2+f1/fs))-fix(N2*(1/2-f1/fs))-1+1);
yyfs44=ones(1,fix(N2*(1/2+f2/fs))-fix(N2*(1/2+f1/fs))-1+1);
yyfs55=zeros(1,N2-fix(N2*(1/2+f2/fs))-1+1-1);
yyfsm2=[yyfs11 yyfs22 1 yyfs33 yyfs44 yyfs55];   % 851111
%figure(5);plot(fy,yyfsm2);

%%
I_f=zeros(nd,N);
Q_f=zeros(nd,N);
parfor i=1:nd
yy_sin=fftshift(fft(data_sin(i,:),N)).*yyfsm2;
yy_cos=fftshift(fft(data_cos(i,:),N)).*yyfsm2;
I_f(i,:) = ifft(ifftshift((yy_sin)),N,'symmetric');
Q_f(i,:) = ifft(ifftshift((yy_cos)),N,'symmetric');
end

%%
signal=I_f-1i*Q_f;
clear data_sin data_cos
%%
Amp=zeros(nd,N);
parfor i=1:nd
Amp(i,:)=abs(signal(i,:));
end 
figure;plot(Amp(100:120,:)');title('Amp');
figure;surf((Amp'),'EdgeColor','None');
%% phase
phase=zeros(nd,N);
for i=1:nd-1
     phase(i,:)=angle(signal(i,:));
end
%% amp_diff
step=6;
Amp_diff=zeros(nd-step,N);
for i=1:(nd-step)
    ip=i+step;
    Amp_diff(i,:)=Amp(i,:)-Amp(i+step,:);
end
figure;plot(Amp_diff');title('Amp_diff');
%% phase_diff
phase_unwrap=zeros(nd-1,N);
for i=1:nd-1
    phase_unwrap(i,:)=unwrap(phase(i,:));
end
%figure;plot(phase_unwrap(1:4,:)');
phase_diff=zeros(nd-1,N);

for i=55:6250
    phase_diff(:,i)=phase_unwrap(:,i)-phase_unwrap(:,i-54);
end
%figure;surf(unwrap(phase_diff),'EdgeColor','None');
figure;plot(unwrap(phase_diff(:,4051)));