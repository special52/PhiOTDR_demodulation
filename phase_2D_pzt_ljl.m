%%
%----demodulated phase named as 'new_phase'---%
%----demoddulated amplitude named as 'aresult'---%

z=0.5*t_daq*v_optfiber*t_x;

%%
%------Amplitude positioning------%
varonce= std(aresult);

figure(12)
plot(z,varonce/max(varonce));
xlabel('Distance(km)');ylabel('Amplitude rms');title('Amplitude positioning');

%%
[max_var,max_index] = max(varonce(1:length(varonce)));
pos_signal=max_index*0.5*t_daq*v_optfiber*t_x;

%%
%-----Phase positioning---------%

differetial_phase_1=zeros(trace_N,pos_N);
for ii = 1:pos_N-250
%    differetial_phase_1(:,ii)= new_phase1(:,ii)- new_phase1(:,ii+2); % without the second low pass filter 
    differetial_phase_1(:,ii)= new_phase2(:,ii)- new_phase2(:,ii+250); % with the second low pass filter
end

figure(13)
plot(z, differetial_phase_1);xlabel('Distance(km)');ylabel('Phase difference');title('Phase positioning');

figure(14)
plot(t_x, differetial_phase_1);xlabel('Point');ylabel('Phase difference');title('Phase positioning');