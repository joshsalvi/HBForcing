%% WHITE NOISE

warning off
clear q_stim psd_q1_stim psd_q2_stim qfft_stim


for j = 1:a
for i = 1:(logdata.data(1+ramp_loc1-1,iter_loc))
    N=1;
    clear z;
    z = smooth(Xd_pulse(:,i,j),length(Xd_pulse(:,i,j))/N);
    Xd_pulse_center(:,i,j)=Xd_pulse(:,i,j)-z;
end
end


time = length(Xd_pulse_center);    
T=1/Fs;
NFFT = 2^3*2^nextpow2(length(Xd_pulse_center));

f = (Fs/2)*linspace(0,1,NFFT/2+1);
effective_delta_f = Fs/NFFT;
actual_delta_f = Fs/numel(time);



% STIMULUS, X
for j = 1:a
for i = 1:(logdata.data(1+ramp_loc1-1,iter_loc))
    
    
h_stim=psd(spectrum.periodogram('rectangular'),Xd_pulse_center(:,i,j),'Fs',Fs,'NFFT',NFFT);

L = length(Xd_pulse_center(:,i,j));

clear y y0 wn

%{
hfft  = fft(Xd_pulse_center(Xd_pulse_center(:,i,j)~=0,p,i,j),NFFT)/L;
qfft_stim(:,2)=2*hfft(1:NFFT/2+1);
qfft_stim(:,1) = Fs/2*linspace(0,1,NFFT/2+1);
%}
y0 = Xd_pulse_center(:,i,j);     %windowing required for short signals
%wn = hann(length(y0));
wn = rectwin(length(y0));
for i2 = 1:length(wn)
    y(i2) = wn(i2)*y0(i2);
end
if isempty(y0)
    break;
end
hfft  = fft(y,NFFT)/L;
qfft_stim(:,2)=hfft(1:NFFT/2+1);
qfft_stim(:,1) = Fs/2*linspace(0,1,NFFT/2+1);

clear q_stim;
q_stim(:,1)=h_stim.Frequencies;
q_stim(:,2)=h_stim.Data;


psd_q1_stim(:,i,j)=q_stim(:,1);    %frequency PRE
psd_q2_stim(:,i,j)=q_stim(:,2);    %psd PRE

fft_q1_stim(:,i,j)=qfft_stim(:,1);    %frequency PRE
fft_q2_stim(:,i,j)=qfft_stim(:,2);    %psd PRE


end
end

