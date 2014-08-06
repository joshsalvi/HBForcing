for i = 1:numfreq
for j = 1:5
for k = 1:3
stim_split(:,i,j,k) = stim(1+(k-1)*20000:1+k*20000,i,1,j);
end
end
end

for i = 1:numfreq
for j = 1:5
for k = 1:3
stim_split_z(:,i,j,k) = smooth(stim_split(:,i,j,k),length(stim_split(:,i,j,k))/1);
stim_split_center(:,i,j,k)=stim_split(:,i,j,k)-stim_split_z(:,i,j,k);
end
end
end

%% FREQUENCY ARRAY

warning off
clear q_stim psd_q1_stim psd_q2_stim qfft_stim

time = length(stim_split_center);    
T=1/Fs;
NFFT = 2^3*2^nextpow2(length(stim_split_center));

f = (Fs/2)*linspace(0,1,NFFT/2+1);
effective_delta_f = Fs/NFFT;
actual_delta_f = Fs/numel(time);



% STIMULUS, X
for j = 1:a
for i = 1:size(stim_split_center,3)
for p = 1:numfreq
    
    
h_stim=psd(spectrum.periodogram('rectangular'),stim_split_center(stim_split_center(1:length(stim_split_center),p,i,j)~=0,p,i,j),'Fs',Fs,'NFFT',NFFT);

L = length(stim_split_center(stim_split_center(1:length(stim_split_center),p,i,j)~=0,p,i,j));

clear y y0 wn

%{
hfft  = fft(stim_split_center(stim_split_center(1:stim_time(p),p,i,j)~=0,p,i,j),NFFT)/L;
qfft_stim(:,2)=2*hfft(1:NFFT/2+1);
qfft_stim(:,1) = Fs/2*linspace(0,1,NFFT/2+1);
%}
y0 = stim_split_center(stim_split_center(1:length(stim_split_center),p,i,j)~=0,p,i,j);     %windowing required for short signals
wn = hann(length(y0));
%wn = rectwin(length(y0));
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

if p>1
    q_stim(length(psd_q1_stim),:)=0;
end

psd_q1_stim(:,p,i,j)=q_stim(:,1);    %frequency PRE
psd_q2_stim(:,p,i,j)=q_stim(:,2);    %psd PRE

q = find((freq_stim(p)-1)<psd_q1_stim(:,p,i,j)&1.05*(freq_stim(p)+1)>psd_q1_stim(:,p,i,j));
q = q-2:1:q+2;
q2 = q(find(psd_q2_stim(q,p,i,j)==max(psd_q2_stim(q,p,i,j))));

q3 = findnearest(qfft_stim(:,1),freq_stim(p));

stim_amp1(p,i,j)=sqrt(psd_q1_stim(q2,p,i,j)*psd_q2_stim(q2,p,i,j));       % Output peak amplitude @ stimulus frequency
stim_pow1(p,i,j)=psd_q2_stim(q2,p,i,j);

stim_angle=(angle(qfft_stim(:,2)));
stim_angle1(p,i,j)=stim_angle(q3);
stim_fft1(p,i,j)=2*qfft_stim(q3,2);     % multiply by 2, two-sided fft


end
end
end
