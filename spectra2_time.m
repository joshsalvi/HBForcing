%% FREQUENCY ARRAY

warning off
clear q_stim psd_q1_stim psd_q2_stim qfft_stim

time = length(stim_center);    
T=1/Fs;
NFFT = 2^3*2^nextpow2(length(stim_center));

f = (Fs/2)*linspace(0,1,NFFT/2+1);
effective_delta_f = Fs/NFFT;
actual_delta_f = Fs/numel(time);



% STIMULUS, X
for j = 1:a
for i = 1:(logdata.data(1+ramp_loc1-1,iter_loc))
for p = 1:numfreq
    
    
h_stim=psd(spectrum.periodogram('rectangular'),stim_center(stim_center(1:stim_time(p),p,i,j)~=0,p,i,j),'Fs',Fs,'NFFT',NFFT);

L = length(stim_center(stim_center(1:stim_time(p),p,i,j)~=0,p,i,j));

clear y y0 wn

%{
hfft  = fft(stim_center(stim_center(1:stim_time(p),p,i,j)~=0,p,i,j),NFFT)/L;
qfft_stim(:,2)=2*hfft(1:NFFT/2+1);
qfft_stim(:,1) = Fs/2*linspace(0,1,NFFT/2+1);
%}
y0 = stim_center(stim_center(1:stim_time(p),p,i,j)~=0,p,i,j);     %windowing required for short signals
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


% FIBER, DELTA
for j = 1:a
for i = 1:(logdata.data(1+ramp_loc1-1,iter_loc))
for p = 1:numfreq
    
    
h_stim=psd(spectrum.periodogram('hann'),Delta_in_center(Delta_in_center(1:stim_time(p),p,i,j)~=0,p,i,j),'Fs',Fs,'NFFT',NFFT);

L = length(Delta_in_center(Delta_in_center(1:stim_time(p),p,i,j)~=0,p,i,j));
%{
hfft  = fft(Delta_in_center(Delta_in_center(1:round(cycles/freq_stim(p)*1e4),p,i,j)~=0,p,i,j),NFFT)/L;
qfft_Delta(:,2)=2*hfft(1:NFFT/2+1);
qfft_Delta(:,1) = Fs/2*linspace(0,1,NFFT/2+1);
%}
clear y y0 wn

y0 = Delta_in_center(Delta_in_center(1:stim_time(p),p,i,j)~=0,p,i,j);
%wn = hann(length(y0));       %windowing required for short signals
wn = rectwin(length(y0));
for i2 = 1:length(wn)
    y(i2) = wn(i2)*y0(i2);
end
if isempty(y0)
    break;
end
hfft  = fft(y,NFFT)/L;
qfft_Delta(:,2)=hfft(1:NFFT/2+1);
qfft_Delta(:,1) = Fs/2*linspace(0,1,NFFT/2+1);


clear q_stim;
q_stim(:,1)=h_stim.Frequencies;
q_stim(:,2)=h_stim.Data;

if p>1
    q_stim(length(psd_q1_stim),:)=0;
end

psd_q1_stim(:,p,i,j)=q_stim(:,1);    %frequency PRE
psd_q2_stim(:,p,i,j)=q_stim(:,2);    %psd PRE

q = find((freq_stim(p)-1)<psd_q1_stim(:,p,i,j)&1.05*(freq_stim(p)+1)>psd_q1_stim(:,p,i,j));
%q = q-2:1:q+2;
q2 = q(find(psd_q2_stim(q,p,i,j)==max(psd_q2_stim(q,p,i,j))));
q3 = findnearest(qfft_Delta(:,1),freq_stim(p));

Delta_in_amp1(p,i,j)=sqrt(psd_q1_stim(q2,p,i,j)*psd_q2_stim(q2,p,i,j));       % Output peak amplitude @ stimulus frequency
Delta_in_pow1(p,i,j)=psd_q2_stim(q2,p,i,j);

Delta_in_angle=(angle(qfft_Delta(:,2)));
Delta_in_angle1(p,i,j)=Delta_in_angle(q3);
Delta_in_fft1(p,i,j)=2*qfft_Delta(q3,2);    % multiply by 2, two-sided fft

end
end
end



% Xc, COMMAND
for j = 1:a
for i = 1:(logdata.data(1+ramp_loc1-1,iter_loc))
for p = 1:numfreq
    
    
h_stim=psd(spectrum.periodogram('rectangular'),input(input(1:stim_time(p),p,i,j)~=0,p,i,j),'Fs',Fs,'NFFT',NFFT);

L = length(input(Delta_in_center(1:stim_time(p),p,i,j)~=0,p,i,j));
%{
hfft  = fft(Delta_in_center(Delta_in_center(1:round(cycles/freq_stim(p)*1e4),p,i,j)~=0,p,i,j),NFFT)/L;
qfft_Delta(:,2)=2*hfft(1:NFFT/2+1);
qfft_Delta(:,1) = Fs/2*linspace(0,1,NFFT/2+1);
%}
clear y y0 wn

y0 = input(input(1:stim_time(p),p,i,j)~=0,p,i,j);
%{
wn = hann(length(y0));
for i2 = 1:length(wn)
    y(i2) = wn(i2)*y0(i2);
end
if isempty(y0)
    break;
end
%}
y=y0;

hfft  = fft(y,NFFT)/L;
qfft_Xc(:,2)=hfft(1:NFFT/2+1);
qfft_Xc(:,1) = Fs/2*linspace(0,1,NFFT/2+1);


clear q_stim;
q_stim(:,1)=h_stim.Frequencies;
q_stim(:,2)=h_stim.Data;

if p>1
    q_stim(length(psd_q1_stim),:)=0;
end

psd_q1_stim(:,p,i,j)=q_stim(:,1);    %frequency PRE
psd_q2_stim(:,p,i,j)=q_stim(:,2);    %psd PRE

q = find((freq_stim(p)-1)<psd_q1_stim(:,p,i,j)&1.05*(freq_stim(p)+1)>psd_q1_stim(:,p,i,j));
%q = q-2:1:q+2;
q2 = q(find(psd_q2_stim(q,p,i,j)==max(psd_q2_stim(q,p,i,j))));
q3 = findnearest(qfft_Xc(:,1),freq_stim(p));

Xc_amp1(p,i,j)=sqrt(psd_q1_stim(q2,p,i,j)*psd_q2_stim(q2,p,i,j));       % Output peak amplitude @ stimulus frequency
Xc_pow1(p,i,j)=psd_q2_stim(q2,p,i,j);

Xc_angle=(angle(qfft_Xc(:,2)));
Xc_angle1(p,i,j)=Xc_angle(q3);
Xc_fft1(p,i,j)=2*qfft_Xc(q3,2);     % multiply by 2, two-sided fft

end
end
end

warning on

%stim_angle=angle(stim_fft1)*180/pi;
%Delta_angle=angle(Delta_in_fft1)*180/pi;
phase_diff1=(angle(Delta_in_fft1)-angle(stim_fft1))*180/pi;

for j = 1:a
for i = 1:(logdata.data(1+ramp_loc1-1,iter_loc))
for p = 1:numfreq
    if phase_diff1(p,i,j)>180
        phase_diff1(p,i,j)=phase_diff1(p,i,j)-360;
    else if phase_diff1(p,i,j)<-180
            phase_diff1(p,i,j)=phase_diff1(p,i,j)+360;
        end
    end
end
end
end


for j = 1:a
for i = 1:(logdata.data(1+ramp_loc1-1,iter_loc))
for p = 1:numfreq
    stim_angle2(p,i,j) = atan2(real(stim_fft1(p,i,j)),imag(stim_fft1(p,i,j)));
    Delta_in_angle2(p,i,j) = atan2(real(Delta_in_fft1(p,i,j)),imag(Delta_in_fft1(p,i,j)));
end
end
end




%% AMPLITUDE ARRAY

warning off
clear q_stim psd_q1_stim psd_q2_stim qfft_stim

time = length(stim_center);    
T=1/Fs;
NFFT = 2^3*2^nextpow2(length(stim_center));

f = (Fs/2)*linspace(0,1,NFFT/2+1);
effective_delta_f = Fs/NFFT;
actual_delta_f = Fs/numel(time);



% STIMULUS, X
for j = 1:a
for i = 1:(logdata.data(1+ramp_loc1-1,iter_loc))
for p = 1:numfreq
    
    
h_stim=psd(spectrum.periodogram('rectangular'),stim_center(stim_center(1:stim_time(p),p,i,j)~=0,p,i,j),'Fs',Fs,'NFFT',NFFT);

L = length(stim_center(stim_center(1:stim_time(p),p,i,j)~=0,p,i,j));

clear y y0 wn

%{
hfft  = fft(stim_center(stim_center(1:stim_time(p),p,i,j)~=0,p,i,j),NFFT)/L;
qfft_stim(:,2)=2*hfft(1:NFFT/2+1);
qfft_stim(:,1) = Fs/2*linspace(0,1,NFFT/2+1);
%}
y0 = stim_center(stim_center(1:stim_time(p),p,i,j)~=0,p,i,j);
wn = hann(length(y0));
for i2 = 1:length(wn)
    y(i2) = wn(i2)*y0(i2);
end
if isempty(y0)
    break;
end
hfft  = fft(y,NFFT)/L;
qfft_stim(:,2)=2*hfft(1:NFFT/2+1);
qfft_stim(:,1) = Fs/2*linspace(0,1,NFFT/2+1);

clear q_stim;
q_stim(:,1)=h_stim.Frequencies;
q_stim(:,2)=h_stim.Data;

if p>1
    q_stim(length(psd_q1_stim),:)=0;
end

psd_q1_stim(:,p,i,j)=q_stim(:,1);    %frequency PRE
psd_q2_stim(:,p,i,j)=q_stim(:,2);    %psd PRE

q = find((freq_stim(1)-1)<psd_q1_stim(:,p,i,j)&1.05*(freq_stim(1)+1)>psd_q1_stim(:,p,i,j));
q = q-2:1:q+2;
q2 = q(find(psd_q2_stim(q,p,i,j)==max(psd_q2_stim(q,p,i,j))));

q3 = findnearest(qfft_stim(:,1),freq_stim(1));

stim_amp1(p,i,j)=sqrt(psd_q1_stim(q2,p,i,j)*psd_q2_stim(q2,p,i,j));       % Output peak amplitude @ stimulus frequency
stim_pow1(p,i,j)=psd_q2_stim(q2,p,i,j);

stim_angle=unwrap(angle(qfft_stim(:,2)));
stim_angle1(p,i,j)=stim_angle(q3);
stim_fft1(p,i,j)=2*qfft_stim(q3,2);


end
end
end


% FIBER, DELTA
for j = 1:a
for i = 1:(logdata.data(1+ramp_loc1-1,iter_loc))
for p = 1:numfreq
    
    
h_stim=psd(spectrum.periodogram('hann'),Delta_in_center(Delta_in_center(1:stim_time(p),p,i,j)~=0,p,i,j),'Fs',Fs,'NFFT',NFFT);

L = length(Delta_in_center(Delta_in_center(1:stim_time(p),p,i,j)~=0,p,i,j));

clear y y0 wn qfft_Delta
%{
hfft  = fft(Delta_in_center(Delta_in_center(1:stim_time(p),p,i,j)~=0,p,i,j),NFFT)/L;
qfft_Delta(:,2)=2*hfft(1:NFFT/2+1);
qfft_Delta(:,1) = Fs/2*linspace(0,1,NFFT/2+1);
%}


y0 = Delta_in_center(Delta_in_center(1:stim_time(p),p,i,j)~=0,p,i,j);
wn = hann(length(y0));
for i2 = 1:length(wn)
    y(i2) = wn(i2)*y0(i2);
end
if isempty(y0)
    break;
end
hfft  = fft(y,NFFT)/L;
qfft_Delta(:,2)=2*hfft(1:NFFT/2+1);
qfft_Delta(:,1) = Fs/2*linspace(0,1,NFFT/2+1);


clear q_stim;
q_stim(:,1)=h_stim.Frequencies;
q_stim(:,2)=h_stim.Data;

if p>1
    q_stim(length(psd_q1_stim),:)=0;
end

psd_q1_stim(:,p,i,j)=q_stim(:,1);    %frequency PRE
psd_q2_stim(:,p,i,j)=q_stim(:,2);    %psd PRE

q = find((freq_stim(1)-1)<psd_q1_stim(:,p,i,j)&1.05*(freq_stim(1)+1)>psd_q1_stim(:,p,i,j));
%q = q-2:1:q+2;
q2 = q(find(psd_q2_stim(q,p,i,j)==max(psd_q2_stim(q,p,i,j))));
q3 = findnearest(qfft_Delta(:,1),freq_stim(1));

Delta_in_amp1(p,i,j)=sqrt(psd_q1_stim(q2,p,i,j)*psd_q2_stim(q2,p,i,j));       % Output peak amplitude @ stimulus frequency
Delta_in_pow1(p,i,j)=psd_q2_stim(q2,p,i,j);

Delta_in_angle=unwrap(angle(0.5*qfft_Delta(:,2)));
Delta_in_angle1(p,i,j)=Delta_in_angle(q3);
Delta_in_fft1(p,i,j)=2*qfft_Delta(q3,2);

end
end
end

warning on

%stim_angle=angle(stim_fft1)*180/pi;
%Delta_angle=angle(Delta_in_fft1)*180/pi;
phase_diff1=(Delta_in_angle1-stim_angle1)*180/pi;

while max(max(max(phase_diff1)))>180 & min(min(min(phase_diff1)))<-180
for j = 1:a
for i = 1:(logdata.data(1+ramp_loc1-1,iter_loc))
for p = 1:numfreq
    if phase_diff1(p,i,j)>180
        phase_diff1(p,i,j)=phase_diff1(p,i,j)-360;
    else if phase_diff1(p,i,j)<-180
            phase_diff1(p,i,j)=phase_diff1(p,i,j)+360;
        end
    end
end
end
end
end

%%


% STIMULUS, X
for j = 1:a
for i = 1:(logdata.data(1+ramp_loc1-1,iter_loc))
for p = 1:numfreq
    
    [S,F,T,P] = spectrogram(stim_center(stim_center(1:stim_time(p),p,i,j)~=0,p,i,j),length(stim_center(stim_center(1:stim_time(p),p,i,j)~=0,p,i,j))/1,0,NFFT,Fs);
    psd_q1_stim(:,p,i,j)=F;
    %psd_q2_stim(:,p,i,j)=P;
    q = find((freq_stim(p)-1)<F&1.05*(freq_stim(p)+1)>F);
    
    stim_amp_STFT(p,i,j)=abs(max(P(q)));
    

    
end
end
end

%% Calculate force

k = [182 182 182 182 182 182 182 182];
k = k*1e-6;


for i = 1:length(amp_stim)
    for j = 1:length(k)
        Fcs(i,j) = amp_stim(i)*1e-9*k(j);
    end
end

for j = 1:a
for i = 1:(logdata.data(1+ramp_loc1-1,iter_loc))
for p = 1:numfreq
    sensitivity_ALL(p,i,j) = abs(stim_fft1(p,i,j))/Fcs(p,j)*1e-9;
end
end
end

% INPUT
ksf = k;       % Fiber stiffness (N/m)
gsf = 200e-9;       % Fiber damping (Ns/m)
amp = amp_stim;         % stimulus amplitude (m)

% Force = (ksf + i*w*gsf)*amp

for j = 1:a
for i = 1:(logdata.data(1+ramp_loc1-1,iter_loc))
for p = 1:numfreq
%sensitivity(p,i,j)=stim_amp(p,i,j)/(sqrt(ksf(j)^2+(gsf*freq_stim*2*pi)^2)*amp_stim(p)*1e-9)*1e-9;
force_in(p,i,j) = (ksf(j) + sqrt(-1)*freq_stim(1)*gsf)*Delta_in_fft(p,i,j);
res_func(p,i,j) = stim_fft1(p,i,j)/force_in(p,i,j);
sensitivity_resfunc(p,i,j) = norm(res_func(p,i,j));
phase_resfunc(p,i,j) = -angle(res_func(p,i,j));
end
end
end

%% Noise floor

N = 20;
freq_stim=8;
clear y
for i = 1:N
    y(:,i) = Smoothedphotodiodeoutputnm1(1+(i-1)*length(Smoothedphotodiodeoutputnm1)/N:i*length(Smoothedphotodiodeoutputnm1)/N);
end



time = length(y2);    
T=1/Fs;
NFFT = 2^3*2^nextpow2(length(y2));

f = (Fs/2)*linspace(0,1,NFFT/2+1);
effective_delta_f = Fs/NFFT;

for i = 1:N
    
y2=y(:,i);
L = length(y2);

clear y0 wn

%{
hfft  = fft(stim_center(stim_center(1:stim_time(p),p,i,j)~=0,p,i,j),NFFT)/L;
qfft_stim(:,2)=2*hfft(1:NFFT/2+1);
qfft_stim(:,1) = Fs/2*linspace(0,1,NFFT/2+1);
%}
y0 = y2;
wn = hann(length(y0));
for i2 = 1:length(wn)
    y(i2) = wn(i2)*y0(i2);
end
if isempty(y0)
    break;
end
hfft  = fft(y,NFFT)/L;
qfft_stim(:,2)=2*hfft(1:NFFT/2+1);
qfft_stim(:,1) = Fs/2*linspace(0,1,NFFT/2+1);

q3 = findnearest(qfft_stim(:,1),freq_stim);

stim_fft1(i)=2*qfft_stim(q3,2);
end

%% SPLIT SIGNAL


warning off
clear q_stim psd_q1_stim psd_q2_stim qfft_stim

time = length(stim_center);    
T=1/Fs;
NFFT = 2^3*2^nextpow2(length(stim_center));

f = (Fs/2)*linspace(0,1,NFFT/2+1);
effective_delta_f = Fs/NFFT;
actual_delta_f = Fs/numel(time);



% STIMULUS, X
for j = 1:a
for i = 1:(logdata.data(1+ramp_loc1-1,iter_loc))
for p = 1:numfreq
clear y3 y4    
y3 = stim_center(stim_center(1:stim_time(p),p,i,j)~=0,p,i,j);
y4 = y3(length(y3)*0.2:length(y3));
h_stim=psd(spectrum.periodogram('rectangular'),y4,'Fs',Fs,'NFFT',NFFT);

L = length(y4);

clear y y0 wn

%{
hfft  = fft(stim_center(stim_center(1:stim_time(p),p,i,j)~=0,p,i,j),NFFT)/L;
qfft_stim(:,2)=2*hfft(1:NFFT/2+1);
qfft_stim(:,1) = Fs/2*linspace(0,1,NFFT/2+1);
%}
y0 = y4;
wn = hann(length(y0));
for i2 = 1:length(wn)
    y(i2) = wn(i2)*y0(i2);
end
if isempty(y0)
    break;
end
hfft  = fft(y,NFFT)/L;
qfft_stim(:,2)=2*hfft(1:NFFT/2+1);
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

stim_angle=angle(0.5*qfft_stim(:,2));
stim_angle1(p,i,j)=stim_angle(q3);
stim_fft1(p,i,j)=2*qfft_stim(q3,2);


end
end
end


% FIBER, DELTA
for j = 1:a
for i = 1:(logdata.data(1+ramp_loc1-1,iter_loc))
for p = 1:numfreq
    
clear y3 y4    
y3 = Delta_in_center(Delta_in_center(1:stim_time(p),p,i,j)~=0,p,i,j);
y4 = y3(1:length(y3)*0.8);
h_stim=psd(spectrum.periodogram('rectangular'),y4,'Fs',Fs,'NFFT',NFFT);

L = length(y4);

clear y y0 wn

%{
hfft  = fft(stim_center(stim_center(1:stim_time(p),p,i,j)~=0,p,i,j),NFFT)/L;
qfft_stim(:,2)=2*hfft(1:NFFT/2+1);
qfft_stim(:,1) = Fs/2*linspace(0,1,NFFT/2+1);
%}
y0 = y4;
wn = hann(length(y0));
for i2 = 1:length(wn)
    y(i2) = wn(i2)*y0(i2);
end
if isempty(y0)
    break;
end
hfft  = fft(y,NFFT)/L;

qfft_Delta(:,2)=2*hfft(1:NFFT/2+1);
qfft_Delta(:,1) = Fs/2*linspace(0,1,NFFT/2+1);


clear q_stim;
q_stim(:,1)=h_stim.Frequencies;
q_stim(:,2)=h_stim.Data;

if p>1
    q_stim(length(psd_q1_stim),:)=0;
end

psd_q1_stim(:,p,i,j)=q_stim(:,1);    %frequency PRE
psd_q2_stim(:,p,i,j)=q_stim(:,2);    %psd PRE

q = find((freq_stim(p)-1)<psd_q1_stim(:,p,i,j)&1.05*(freq_stim(p)+1)>psd_q1_stim(:,p,i,j));
%q = q-2:1:q+2;
q2 = q(find(psd_q2_stim(q,p,i,j)==max(psd_q2_stim(q,p,i,j))));
q3 = findnearest(qfft_Delta(:,1),freq_stim(p));

Delta_in_amp1(p,i,j)=sqrt(psd_q1_stim(q2,p,i,j)*psd_q2_stim(q2,p,i,j));       % Output peak amplitude @ stimulus frequency
Delta_in_pow1(p,i,j)=psd_q2_stim(q2,p,i,j);

Delta_in_angle=angle(0.5*qfft_Delta(:,2));
Delta_in_angle1(p,i,j)=Delta_in_angle(q3);
Delta_in_fft1(p,i,j)=2*qfft_Delta(q3,2);

end
end
end

warning on

%stim_angle=angle(stim_fft1)*180/pi;
%Delta_angle=angle(Delta_in_fft1)*180/pi;
phase_diff1=(Delta_in_angle1-stim_angle1)*180/pi;

while max(max(max(phase_diff1)))>180 & min(min(min(phase_diff1)))<-180
for j = 1:a
for i = 1:(logdata.data(1+ramp_loc1-1,iter_loc))
for p = 1:numfreq
    if phase_diff1(p,i,j)>180
        phase_diff1(p,i,j)=phase_diff1(p,i,j)-360;
    else if phase_diff1(p,i,j)<-180
            phase_diff1(p,i,j)=phase_diff1(p,i,j)+360;
        end
    end
end
end
end
end

%% Phase locking

warning off
  
  freq =freq_stim;
  i = 1;
  t = linspace(0,3,length(Delta_in_center(Delta_in_center(1:stim_time(p),1,1,1)~=0,p,i,j)));
  for i = 1:length(freq)
      yw(:,i) = sin(freq(i)*2*pi*t);
  end
  

  
for j = 1:a
for i = 1:(logdata.data(1+ramp_loc1-1,iter_loc))
for p = 1:numfreq
    
clear y3 y4    
y3 = yw(:,p);
y4 = y3(1:length(y3)*0.8);
h_stim=psd(spectrum.periodogram('rectangular'),y4,'Fs',Fs,'NFFT',NFFT);

L = length(y4);

clear y y0 wn
  

y0 = y4;
wn = hann(length(y0));
for i2 = 1:length(wn)
    y(i2) = wn(i2)*y0(i2);
end
if isempty(y0)
    break;
end
hfft  = fft(y,NFFT)/L;
qfft_Delta(:,2)=hfft(1:NFFT/2+1);
qfft_Delta(:,1) = Fs/2*linspace(0,1,NFFT/2+1);


clear q_stim;
q_stim(:,1)=h_stim.Frequencies;
q_stim(:,2)=h_stim.Data;

if p>1
    q_stim(length(psd_q1_stim),:)=0;
end

psd_q1_stim(:,p,i,j)=q_stim(:,1);    %frequency PRE
psd_q2_stim(:,p,i,j)=q_stim(:,2);    %psd PRE

q = find((freq_stim(p)-1)<psd_q1_stim(:,p,i,j)&1.05*(freq_stim(p)+1)>psd_q1_stim(:,p,i,j));
%q = q-2:1:q+2;
q2 = q(find(psd_q2_stim(q,p,i,j)==max(psd_q2_stim(q,p,i,j))));
q3 = findnearest(qfft_Delta(:,1),freq_stim(p));

input_amp1(p,i,j)=sqrt(psd_q1_stim(q2,p,i,j)*psd_q2_stim(q2,p,i,j));       % Output peak amplitude @ stimulus frequency
input_pow1(p,i,j)=psd_q2_stim(q2,p,i,j);

input_angle=unwrap(angle(qfft_Delta(:,2)));
input_angle1(p,i,j)=input_angle(q3);
input_fft1(p,i,j)=2*qfft_Delta(q3,2);

end
end
end

input_angle1=input_angle1*180/pi;

for j = 1:a
for i = 1:(logdata.data(1+ramp_loc1-1,iter_loc))
for p = 1:numfreq
    if input_angle1(p,i,j)>180
        input_angle1(p,i,j)=input_angle1(p,i,j)-360;
    else if input_angle1(p,i,j)<-180
            input_angle1(p,i,j)=input_angle1(p,i,j)+360;
        end
    end
end
end
end

% PHASE DIFFERENCE HERE
phase_diff_fft1 = stim_angle1-input_angle1;

for j = 1:a
for i = 1:(logdata.data(1+ramp_loc1-1,iter_loc))
for p = 1:numfreq
    if phase_diff_fft1(p,i,j) > 180
        phase_diff_fft1(p,i,j) = phase_diff_fft1(p,i,j) - 360;
    else if phase_diff_fft1(p,i,j) < -180
            phase_diff_fft1(p,i,j) = phase_diff_fft1(p,i,j) + 360;
        end
    end
end
end
end


