warning off
clear q_stim psd_q1_stim psd_q2_stim qfft_stim

time = length(stim_center);    
T=1/Fs;
NFFT = 2^3*2^nextpow2(length(stim_center));

f = (Fs/2)*linspace(0,1,NFFT/2+1);

% Phase difference
for j = 1:a
for i = 1:(logdata.data(1+ramp_loc1-1,iter_loc))
for p = 1:numfreq
    
                                                                    % X,
                                                                    % HAIR
                                                                    % BUNDLE
    L = length(stim_center(stim_center(1:round(cycles/freq_stim(p)*1e4),p,i,j)~=0,p,i,j));
    
    clear y0;
    y0 = stim_center(stim_center(1:round(cycles/freq_stim(p)*1e4),p,i,j)~=0,p,i,j);
    
    %%%
    wn = hann(length(y0));
    for i2 = 1:length(wn)
        y(i2) = wn(i2)*y0(i2);
    end
    if isempty(y0)
        break;
    end
    %%%
    
    %y=y0;
    
    hfft  = fft(y,NFFT)/L;
    qfft_stim1(:,2,p,i,j)=2*hfft(1:NFFT/2+1);
    qfft_stim1(:,1,p,i,j) = Fs/2*linspace(0,1,NFFT/2+1);
    
    q3 = findnearest(qfft_stim1(:,1,p,i,j),freq_stim(p));
    
    stim_fft_phase(p,i,j)=qfft_stim1(q3,2,p,i,j);
    
    
    clear y1 y2 hfft;
    y1 = input(1:length(y),p,i,j);                                  % Xc, COMMAND SIGNAL
    
    %%%
    wn = hann(length(y0));
    for i2 = 1:length(wn)
        y2(i2) = wn(i2)*y1(i2);
    end
    if isempty(y0)
        break;
    end
    %%%
    
    %y2=y1;
    
    hfft = fft(y2,NFFT)/L;
    qfft_stim2(:,2,p,i,j)=2*hfft(1:NFFT/2+1);
    qfft_stim2(:,1,p,i,j) = Fs/2*linspace(0,1,NFFT/2+1);
    
    Xc_fft_phase(p,i,j)=qfft_stim2(q3,2,p,i,j);
   
    
    clear yc yc2 hfft t;
    t=linspace(0,cycles/freq_stim(p),length(y0));
    yc = -sin(2*pi*freq_stim(p)*t);                                 % CONSTRUCTED SIGNAL
    
    %%%
    wn = hann(length(y0));
    for i2 = 1:length(wn)
        yc2(i2) = wn(i2)*yc(i2);
    end
    if isempty(y0)
        break;
    end
    %%%
    %yc2=yc;

    hfft = fft(yc2,NFFT)/L;
    qfft_stim3(:,2,p,i,j)=2*hfft(1:NFFT/2+1);
    qfft_stim3(:,1,p,i,j) = Fs/2*linspace(0,1,NFFT/2+1);
    
    yc_fft_phase(p,i,j)=qfft_stim3(q3,2,p,i,j);    
end
end
end


stim_Xc_phase2 = (atand(imag(Xc_fft_phase)./real(Xc_fft_phase))-atand(imag(stim_fft_phase)./real(stim_fft_phase)));


for i = [1 3 4 5 7]
figure
o=i;
subplot(4,1,1);plot(freq_stim(1:9),imag(stim_fft_phase(:,1,o))./real(stim_fft_phase(:,1,o)),'k');hold on;plot(freq_stim(1:9),imag(Xc_fft_phase(:,1,o))./real(Xc_fft_phase(:,1,o)),'r');
subplot(4,1,2);plot(freq_stim(1:9),atan(imag(stim_fft_phase(:,1,o))./real(stim_fft_phase(:,1,o))),'k');hold on;plot(freq_stim(1:9),atan(imag(Xc_fft_phase(:,1,o))./real(Xc_fft_phase(:,1,o))),'r');
subplot(4,1,3);plot(freq_stim(1:9),atan2(imag(stim_fft_phase(:,1,o)),real(stim_fft_phase(:,1,o))),'k');hold on;plot(freq_stim(1:9),atan2(imag(Xc_fft_phase(:,1,o)),real(Xc_fft_phase(:,1,o))),'r');
subplot(4,1,4);plot(freq_stim(1:9),angle(stim_fft_phase(:,1,o)),'k');hold on;plot(freq_stim(1:9),angle(Xc_fft_phase(:,1,o)),'r');
end

