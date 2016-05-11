function plothbvssens(filename,stiffforce)
% This function calculates the sensitivity and vector strength across
% various sets of control parameters.
%
% plothbvssens(filename,stiffforce)
%
% stiffforce: 1=stiffness;2=force
%
% Joshua D. Salvi
% jsalvi@rockefeller.edu

clc;
setfiguredefaults
close all

% Load the data file
disp('Loading..');
load(filename);
disp('Done.');

analyzeyn=1;

while analyzeyn == 1
    
close all;
% Cell to analyze?
ksf = input('Fiber stiffness (µN/m): ');
disp(['Raw indices: ' num2str(raw)]);
disp(['Averaged indices: ' num2str(nonraw)]);
disp('Select an index for analysis...');
ind = input('Index: ');
if isempty(intersect(ind,raw)) == 0
    disp('Selected RAW index.');
    mm = 2;
elseif isempty(intersect(ind,nonraw)) == 0
    disp('Selected AVERAGED index.');
    mm = 1;
else
    disp('Did not selected appropriate index.');
    return;
end

if mm ==2
    
% Time vector
clear dt tvec
dt = 1/Fs;
tvec = 0:dt:length(Xd_pulse{1,ind})*dt-dt;


% Plot the traces
sizeX = size(Xd_pulse{1,ind});  % set figure defaulta
setfiguredefaults(sizeX(2));    

for j = 1:size(Xd_pulse,1)      % number of averages
    pulseind0(j) = isempty(Xd_pulse{j,ind});
end
pulseind = find(pulseind0==0); clear pulseind0
pulseL = length(pulseind);
XsegL = floor(length(tvec));

hrt(1)=figure(1);       % plot X
for j = 1:pulseL
    subplot(1,pulseL,j);
    for i = 1:sizeX(2)
        plot(tvec(1:length(Xd_pulse{j,ind+1}(:,1))),Xd_pulse{j,ind}(:,i)-smooth(Xd_pulse{j,ind}(:,i),XsegL/10));hold all;
    end
    grid on;axis tight;
    set(1,'WindowStyle','docked')
    grid on;
end

hrt(2)=figure(2);       % plot Fe
for j = 1:pulseL
    subplot(1,pulseL,j);
    sizeF = size(Fe_pulse{j,ind+1});
    for i = 1:sizeF(2)
        plot(tvec(1:length(Fe_pulse{j,ind+1}(:,1))),Fe_pulse{j,ind+1}(:,i)-mean(Fe_pulse{j,ind+1}(:,i)));hold all;
    end
    grid on;axis tight;
    set(2,'WindowStyle','docked')
    grid on;
end

% FFT Preliminaries
NFFT = (2^4)*2^nextpow2(numel(tvec));
nw = 1;     % one window

welchwin = round(XsegL);
NPSD = floor(NFFT/nw);
noverlap = 0;       % zero overlap
winfunc = rectwin(welchwin);    % window/taper
f = Fs/2*linspace(0,1,NFFT/2+1);
% Test function
freq = 1;
Xsine = sin(2*pi*freq.*tvec);
Xsinefft = fft(Xsine,NFFT)./XsegL; Xsinefft = Xsinefft(1:NFFT/2+1);
winpeaknorm = sqrt(max(abs(Xsinefft)).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL;

% Calculate FFT for Fe
j=1;
sizeF = size(Fe_pulse{j,ind+1});
for j = 1:pulseL
    for i = 1:sizeF(2)
        Fe_fft{j,i} = fft(Fe_pulse{j,ind+1}(:,1)-mean(Fe_pulse{j,ind+1}(:,1)),NFFT)./XsegL; Fe_fft{j,i} = Fe_fft{j,i}(1:NFFT/2+1);
        Fepk00 = find(abs(Fe_fft{j,i}(:)) == max(abs(Fe_fft{j,i}(:))));
        Fepk0(j,i) = Fepk00(1);
        Fepk(j,i) = Fe_fft{j,i}(Fepk0(j,i));
        freqstim(j,i) = f(Fepk0(j,i));
    end
end
% Calculate FFT for Xd
sizeX = size(Xd_pulse{1,ind});
for j = 1:pulseL
    for i = 1:sizeX(2)
        Xd_fft{j,i} = fft(Xd_pulse{j,ind}(:,i)-smooth(Xd_pulse{j,ind}(:,i),XsegL/10),NFFT)./XsegL; Xd_fft{j,i} = Xd_fft{j,i}(1:NFFT/2+1);
        Xdpk0(j,i) = find(abs(Xd_fft{j,i}) == max(abs(Xd_fft{j,i})));
        Xdpk(j,i) = max(Xd_fft{j,i}(Fepk0(j,1)-10:Fepk0(j,1)+10));
    end
end
% Calculate FFT for Xo
sizeX = size(Xo_pulse{1,ind});
for j = 1:pulseL
    for i = 1:sizeX(2)
        Xo_fft{j,i} = fft(Xo_pulse{j,ind}(:,i)-smooth(Xo_pulse{j,ind}(:,i),XsegL/10),NFFT)./XsegL; Xo_fft{j,i} = Xo_fft{j,i}(1:NFFT/2+1);
        Xopk0(j,i) = find(abs(Xo_fft{j,i}) == max(abs(Xo_fft{j,i})));
        Xopk(j,i) = max(Xo_fft{j,i}(Fepk0(j,1)-10:Fepk0(j,1)+10));
    end
end

% Plot the individual results
if stiffforce == 1
    mu = logdata.data(2,138:138+sizeX(2)-1);
else
    mu = logdata.data(2,138+sizeX(2):138+2*sizeX(2)-1);
end
setfiguredefaults(sizeX(2));
hrt(3) = figure(3);
for j = 1:pulseL
    subplot(1,pulseL,j)
    plot(mu,abs(sqrt(Xdpk(j,:)).*1e-9)./(abs(sqrt(Xopk(j,:))).*ksf.*1e-6.*1e-9) .* 1e-3);title(num2str(j));
    xlabel('Control parameter');ylabel('Xd');
    set(3,'WindowStyle','docked')
end



% Calculate vector strength
for j = 1:pulseL
    for i = 1:sizeX(2)
        [VS(j,i), rayleigh_p(j,i), rayleigh_stat(j,i)] = vscalc2(Xd_pulse{j,ind}(:,i)-smooth(Xd_pulse{j,ind}(:,i),XsegL/10),Xo_pulse{j,ind}(:,i),1,1e-10);
        [VS2(j,i), rayleigh_p2(j,i), rayleigh_stat2(j,i)] = vscalc2(Xd_pulse{j,ind}(:,i)-smooth(Xd_pulse{j,ind}(:,i),XsegL/10),Fe_pulse{j,ind+1}(:,1),1,1e-10);
    end
end

hrt(4) = figure(4);
for j = 1:pulseL
    subplot(1,pulseL,j)
    plot(mu,VS(j,:));title(num2str(j));
    set(4,'WindowStyle','docked')
end
%[pxxf{j,k,l},fxx{j,k,l}]=pwelch(Xo_pulse{k,raw(j)}(:,l),[],[],[],Fs);

analyzedind = input('Analyze which indices?:  ');

for j = 1:size(Xdpk,2)
    Xdpkmean(j) = mean(abs(Xdpk(analyzedind,j)));
    Fepkmean(j) = mean(abs(Fepk(analyzedind)));
    Xopkmean(j) = mean(Xopk(analyzedind,j));
    sensmean(j) = abs(sqrt(Xdpkmean(j))*1e-9) / abs(sqrt(mean(Xopk(analyzedind,j))).*ksf.*1e-6.*1e-9) * 1e-3;
    sensmean2(j) = abs(sqrt(Xdpkmean(j))*1e-9) / abs(sqrt(Fepkmean(j))*1e-12) *1e-3;
    Xdpksem(j) = std(abs(sqrt(Xdpk(analyzedind,j))))/sqrt(length(analyzedind));
    Xopksem(j) = std(abs(sqrt(Xopk(analyzedind,j))))/sqrt(length(analyzedind));
    Fepksem(j) = std(abs(sqrt(Fepk(analyzedind))))/sqrt(length(analyzedind));
    senssem(j) = Xdpksem(j)*1e-9 / abs(sqrt(mean(Xopk(analyzedind,j))).*ksf.*1e-6.*1e-9) * 1e-3 ./sqrt(length(analyzedind));
    senssem2(j) = abs(sqrt(Xdpksem(j))*1e-9) /  abs(sqrt(Fepkmean(j))*1e-12) * 1e-3 ./sqrt(length(analyzedind));
    VSmean(j) = mean(VS(analyzedind,j));
    VSsem(j) = std(VS(analyzedind,j))/sqrt(length(analyzedind));
    VSmean2(j) = mean(VS2(analyzedind,j));
    VSsem2(j) = std(VS2(analyzedind,j))/sqrt(length(analyzedind));
end

hrt(5) = figure(5);
subplot(2,1,1);
errorbar(mu,sensmean,senssem); xlabel('Control parameter'); ylabel('Sensitivity (km/N)');
subplot(2,1,2);
errorbar(mu,VSmean,VSsem); xlabel('Control parameter'); ylabel('Vector strength');
set(5,'WindowStyle','docked')

hrt(6) = figure(6);
subplot(2,1,1);
errorbar(mu,sensmean2,senssem2); xlabel('Control parameter'); ylabel('Sensitivity (km/N)');
subplot(2,1,2);
errorbar(mu,VSmean2,VSsem2); xlabel('Control parameter'); ylabel('Vector strength');
set(6,'WindowStyle','docked')

end

saveyn = input('Save the figures? (1=yes): ');
hrt(hrt==0)=[];
if saveyn==1
    disp('Saving...');
    filename2 = filename(1:end-18);
    savefig(hrt,sprintf('%s%s%s%s',filename2,'FX-Figures-ind',num2str(ind),'.fig'));
    disp('Finished.');
else
    disp('Not saved.');
end
saveyn = input('Save the data? (1=yes): ');
if saveyn==1
    disp('Saving...');
    filename2 = filename(1:end-18);
    save([filename2 'FXdata-ind' num2str(ind) '.mat'],'VS','VS2','VSmean','VSmean2','VSsem','VSsem2','senssem','sensmean','Xdpkmean','Xopkmean','Fepkmean','Xdpksem','Xopksem','Fepksem','Xd_fft','Xo_fft','Fe_fft','freqstim','ind','mu','rayleigh_p','rayleigh_stat','Xopk','Xdpk','Fepk','ksf','analyzedind');
    disp('Finished.');
else
    disp('Not saved.');
end

analyzeyn = input('Analyze again? (1=yes):  ');
end

end
