% input noise levels
noiselevel = [0.01 0.5 0.75 1 1.5 2 0.9 3 5 6 7 8];
sizeX = size(Xd_pulse);
sizeX0 = size(Xd_pulse{1,raw(1)},2);
ksf = [300 300 300 300 300];
ksf = ksf.*1e-6;

for j = 1:length(raw)
for k = 1:sizeX(1)
for l = 1:sizeX0
    x = Xd_pulse{k,raw(j)}(:,l);
    x = x - mean(x);
    L = length(x);
    NFFT = 2^6*2^nextpow2(length(x));
    fxx{j,k,l} = Fs/2*linspace(0,1,NFFT/2+1);
    pxx{j,k,l} = fft(x,NFFT)./L; pxx{j,k,l} = pxx{j,k,l}(1:NFFT/2+1);
end
end
end
disp('1')
m=1;
for j = 1:length(raw)
for k = 1:sizeX(1)
for l = 1:sizeX0
q=findnearest(fxx{j,k,l},10);
pow1(j,k,l) = max(pxx{j,k,l}(q-3:q+3));
ampl(j,k,l) = sqrt(pow1(j,k,l));
end
end
end
disp('2')
[a b]=sort(noiselevel);
for k = 1:sizeX0;
    for j = 1:sizeX(1);
        m=1;
        for l = b;
ampl2(m,j,k) = ampl(l,j,k);m = m+1;
end
end
end



for j = 1:length(raw)
for k = 1:sizeX(1)
for l = 1:sizeX0
     x = Xo_pulse{k,raw(j)}(:,l);
     x = x - mean(x);
    L = length(x);
    NFFT = 2^6*2^nextpow2(length(x));
    fxx{j,k,l} = Fs/2*linspace(0,1,NFFT/2+1);
    pxxf{j,k,l} = fft(x,NFFT)./L; pxx{j,k,l} = (pxx{j,k,l}(1:NFFT/2+1));
end
end
end;
disp('3');

for j = 1:length(raw)
for k = 1:sizeX(1)
for l = 1:sizeX0
q=findnearest(fxx{j,k,l},10);
powf(j,k,l) = max(pxxf{j,k,l}(q));
amplf(j,k,l) = sqrt(powf(j,k,l));
end
end
end
for k = 1:sizeX0;for j = 1:sizeX(1);m=1;for l = b;
amplf2(m,j,k) = amplf(l,j,k);m=m+1;
end
end
end
disp('4')
noiselevel2=sort(noiselevel);
for j = 1:length(raw)
for k = 1:sizeX(1)
for l = 1:sizeX0
sens(j,k,l) = ampl2(j,k,l)*1e-9/(amplf2(j,k,l)*1e-9*ksf(k));
end
end
end
for j = 1:length(raw)
for k = 1:sizeX(1)
for l = 1:sizeX0
vs(j,k,l) = vscalc2(Xd_pulse{k,raw(j)}(:,l)-mean(Xd_pulse{k,raw(j)}(:,l)),Xo_pulse{k,raw(j)}(:,l)-mean(Xo_pulse{k,raw(j)}(:,l)),1,1e-10);
end
end
end
disp('5')
for k = 1:sizeX0;for j = 1:sizeX(1);m=1;for l = b;
vs2(m,j,k) = vs(l,j,k);m=m+1;
end
end
end


for j = 1:length(raw)
for k = 1:sizeX(1)
for l = 1:sizeX0
vssine(j,k,l) = vscalc2(Xd_pulse{k,raw(j)}(:,l)-mean(Xd_pulse{k,raw(j)}(:,l)),5*sin(2*10*pi*tvec{1,1}(1:length(Xd_pulse{k,raw(j)}(:,l)))),1,1e-10);
end
end
end
for k = 1:sizeX0;for j = 1:sizeX(1);m=1;for l = b;
vssine2(m,j,k) = vs2(l,j,k);m=m+1;
end
end
end
disp('done')
%%

close all;

sizeS = size(sens,1);
jj=[1:5];
inds = 1:sizeS;
inds = [4:sizeS-2 sizeS];

for k = 1:sizeX0
    
    subplot_tight(2,sizeX0,k);hold on;errorbar((noiselevel2(inds).^2)./2,mean(abs(sens(inds,jj,k)),2)./sqrt(length(jj))./1e3,std(abs(sens(inds,jj,k))./sqrt(length(jj))./1e3,[],2),'k')
    set(gca,'XScale','log');ylabel('sensitivity');
    subplot_tight(2,sizeX0,k+sizeX0);hold on;errorbar((noiselevel2(inds).^2)./2,mean(vssine2(inds,jj,k),2)./sqrt(length(jj)),std(vssine2(inds,jj,k)./sqrt(10),[],2),'k')
    set(gca,'XScale','log');ylabel('VS');xlabel('Control parameter')
end
%% plot all on individual plots
close all;

sizeS = size(sens,1);
inds = 1:sizeS;
inds = [1:sizeS];
for jj = 1:sizeX(1)
    figure(jj);
for k = 1:sizeX0
    subplot_tight(2,sizeX0,k);hold on;errorbar((noiselevel2(inds).^2)./2,abs(mean(sens(inds,jj,k),2))./sqrt(length(jj))./1e3,std(abs(sens(inds,jj,k))./sqrt(length(jj))./1e3,[],2),'k')
    set(gca,'XScale','log');ylabel('sensitivity');xlabel('Noise level (nm^2)')
    subplot_tight(2,sizeX0,k+sizeX0);hold on;errorbar((noiselevel2(inds).^2)./2,mean(vs2(inds,jj,k),2)./sqrt(length(jj)),std(vs2(inds,jj,k)./sqrt(length(jj)),[],2),'k')
    set(gca,'XScale','log');ylabel('VS');xlabel('Noise level (nm^2)')
end
end
%%
close all;

sizeS = size(sens,1);
jj=[1:5];
inds = 1:sizeS;
inds = [5:sizeS];

for k = 1:sizeX0
    subplot_tight(2,sizeX0,k);hold on;errorbar((noiselevel2(inds).^2)./2,abs(mean(sens(inds,jj,k),2))./sqrt(length(jj))./1e3,std(abs(sens(inds,jj,k))./sqrt(length(jj))./1e3,[],2),'k')
    set(gca,'XScale','log');ylabel('sensitivity');xlabel('Noise level (nm^2)')
    subplot_tight(2,sizeX0,k+sizeX0);hold on;errorbar((noiselevel2(inds).^2)./2,mean(vs2(inds,jj,k),2)./sqrt(length(jj)),std(vs2(inds,jj,k)./sqrt(length(jj)),[],2),'k')
    set(gca,'XScale','log');ylabel('VS');xlabel('Noise level (nm^2)')
end

%% choose the plots to overlay
close all

indsanalyze = [1:5 ];
jj= [3:5];
inds = 1:sizeS;
inds = [4:sizeS];
setfiguredefaults(length(indsanalyze));
figure;
for j = 1:length(indsanalyze)
    subplot(2,1,1);errorbar((noiselevel2(inds).^2)./2,abs(mean(sens(inds,jj,indsanalyze(j)),2))./sqrt(length(jj))./1e3,std(abs(sens(inds,jj,indsanalyze(j)))./sqrt(length(jj))./1e3,[],2));hold all;
    set(gca,'XScale','log');ylabel('sensitivity');xlabel('Noise level (nm^2)')
    subplot(2,1,2);errorbar((noiselevel2(inds).^2)./2,mean(vs2(inds,jj,indsanalyze(j)),2)./sqrt(length(jj)),std(vs2(inds,jj,indsanalyze(j))./sqrt(length(jj)),[],2));hold all;
    set(gca,'XScale','log');ylabel('VS');xlabel('Noise level (nm^2)')
end
subplot(2,1,1);title(['Avgs: ' num2str(jj) '; Indices: ' num2str(indsanalyze) '; ksf_m_i_n= ' num2str(ksf(1).*1e6) '; ksf_m_a_x= ' num2str(ksf(end).*1e6)]);

%%
save('/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/2015-08-11.01/Ear 1/Cell 6/SRdata.mat',...
    'ampl','ampl2','amplf','amplf2','noiselevel','noiselevel2','Fs','inds','indsanalyze',...
    'jj','ksf','sens','vs2','vs','vssine','vssine2','-v7.3');
