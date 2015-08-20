%% set initial conditions
close all

mu = linspace(0,50,20);

inds = [1:20];

ws=1;
ws2=5000;
indstart=1;
indend=15e3;
for j = 1:length(Xd_pulse)
x = Xd_pulse{j} - smooth(Xd_pulse{j},ws2);
x = smooth(x,ws);
x = x(indstart:indend);
[dip(j),dipp(j)]=HartigansDipSignifTest(x,1);
RMSmag(j) = std(x);

end

subplot(2,1,1);plot(mu(inds),dip(inds));axis tight
set(gca,'XDir','reverse')
xlabel('Control parameter');ylabel('Dip statistic');
subplot(2,1,2);plot(mu(inds),RMSmag(inds));axis tight
set(gca,'XDir','reverse')
xlabel('Control parameter');ylabel('RMS mag.');

%% bootstrapping

Nboot = 1e3;

for j = 1:length(Xd_pulse)
x = Xd_pulse{j} - smooth(Xd_pulse{j},ws2);
x = smooth(x,ws);
x = x(indstart:indend);
[dip(j),dipp(j)]=HartigansDipSignifTest(x,Nboot);
RMSmag(j) = std(x);
RMSmagboot{j} = bootstrp(Nboot,@(x) std(x), x);
dipboot{j} = bootstrp(Nboot,@(x) HartigansDipSignifTest(x,1), x);
RMSmagsem(j) = std(RMSmagboot{j});
dipsem(j) = std(dipboot{j});
disp(num2str(j));
end

subplot(2,1,1);errorbar(mu(inds),dip(inds),dipsem(inds));axis tight
set(gca,'XDir','reverse')
xlabel('Control parameter');ylabel('Dip statistic');
subplot(2,1,2);errorbar(mu(inds),RMSmag(inds),RMSmagsem(inds));axis tight
set(gca,'XDir','reverse')
xlabel('Control parameter');ylabel('RMS mag.');

%% save
save('/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/2015-08-10.01/Ear 1/Cell 1/dipRMSdata-Xd_pulse{1,3}.mat',...
    'Xd_pulse','indstart','indend','ws','ws2','Fs','mu','dip','dipp','RMSmag','RMSmagboot','dipboot',...
    'RMSmagsem','dipsem','Nboot','inds','-v7.3');
