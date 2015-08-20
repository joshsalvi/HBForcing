analyzedind=[2:4]
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
inds=[1     2     5    7   8  9 10  11 12      14        16    18     20   ];
subplot(2,1,1);errorbar(mu(inds).*2,sensmean(inds),senssem(inds)./sqrt(4),'r');axis tight
set(gca,'XDir','reverse')
xlabel('Control parameter');ylabel('Sensitivity (km/N)');
subplot(2,1,2);errorbar(mu(inds).*2,VSmean(inds),VSsem(inds)./sqrt(4),'r');axis tight
set(gca,'XDir','reverse')
xlabel('Control parameter');ylabel('Vector strength');
