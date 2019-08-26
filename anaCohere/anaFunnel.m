d = loadFiles('funnel'); 
%%
eps = d.scan.data.pulsegroups.varpar'; 
data = d.data{1}; 
b = scanRng(d.scan,1); 
dataNorm = [];
for i = 1:size(data,1)    
    p = polyfit(eps,data(i,:),1); 
    dataNorm(i,:) = data(i,:)-p(1)*eps-p(2);    
end
figure(123); clf; 
imagesc(eps,b*6,dataNorm); 
xlabel('\epsilon (mV)'); ylabel('Exchange (GHz)'); colorbar; 
a = gca; a.YDir = 'Normal'; 
%%
data = squeeze(d(1).data{1}); 
figure(400); clf; hold on;
mi = [];
for i = 1:5:size(data,1)-5
    dataSmooth = nanmean(data(i:i+5,40:100));
    
    plot(dataSmooth+0.04*i)
end
figure(500); clf; 
plot(mi)
%%
figure(145); clf; 
[~,mi] = max(dataNorm,[],2); 
epsSTP = eps(mi); 
inds = 1:40; 
e = epsSTP(inds); 
epsClean = removeOutliers(epsSTP(inds));
plot(b(inds)*6,e); 