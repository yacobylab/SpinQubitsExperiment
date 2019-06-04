d = procPlsData; 
%%
d = procPlsData('','none') 
%%
data = squeeze(nanmean(d.data{1})); 

figure(400); clf; 
imagesc(d.xv{1},d.tv,data,'ButtonDownFcn',@btn);
xlabel('Time (ns)'); ylabel('\epsilon'); 
a = gca; a.YDir = 'normal'; 
colorbar; 
%% 
plotChrg('sens stop noppt',f); 
%%
figure(3); hold on; 
for i = 1:length(d)
exch=8*d(i).scan.data.pulsegroups(4).dict{2}{1}.exch.val*1e-3;
plot(exch(1),exch(2),'.','MarkerSize',20)
end
%%
fidelity = [d.fidelity];
meanvals = [d.meanvals];
meanvals = reshape(meanvals,2,length(meanvals)/2); 
plot(fidelity)