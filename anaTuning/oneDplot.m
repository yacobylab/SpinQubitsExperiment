filtList = {'stp','tl','loadPos','load_','lead'}; 
%filt = 'stp'; 
%filt = 'lead'; 
%filt = 'loadPos'; 
%filt = 'tl'; 
%filt = 'load_'; 
%filt = 'loadTime'; 
d = loadFiles;
%%
nPlots = ceil(length(d)/12); 
axesList = [];
plotSpace = {4,3, [0.063, 0.073], [0.06 0.045], [0.06, 0.1]};            
for i = 1:nPlots 
    figure(12+i); clf; 
    axesList = [axesList, tight_subplot(plotSpace{:})]; 
end
 for i = 1:length(d) 
     plot(axesList(i),nanmean(d(i).data{1})); 
 end 
 %%
 nPlots = ceil(length(d)/12); 
axesList = [];
plotSpace = {4,3, [0.063, 0.073], [0.06 0.045], [0.06, 0.1]};            
for i = 1:nPlots 
    figure(12+i); clf; 
    axesList = [axesList, tight_subplot(plotSpace{:})]; 
end
 for i = 1:length(d) 
     plot(axesList(i),squeeze(nanmean(d(i).data{1}))'); 
 end 