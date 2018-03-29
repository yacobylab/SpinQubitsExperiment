function colorPulsePlot(num)

filename = sprintf('Z:/qDots/data/data_2015_11_05/sm_RamseyL_%02d',num);
load(filename);

ramseyData = data{1}; 
meanData = squeeze(nanmean(ramseyData)); 

figure(15); clf; 
imagesc(meanData)

%%
figure(14); clf; 
plot(meanData(9,:))

%% 
for i = 1:length(scan.data.pulsegroups) 
    eps(i) = scan.data.pulsegroups(i).params(2); 
end

%%
figure(14); clf; 
imagesc(1:200,eps(1:10),meanData(1:10,:))

end
