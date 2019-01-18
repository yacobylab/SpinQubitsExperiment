function colorPulsePlot(file,num) 
% Plot 2D pulsed data in color plot with no processing.
% function colorPulsePlot(file,num) 
% Plot line num as a line plot. 9 by default. 

if ~exist('file','var') || isempty(file)
    d = loadFiles; 
else 
    d = load(file); 
end
if ~exist('num','var') || isempty(num), num = 9; end
ramseyData = d.data{1};
meanData = squeeze(nanmean(ramseyData));

figure(15); clf;
imagesc(meanData)

figure(14); clf;
plot(meanData(9,:))

%for i = 1:length(scan.data.pulsegroups)
%    eps(i) = scan.data.pulsegroups(i).params(2);
%end
end