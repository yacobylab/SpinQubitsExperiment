d = procPlsData; 
%%
d = procPlsData('','2d noscale') 
%% Plot the proc pls data from a single scan, allowing one to click on it. 
i = 1;
data = squeeze(nanmean(d(i).data{1})); 

figure(400); clf; 
imagesc(d(i).xv{1},d(i).tv,data,'ButtonDownFcn',@btn);
xlabel('Time (ns)'); ylabel('\epsilon'); 
a = gca; a.YDir = 'normal'; 
colorbar; 
%%
[f,fp] = getFiles; 
f = fullfile(fp,f); 
%% Load the most recent charge scan
plotChrg('sens stop noppt',f); 
%% This cell can be used to plot the direction of exchange for the different pulses. 
figure(3); hold on;
for i = 1:length(d)
    exch=6*d(i).scan.data.pulsegroups(4).dict{2}{1}.exch.val*1e-3;
    plot(exch(1),exch(2),'.','MarkerSize',20)
end
%% Plot the fidelity across a bunch of scans. 
fidelity = [d.fidelity];
meanvals = [d.meanvals];
meanvals = reshape(meanvals,2,length(meanvals)/2); 
figure(1020); 
plot(fidelity)

%% 
d = loadFiles('*load*'); 
figure(44); 
imagesc(squeeze(nanmean(d(1).data{1})),'ButtonDownFcn',@btn)
a = gca; a.YDir = 'normal'; 
%% Just zoom
d = loadFiles('*pulsedZoom*');
%d = loadFiles('*zoom*');
%%
figure(444); clf 

subplot(2,2,1); 
data = squeeze(d(1).data{1}); 
coeff=fitPlane(data);
[mx,my]=meshgrid(1:size(data,2),1:size(data,1));
d1=data-mx*coeff(1)-my*coeff(2)-coeff(3);

imagesc(d1); 
a = gca; a.YDir = 'Normal'; 

subplot(2,2,2); 
data = squeeze(d(2).data{1}); 
coeff=fitPlane(data);
[mx,my]=meshgrid(1:size(data,2),1:size(data,1));
d2=data-mx*coeff(1)-my*coeff(2)-coeff(3);

imagesc(d2); 
a = gca; a.YDir = 'Normal'; 
subplot(2,2,3); 
imagesc(d2-d1); 
a = gca; a.YDir = 'Normal'; 
%%
for i = 1:2
%data = squeeze(nanmean(d(i).data{1}));
data = squeeze(d(i).data{1});
figure(100); subplot(1,2,i)
imagesc(-data); 
a = gca; a.YDir = 'normal';
end
%%
data = squeeze(nanmean(d(1).data{1}));
figure(100); clf;
imagesc(data); 
a = gca; a.YDir = 'normal';

%%
i = 1;
data = squeeze(d(i).data{1}); 

figure(400); clf; 
imagesc(d(i).xv{1},d(i).tv,data,'ButtonDownFcn',@btn);
xlabel('Time (ns)'); ylabel('\epsilon'); 
a = gca; a.YDir = 'normal'; 
colorbar; 
%%

%%
data = squeeze(d(1).data{1}); 
figure(400); clf; hold on;
mi = [];
for i = 1:10:size(data,1)-10
    dataSmooth = nanmean(data(i:i+9,1:40));
    [~,mi(end+1)] = max(dataSmooth); 
    plot(dataSmooth+0*i)
end
%% 
d= loadFiles; 
%%
close all; 
imagesc(data)