function chrgTemplate(tsize)
% Generate templates for automated triple point detection.
% function chrgTemplate(tsize)
% tsize: def is 10
% sets tuneData.chrg.imgl, imgr to be a square of size tsize around triple
% point
% Loads file and asks user to click triple points. 

global tuneData;
if ~exist('tsize','var'), tsize=10; end
[file,fpath] = uigetfile(sprintf('%s/sm_chrg*', tuneData.dir));
d=load([fpath file]); data = d.data{1}; 
dataDiff = diff(data, [], 2);
dataDiff = dataDiff- median(dataDiff(:)); % Clean up data by subtracting off median. 
dataDiff = dataDiff .* sign(mean(dataDiff(:)));
m=nanmean(dataDiff(:)); s=nanstd(dataDiff(:));
dataDiff(abs(dataDiff)-m>8*s)=NaN; % Remove outliers

xvals = linspace(d.scan.loops(1).rng(1),d.scan.loops(1).rng(2),size(d.data{1},2)); xvals = xvals(1:end-1)/2+xvals(2:end); 
yvals = linspace(d.scan.loops(2).rng(1),d.scan.loops(2).rng(2),size(d.data{1},1));
dx = abs(xvals(2)-xvals(1)); dy = abs(yvals(2)-yvals(1)); 
figure(1); clf;
imagesc(xvals,yvals,dataDiff);
axis xy;
fprintf('Please click on right, then left triple points.\n');
[xp,yp] = ginput(2);
spx=abs(xvals-xp(1)) < tsize*dx;
spy=abs(yvals-yp(1)) < tsize*dy;
tuneData.chrg.imgr=dataDiff(spy,spx);
spx=abs(xvals-xp(2)) < tsize*dx;
spy=abs(yvals-yp(2)) < tsize*dy;
tuneData.chrg.imgl=dataDiff(spy,spx);
end