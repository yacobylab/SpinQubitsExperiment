function betterCbar(f,nFile,ncol,opts) 
% Takes in a figure with subplots and determines how many colorbars are
% necessary. Tries to make them small, quickly. 
% function betterCbar(n,nFile,ncol)  
if ~exist('opts','var'), opts = ''; end 
cBarPos = [0.03 0.015 0.01 0.008,0.008,0.007]; 
pos = [0.85 0.375 0.208 0.2]; %85/94
if ~isgraphics(f,'figure')  
    set(0,'CurrentFigure',f); f = gcf;
end
lims = [f.Children.CLim]; minLims = lims(1:2:end); maxLims=lims(2:2:end); % Find ranges of all plots in figure. 
badAxes = minLims==0 & maxLims ==1; % This is the default for an empty figure. 
minLims(badAxes)=NaN; maxLims(badAxes)=NaN; % Ignore the limits of empty figures. 
% get some stats on the limits. 
minMed = nanmedian(minLims); maxMed = nanmedian(maxLims);
diffMax = abs(maxLims - maxMed); diffMin = abs(minLims - minMed);
rangeClim = min(maxLims-minLims);
minInds = diffMin > rangeClim / 2; maxInds = diffMax > rangeClim/2;
badLim = minInds | maxInds;
if length(find(badLim))+length(find(badAxes)) >= length(f.Children)/2 || ~isopt(opts,'skip')
    nAx = length(find(isgraphics(f.Children,'axes'))); 
    for i = 1:nAx
        axesInds = find(isgraphics(f.Children,'axes')); 
        a=f.Children(axesInds(i));  axPos = a.Position;        
        colorbar(a,'Position',[axPos(1)+axPos(3)-cBarPos(ncol), axPos(2), 0.015, axPos(4)], 'Location','manual');
        a.Position(3) = axPos(3)-cBarPos(ncol)-0.008;
    end
else
    goodLim = find(~badLim);
    badLim = find(badLim);
    minLims(badLim)=[]; maxLims(badLim)=[];
    clim = [min(minLims), max(maxLims)];
    [f.Children(goodLim).CLim]=deal(clim);
    a=f.Children(goodLim(1));
    colorbar(a, 'Limits',clim, 'Position',[0.9300 0.2500 0.0250 0.4500]);    
    for i = 1:length(badLim)
        axesInds = find(isgraphics(f.Children,'axes'));
        a=f.Children(axesInds(badLim(i)));  axPos = a.Position;
        colorbar(a,'Position',[axPos(1) + 0.218, axPos(2), 0.015, axPos(4)], 'Location','manual');
        a.Position(3) = 0.212;
    end
end
end
% 16 x 14 x 19