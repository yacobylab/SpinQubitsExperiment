function [blTriple,trTriple] = atCorrelate(scan,data)
% function [xp,yp] = at_correlate(runnumber,tsize)
% run number and tsize can be the data/scan struct or a run number tsize. 
%  guess where the triple points are.
%  used by atChargeana 

global tuneData;
debug = false; 
xvals = linspace(scan.loops(1).rng(1),scan.loops(1).rng(2),size(data,2)); xvals= xvals(2:end)/2+xvals(1:end-1)/2; 
yvals = linspace(scan.loops(2).rng(1),scan.loops(2).rng(2),size(data,1));
[x,y] = meshgrid(xvals,yvals);    
if isempty(tuneData.chrg.imgr) 
    return
end
crossCorr=normxcorr2(tuneData.chrg.imgr,data); % This gives the cross correlation of the charge scan data w/ an image of TR triple point
[~, imax] = max(crossCorr(:));
[yrpeak, xrpeak] = ind2sub(size(crossCorr),imax(1)); % Point of max cross correlation. 
xrpeak=xrpeak-round(size(tuneData.chrg.imgr,2)/2); % Padding added to crossCorr -- see help for normxcorr2. 
yrpeak=yrpeak-round(size(tuneData.chrg.imgr,1)/2);
if debug
    figure(56); clf; subplot(2,2,1);
    imagesc(xvals,yvals,crossCorr);
    axis xy;
end
crossCorr=normxcorr2(tuneData.chrg.imgl,data); %Now left point 
[~, imax] = max(crossCorr(:));
[ylpeak, xlpeak] = ind2sub(size(crossCorr),imax(1));
xlpeak=xlpeak-round(size(tuneData.chrg.imgl,2)/2);
ylpeak=ylpeak-round(size(tuneData.chrg.imgl,1)/2);
if debug
    subplot(2,2,2);
    imagesc(xvals,yvals,crossCorr);
    axis xy;
    subplot(2,2,3);
    imagesc(xvals,yvals,data); axis xy; hold on;
end
if yrpeak < 0 || xrpeak < 0 || ylpeak < 0 || xlpeak < 0 || yrpeak > size(x,1) || ylpeak >size(x,1) || xrpeak > size(x,2) || xlpeak > size(x,2)
    warning('Found triple point off edge');
    blTriple = [0,0]; trTriple = [0,0];
else
    plot(x(yrpeak,xrpeak),y(yrpeak,xrpeak),'kx');
    plot(x(ylpeak,xlpeak),y(ylpeak,xlpeak),'k+');
    blTriple = [xvals(xlpeak),yvals(ylpeak)]; trTriple = [xvals(xrpeak),yvals(yrpeak)];    
end
end

