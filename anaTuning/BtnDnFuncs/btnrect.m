function btnrect(src,~)
% Lets you draw a rectangle on 1D plot and fit only what's inside. 
% function btnrect(src,~)
% src is data from figure. 
% For charge sensitivity meas. 

ax = gca;
% Delete current rectangles. 
for i = 1:length(ax.Children)
    isrect(i) = strcmp(ax.Children(i).Type,'rectangle');
end
delete(ax.Children(isrect));

xData = src.Children(1).XData;
yData = src.Children(1).YData;

rect= getrect;
r = rectangle('Parent',ax);
r.Position = rect; 

mask = xData>=rect(1) & xData<=rect(1)+rect(3);
fitX =  xData(mask);
fitY =  yData(mask);
p = polyfit(fitX,fitY,1);
title(sprintf('Slope = %3.3g, Wdth = %3.3g mV',p(1),1e3*rect(3)))
end