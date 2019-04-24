function doubleYPltDiff(src,clk)
% Allows you to plot differentiated date and non together -- useful with CB
% + charge sensing. 
%function hyst(src,clk,cache)
fignum = 48; 
figure(fignum); clf; hold on; 
if src.Parent.Parent.Number == 1
    fig = 3;
else
    fig = 1; 
end
yVal=clk.IntersectionPoint(2);
yData = linspace(src.YData(1),src.YData(end),size(src.CData,1));
[~,yInd] = min(abs(yData - yVal));
yVal=yData(yInd);
xData = linspace(src.XData(1),src.XData(end),size(src.CData,2));
yyaxis left
plot(xData, src.CData(yInd,:))
title(sprintf('Y = %3.3g', yVal)); hold on;
xlabel(src.Parent.XLabel.String);

f = figure(fig); 
figure(fignum);
yyaxis right
axesInds = find(isgraphics(f.Children,'axes')); 

otherPlot = f.Children(axesInds(1)).Children;  
plot(otherPlot.XData,otherPlot.CData(yInd,:)); 
end