function btnFit(src,clk)
% Fit the line you click on to a line. 
% function btnFit(src,clk)

figure(148); clf; 
if clk.Button ==1
    yVal=clk.IntersectionPoint(2);
    yData = linspace(src.YData(1),src.YData(end),size(src.CData,1));
    [~,yInd] = min(abs(yData - yVal));
    yVal=yData(yInd);
    xData = linspace(src.XData(1),src.XData(end),size(src.CData,2));
    plot(xData, src.CData(yInd,:))
    title(sprintf('Y = %3.3f', yVal)); hold on;
    a = gca; a.ButtonDownFcn=@btnrect; 
p = polyfit(xData, src.CData(yInd,:), 1); 
title(sprintf('slope = %3.3f', p(1))); hold on; 
xlabel(src.Parent.XLabel.String);
else
    xVal = clk.IntersectionPoint(1);
    xData = linspace(src.XData(1),src.XData(end),size(src.CData,2));
    [~,xInd] = min(abs(xData - xVal));
    xVal=xData(xInd);
    yData = linspace(src.YData(1),src.YData(end),size(src.CData,1));
    plot(yData, src.CData(:,xInd))
    title(sprintf('X = %3.3f', xVal)); hold on;    
    
    a = gca; a.ButtonDownFcn=@btnrect; 
    p = polyfit(yData, src.CData(:,xInd), 1); 
    title(sprintf('slope = %3.3f', p(1))); hold on; 
    xlabel(src.Parent.XLabel.String);
end


end