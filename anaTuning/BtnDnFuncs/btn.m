function btn(src,clk)
% Plot a line cut of the data. Left click horizontal line, right click vertical.
% function btn(src,clk)
% if d typed while figure active, calculate and print distance between 2 points

f=figure(148); clf; f.Name = 'Linecut';
if clk.Button ==1 % Horizontal
    yVal=clk.IntersectionPoint(2);
    yData = linspace(src.YData(1),src.YData(end),size(src.CData,1));
    [~,yInd] = min(abs(yData - yVal));
    yVal=yData(yInd);
    xData = linspace(src.XData(1),src.XData(end),size(src.CData,2));
    plot(xData, src.CData(yInd,:))
    title(sprintf('Y = %3.3g', yVal)); hold on;
    xlabel(src.Parent.XLabel.String);
else % Vertical
    xVal = clk.IntersectionPoint(1);
    xData = linspace(src.XData(1),src.XData(end),size(src.CData,2));
    [~,xInd] = min(abs(xData - xVal));
    xVal=xData(xInd);
    yData = linspace(src.YData(1),src.YData(end),size(src.CData,1));
    plot(yData, src.CData(:,xInd))
    title(sprintf('X = %3.3g', xVal)); hold on;
    xlabel(src.Parent.YLabel.String);
end

fs = figure(148); fs.WindowKeyPressFcn=@(h,e) keyPressed(e);
end

function keyPressed(e)
switch e.Character
    case 'd'
        m = ginput(2);
        widCB = m(2,1)-m(1,1);
        ax = gca;
        xvals = ax.Children(1).XData;
        [~,xInd1]=min(abs(m(1,1)-xvals));
        [~,xInd2]=min(abs(m(2,1)-xvals));
        yData = ax.Children(1).YData([xInd1,xInd2]);
        fprintf('Distance between points is %3.3f mV, Change in value %3.3g \n',1e3*widCB,yData(1)-yData(2))
end
end