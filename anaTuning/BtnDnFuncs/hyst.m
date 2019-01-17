function hyst(src,clk,cache)
% Plot all the previous data lines in figure.
% function hyst(src,clk,cache)
% Used with plotChrg, where many color plots examined at once. This lets you see line 
% graphs of the same portion of multiple scans. Usefully for seeing how
% small gate changes move things and how much hysteresis there is (depending on scans being compared). 
% cache is the data from other plots. 

figure(48); clf; hold on;
if clk.Button == 1 % x cut
    yVal=clk.IntersectionPoint(2);
    for i = 1:length(cache)
        yData = cache(i).yvals;
        [~,yInd] = min(abs(yData - yVal));
        xData = cache(i).xvals;
        h(i) = plot(xData, cache(i).data(yInd,:),'DisplayName',cache(i).file(end-7:end-4));
    end
    title(sprintf('Y = %3.3f', yVal)); hold on;
    xlabel(src.Parent.XLabel.String);
    legend(h);
else % y cut
    xVal = clk.IntersectionPoint(1);
    for i = 1:length(cache)
        xData = cache(i).xvals;
        [~,xInd] = min(abs(xData - xVal));
        yData = cache(i).yvals;
        h(i) = plot(yData, cache(i).data(:,xInd),'DisplayName',cache(i).file(end-7:end-4));
    end
    title(sprintf('X = %3.3f', xVal)); hold on;
    xlabel(src.Parent.YLabel.String);
    legend(h);
end