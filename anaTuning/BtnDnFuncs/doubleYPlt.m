function doubleYPlt(src,clk,cache)
% Allows you to plot differentiated date and non together -- useful with CB
% + charge sensing. 
%function hyst(src,clk,cache)

figure(48); clf; hold on; 
if clk.Button == 1 % x cut
    yVal=clk.IntersectionPoint(2);
    yData = cache(1).yvals;    
    [~,yInd] = min(abs(yData - yVal));
    xData = cache(1).xvals;
    xDataN = cache(2).xvals;
    yyaxis left 
    h(1) = plot(xData, cache(1).data(yInd,:),'DisplayName',cache(1).fileList(end-7:end-4));          
    title(sprintf('Y = %3.3f', yVal)); hold on;
    xlabel(src.Parent.XLabel.String);
    yyaxis right 
    h(2) = plot(xDataN, cache(2).data(yInd,:),'DisplayName',cache(1).fileList(end-7:end-4));   
else % y cut
    xVal = clk.IntersectionPoint(1);
    for i = 1:length(cache)
        xData = cache(i).xvals;
        [~,xInd] = min(abs(xData - xVal));
        yData = cache(i).yvals;
        h(i) = plot(yData, cache(i).data(:,xInd),'DisplayName',cache(i).fileList(end-7:end-4));
    end
    title(sprintf('X = %3.3f', xVal)); hold on;
    xlabel(src.Parent.YLabel.String);
    legend(h);
end