function fitAxis(a)
% Set plot limits to be the same as data limits. This is particularly
% helpful for log plots, which often have a lot of empty space. Takes axis
% handle as argument, but if not given uses current axis. 
%function fitAxis(a)

if ~exist('a','var'), a = gca; end

for k = 1:length(a)
    XData = []; YData = [];
    for i =1:length(a(k).Children)
        XData = [XData,a(k).Children(i).XData];
        YData = [YData,a(k).Children(i).YData];
    end
    % Next lines protect for case of no data.
    if min(XData)<max(XData)
        a(k).XLim = [min(XData),max(XData)];
    end
    if min(YData)<max(YData)
        a(k).YLim = [min(YData),max(YData)];
    end
end
end