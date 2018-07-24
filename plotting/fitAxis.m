function fitAxis(a)
if ~exist('a','var'), a = gca; end
XData = []; YData = [];
for i =1:length(a.Children)
    XData = [XData,a.Children(i).XData]; 
    YData = [YData,a.Children(i).YData]; 
end
if min(XData)<max(XData)
    a.XLim = [min(XData),max(XData)];
end
if min(YData)<max(YData)
    a.YLim = [min(YData),max(YData)]; 
end
end