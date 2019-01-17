function doublePt(src,clk)
% Click on a pt on a color plot and plot point at same x/y on other figure. 
% function doublePt(src,clk)
% Assumes two plots are figures 1 / 3. 
yVal = clk.IntersectionPoint(2);
xVal = clk.IntersectionPoint(1); 

if src.Parent.Parent.Number == 1
    fig = 3; 
else
    fig = 1; 
end
f = gcf;
axesInds = find(isgraphics(f.Children,'axes')); 
for i =1:length(axesInds)
    hold(f.Children(axesInds(i)),'on')
    plot(f.Children(axesInds(i)),xVal,yVal,'.','MarkerSize',15)
end
set(0,'CurrentFigure',fig);
f = gcf;
axesInds = find(isgraphics(f.Children,'axes')); 
for i =1:length(axesInds)
    hold(f.Children(axesInds(i)),'on')
    plot(f.Children(axesInds(i)),xVal,yVal,'.','MarkerSize',15)
end
end