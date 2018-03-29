function btnrect(src,~)
% Lets you draw a rectangle on data and just fit what's inside. 
% function btnrect(src,~) 

  ax = gca; 
  for i = 1:length(ax.Children) 
      isrect(i) = strcmp(ax.Children(i).Type,'rectangle');
  end
  
  delete(ax.Children(isrect));
  xdata = src.Children(1).XData;
  ydata = src.Children(1).YData;
  
  rect= getrect;
  r = rectangle('Parent',ax);
  r.Position = rect;
  
  fig = ancestor(ax,'figure');
  
  mask = xdata>=rect(1) & xdata<=rect(1)+rect(3);
  fitX =  xdata(mask);
  fitY =  ydata(mask);
  p = polyfit(fitX,fitY,1);
  title(sprintf('Slope = %3.3g, Wdth = %3.3g mV',p(1),1e3*rect(3)))
  
end