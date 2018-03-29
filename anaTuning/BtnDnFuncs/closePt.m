function [xfin]=closePt(data,x,xvals)
% Return the closest point on the graph after clicking. 
for i = 1:length(x) 
    [~,ind] = min(abs(data-x(i))); 
    if exist('xvals','var') && ~isempty(xvals)
        xfin(i) = xvals(ind); 
    else 
        xfin(i) = ind; 
    end
end
end