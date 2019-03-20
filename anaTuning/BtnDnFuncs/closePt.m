function xCls=closePt(data,yCls,xVals)
% Given a 1D array of data, find the index with value closest to yCls. If you give xVals,
% will return xVal of the index. 
% function xfin=closePt(data,x,xvals)
% Primarily used when clicking on graph, to find nearby data points.
% Can give list of yCls values. 

for i = 1:length(yCls) 
    [~,ind] = min(abs(data-yCls(i))); 
    if exist('xVals','var') && ~isempty(xVals)
        xCls(i) = xVals(ind); 
    else 
        xCls(i) = ind; 
    end
end
end