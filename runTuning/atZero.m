function good = atZero(~,data,thresh,npoints,count)
% Check if pinch off has occured. 
% function good = atZero(~,data,thresh,npoints,count)
if ndims(data{1}>1) && size(data{1},2)>1 
    dataGood = data{1}(count(2),:); 
    dataGood = dataGood(~isnan(dataGood));
else
    dataGood = data{1}(~isnan(data{1}));
end
if length(dataGood)>npoints
    testData = dataGood(end-npoints+1:end);
    if all(abs(testData)<thresh)
        good = 0;
    else
        good = 1; 
    end
else
    good = 1;
end
end