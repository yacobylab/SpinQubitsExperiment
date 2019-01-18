function [d,fileList] = sortConfig(d,fileList,channel,chanval,precis)
% function [d,fileList] = sortConfig(d,fileList,channel,chanval,precis)
% Look for files with a certain range of values. 
val = nan(1,length(d));
for i = 1:length(d)
    chanInd = strcmp(channel,d(i).configch);
    if any(chanInd)
        val(i) = d(i).configvals(chanInd);
    end
end
if ~exist('precis','var')
    precis = 3e-3;
end

if ~exist('chanval','var')
    fprintf('List of possible values for %s:',channel); 
    fprintf('%3.3f. ',unique(val))
    chanval = input('Which one do you want?');
    precis = input('What precision?');
end
inds = abs(val-chanval) < precis;
d = d(inds); 
fileList= fileList(inds); 
end