function fileNames = sortFiles(dirName)
% Sort all files from dirName by date. If no directory given, use current directory. 
% function fileNames = sortFiles(dirName)

if ~exist('dirName','var') || isempty(dirName)
    dirInfo = dir; 
else
    dirInfo = dir(dirName); 
end

dirInfo([dirInfo.isdir])=[]; % Remove directories 
dates = [dirInfo(:).datenum]; 

[~,inds] = sort(dates);
fileNames = {dirInfo(inds).name}; % Cell array of names in order by datenum.
end