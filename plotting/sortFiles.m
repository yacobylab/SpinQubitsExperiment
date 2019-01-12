function fileNames = sortFiles(dirName)
% Sort all files from dirName by date. If no directory given, use current directory. 
% function fileNames = sortFiles(dirName)

if ~exist('dirName','var') || isempty(dirName)
    fileFolder = dir; 
else
    fileFolder = dir(dirName); 
end

dates = [fileFolder(:).datenum]; % you may want to eliminate . and .. first.
[~,inds] = sort(dates);
fileNames = {fileFolder(inds).name}; % Cell array of names in order by datenum.
end