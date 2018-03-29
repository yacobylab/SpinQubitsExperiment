function [d,file] = loadFiles(filt,opts)
if ~exist('filt','var') || isempty(filt)
    [file,fpath] = uigetfile('','MultiSelect','on');
else
    [file,fpath]= uigetfile(['*' filt '*'],'MultiSelect','on');
end
if ~exist('opts','var')
    opts = '';
end
if ~iscell(file) && (isempty(file) || all(file ==0))
    return;
end
if ~iscell(file)
    file = {file};
end
if isopt(opts,'cell')
    for i = 1:length(file)
        d{i} = load([fpath file{i}]);
    end
else
    for i = 1:length(file)
        d(i) = load([fpath file{i}]);
    end
end