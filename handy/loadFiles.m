function [d,file] = loadFiles(filt,opts)
% load files, return in struct or cell d with list of names. 
% function [d,file] = loadFiles(filt,opts)
% filt gives filter for file names. 
% opts can be "cell" if the files have different content structures, default is struct. 

if ~exist('filt','var') || isempty(filt)
    [file,fpath] = uigetfile('','MultiSelect','on');
else
    [file,fpath]= uigetfile(['*' filt '*'],'MultiSelect','on');
end
if ~exist('opts','var'), opts = ''; end
if ~iscell(file) && (isempty(file) || all(file ==0))
    return;
end
if ~iscell(file), file = {file}; end
if isopt(opts,'cell')
    for i = 1:length(file)
        d{i} = load([fpath file{i}]);
    end
else
    try % try to load as struct, but load as cell if not possible. 
        for i = 1:length(file)
            d(i) = load([fpath file{i}]);
        end
    catch
        for i = 1:length(file)
            d{i} = load([fpath file{i}]);
        end
    end
end