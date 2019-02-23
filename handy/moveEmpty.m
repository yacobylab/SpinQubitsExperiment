function moveEmpty(perc,cellNum)
% Remove all the scans without data from files selected (through popup) into the dead_files subfolder. 
% function moveEmpty(perc,cellNum)
% if perc given, remove all with < perc of data. Default is 90%.
% Removes filesRequires that the file contain a cell array named 'data' containing the
% data, and that the uncollected data be NaN. Examines first cell by
% default, but can specify with cellNum

if ~exist('perc','var'), perc = 0.9; end
if ~exist('cellNum','var'), cellNum = 1; end
[file,fpath]=getFiles('sm*.mat');
if ~exist('dead_files','folder') 
    mkdir('dead_files'); 
end
for i = 1:length(file)
    d=load([fpath file{i}]);    
    if ~isfield(d,'data') || (iscell(data) && length(find(isnan(d.data{cellNum})))>=perc*numel(d.data{cellNum}))
        movefile([fpath file{i}],'dead_files/')
        fprintf('Moved %s \n',file{i}(4:end-4));
    end
end
end