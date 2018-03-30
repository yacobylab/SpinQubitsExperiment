function moveEmpty(perc)
% remove all the scans without data from the folder into the dead_files subfolder. 
% function moveEmpty(perc)
% if perc given, remove all with < perc of data. Default is 90%.

if ~exist('perc','var')
    perc = 0.9;
end

[file,fpath]=get_files('sm*.mat');
for i = 1:length(file)
    d=load([fpath file{i}]);    
    if ~isfield(d,'data') || length(find(isnan(d.data{1})))>=perc*numel(d.data{1})
        movefile([fpath file{i}],'dead_files/')
        fprintf('Moved %s \n',file{i}(4:end-4));
    end
end
end