function moveEmpty(perc)
% remove all the scans without data from the folder. if perc given, remove
% all with < perc of data. 
% function moveEmpty(perc)
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