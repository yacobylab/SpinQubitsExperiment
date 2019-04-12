function [data,scan,file,time] = loadAna(file) %#ok<*STOUT>
% Tool to simplify loading autotune data
% function [data,scan,file,time] = loadAna(file)
% If star given in filename (e.g. Zoom*) uses that as filter to load zoom
% data. 
% Otherwise, assums filename is exact
global tuneData
if contains(file,'*')
    file = uigetfile([tuneData.dir '/' file]);
    file = [tuneData.dir '/' file]; 
    load(file,'data','scan'); 
    data = data{1};  %#ok<*NODEF>
    fileInfo = dir(file);
    time = fileInfo.datenum;
else
    fileName = [tuneData.dir '/' file];
    if exist(fileName,'file')
        load(fileName,'data','scan'); 
        if exist('data','var')
            data = data{1};  %#ok<*NODEF>
        else
            data = [];
        end
        fileInfo = dir(fileName); 
        time = fileInfo.datenum; 
    else
        data = []; scan = []; time =[];
    end
end    
end
