function [data,scan,file,time] = loadAna(file) %#ok<*STOUT>
global tuneData
if ~isempty(strfind(file,'*'))
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
        data = data{1};  %#ok<*NODEF>
        fileInfo = dir(fileName); 
        time = fileInfo.datenum; 
    else
        data = []; scan = []; time =[];
    end
end    
end
