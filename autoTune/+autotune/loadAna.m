function [data,scan] = loadAna(file)

if ~isempty(strfind(file,'*'))
    file = uigetfile([tuneData.dir '/' file]);
end    
load([tuneData.dir '/' file],'data','scan'); 
end
