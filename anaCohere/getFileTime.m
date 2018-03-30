function time = getFileTime(file)
% Returns time file created. 
% function time = getFileTime(file)
    if ~contains(file,'.mat') 
        file = [file '.mat']; 
    end
    fileInfo = dir(file);
    time = fileInfo.datenum;
end