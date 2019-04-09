function updateFiles
% Check if any new files have been created in since the last time plots were
% added to ppt.
global pptdata
pptdata.dir = dir; % Is this necessary? 
pptdata.qpcDir = dir(ptptdata.qpcFolder); 

pptdata.fileNames = sortFiles; % Cell array of names in order by datenum.
pptdata.qpcfileNames = sortFiles(pptdata.qpcFolder); % Cell array of names in order by datenum.
save(pptdata.file,'pptdata');
end