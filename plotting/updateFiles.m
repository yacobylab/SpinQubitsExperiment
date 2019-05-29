function updateFiles
% Check if any new files have been created in since the last time plots were
% added to ppt.
global pptdata
pptdata.dir = dir; % Is this necessary? 
pptdata.qpcDir = dir(pptdata.qpcFolder); 
pptdata.tuneDir = dir(pptdata.tuneFolder); 
pptdata.fileNames = sortFiles(pptdata.dataFolder); % Cell array of names in order by datenum.
pptdata.qpcfileNames = sortFiles(pptdata.qpcFolder); % Cell array of names in order by datenum.
pptdata.tuneFileNames = sortFiles(pptdata.tuneFolder); % Cell array of names in order by datenum.
save(pptdata.file,'pptdata');
end