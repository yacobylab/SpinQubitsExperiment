function updateFiles 
global pptdata 
pptdata.dir = dir; 
pptdata.qpcDir = dir(pptdata.qpcFolder); 

pptdata.fileNames = sortFiles; % Cell array of names in order by datenum.
pptdata.qpcfileNames = sortFiles(pptdata.qpcFolder); % Cell array of names in order by datenum.
save('Z:/Shannon/Data/pptdata','pptdata'); 
end