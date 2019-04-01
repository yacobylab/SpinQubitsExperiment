function files = fileNumbers(nums,dirName)
if ~exist('dirName','var'), dirName = pwd; end
dirData = dir(dirName); 
fileNames = {dirData.name};
files = {};
for i = 1:length(nums)
    inds=contains(fileNames,sprintf('%04d.mat',nums(i)));
    files = [files, fileNames(inds)]; 
end
end