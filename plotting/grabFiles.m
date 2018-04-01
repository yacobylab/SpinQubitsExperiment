function finFiles = grabFiles(fileNames, fpat, cooldowns,IndRng)
% Grab a set of files that match the pattern from given set of cooldowns
% and file numbers. 
% function finFiles = grabFiles(fileNames, fpat, cooldowns,IndRng)
% IndRng is a set of start and end inds for scans.  
% fpat is the pattern in regexp to look at 
% fileNames is the set of possible files to consider. 
% cooldowns: 

stVals = IndRng(:,1); % First column is start, second is end. If we go to current values, end is Nan. 
if length(IndRng) > 1
    endVals = IndRng(:,2);
else 
    endVals = NaN;
end

[fpatNames,fileNumsC]=regexp(fileNames,fpat,'match','tokens');
fileNumsC = [fileNumsC{:}];
fileNums = str2double([fileNumsC{:}]); % set of file indices (because sm numbers files). 
if exist('cooldowns','var') && ~isempty(cooldowns) && any(cooldowns>0) % only consider positive cooldowns (how we differentiate folders)
    stValsSet = cooldowns(cooldowns > 0);
elseif exist('cooldowns','var') && ~isempty(cooldowns) % if no correct ones, don't return anything. 
    stValsSet=[];
    findCell = {};
else % if none given, do all. 
    stValsSet = 1:length(stVals);
end
for j = stValsSet
    if j < length(stVals) || ~any(isnan(endVals))
        indsCell{j} = fileNums >= stVals(j) & fileNums <= endVals(j);
    else
        indsCell{j} = fileNums >= stVals(j);
    end
    findCell{j} = find(indsCell{j});
end
inds = cell2mat(findCell);
fileList = [fpatNames{:}];
finFiles = fileList(inds);
end