function swap(newActiveSet)
% swap the currect active set for the set named newActiveSet (a string)
%function swap(newActiveSet)

global tuneData;
setNames = {tuneData.alternates.activeSetName};
newInd = find(strcmp(setNames,newActiveSet));
stashInd = find(strcmp(setNames,tuneData.activeSetName));

if length(newInd)~=1
    error('Cannot find correct set %s in alternate set names',newActiveSet);
end

if length(stashInd)~=1
    error('Cannot find current active set (%s) in alternates',tuneData.activeSetName);
end
% make a deep copy of tuneData, make swapCopy of alternates
% modify alternates with current info (and delete alternates.alternates)
tmpData = tuneData.copy();
alts = tuneData.alternates.swapCopy();
alts(stashInd) = tuneData.swapCopy();
tuneData.alternates = autotune.Data; %initialize to empty
alts(stashInd) = tuneData;

% now re-init tunedata with deep copy of correct alternate
tuneData = alts(newInd).copy();

%figure out which props to copy over to tunedata from the original copy
mco = ?autotune.Data; %get propertylist data
copyProps = autotune.Data.doNotSwap;
depProps = {mco.PropertyList([mco.PropertyList.Dependent]).Name};
copyProps = setdiff(copyProps,depProps); % do not copy dependent props
for j= 1:length(copyProps)
   tuneData.(copyProps{j}) = tmpData.(copyProps{j}); 
end
tuneData.alternates = alts;

%fprintf('Stashing set %s\n',tuneData.alternates(stashInd).activeSetName);
%fprintf('New active set: %s\n',newActiveSet);
end