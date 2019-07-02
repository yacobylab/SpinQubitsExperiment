function b=basislookup(basis)
global tuneData
% Give this either a single gate in string form. 
if ischar(basis)
    basis={basis}; 
end
if iscell(basis)
    for j = 1:length(basis)
        b(j)=find(strncmp(basis{j},tuneData.baseNames,length(basis{j}))); 
    end
else 
    fprintf('Please input a cell or a string')
    b=nan; 
end
end