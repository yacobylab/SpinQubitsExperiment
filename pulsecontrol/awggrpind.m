function grp = awggrpind(grp)
% Find group index from name of loaded group (unloaded groups return NaN)
% Can give as char, cell array, or number. 
% function grpind = grpind(grp)

global awgdata;

if ischar(grp)
    grp = {grp};
end
if ~isfield(awgdata(1).pulsegroups,'name')
    names={};
else
    names={awgdata(1).pulsegroups.name};
end
if iscell(grp)   
    for i = 1:length(grp) 
        grp{i} = find(strcmp(grp{i}, names),1,'last');
        if isempty(grp{i}) || grp{i}==0
            grp{i} = nan;         
        end           
    end
    grp = cell2mat(grp);
elseif any(grp > length(awgdata(1).pulsegroups))
    error('Group index too large.');
end
end
