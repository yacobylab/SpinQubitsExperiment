function chgStr = changeConfigGen(d,oldconfig,oldconfigch,file)
% Compare configvals for a set of files.
% function chgStr = changeConfig(d,oldconfig,oldconfigch,file)
% d: struct loaded with file.

if ~isempty(d.configvals)
    configch = d.configch; configvals = d.configvals;
    changedChans = nan(2,length(configch));
    changeList = [];
    for i = 1:length(configch)        
        % Allow for list of configch to change between files
        ind = strcmp(configch{i}, oldconfigch);
        if any(ind)            
            if configvals(i) ~= oldconfig(ind)
                changeList = [changeList i]; % List of change channels
                changedChans(1,i) = configvals(i); % List of new values
                changedChans(2,i) = oldconfig(ind); % List of old values
            end            
        end
    end
    changedChans = configch(changeList);
    
    % Exclude the channels ramped in scan and time. 
    for i = 1:length(d.scan.loops)
        [changedChans, newInds] = setdiff(changedChans, d.scan.loops(i).setchan);
        changeList = changeList(newInds);
    end    
    timeInd = strcmpi(changedChans,'Time');
    changeList(timeInd)=[];
    
    % Write the string.
    if exist('file','var')
        chgStr = sprintf('%s: ', file);
    else 
        chgStr ='';
    end
    for i = 1:length(changeList) 
        chgStr = [chgStr,sprintf('%s: %3.3g to %3.3g. ', configch{changeList(i)}, oldconfig(changeList(i)), configvals(changeList(i)))];
    end
    chgStr = [chgStr,sprintf('\n')]; %#ok<*SPRINTFN>
else
    chgStr =sprintf('%s : No configvals \n',file);
end
end