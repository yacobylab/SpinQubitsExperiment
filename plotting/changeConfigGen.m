function chgStr = changeConfigGen(d,oldconfig,oldconfigch,file)
% Compare configvals for a set of files. 
%function chgStr = changeConfig(d,oldconfig,file)
% d: struct loaded with file. 

chgStr ='';
if ~isempty(d.configvals)
    configch = d.configch; configvals = d.configvals;    
    changedChans = nan(2,length(configch)); 
    changeList = [];    
    for i = 1:length(configch) 
       ind = find(strcmp(configch{i}, oldconfigch));  %#ok<*EFIND>
       if isempty(ind) 
           inds(i) = NaN; 
       else
           if configvals(i)~=oldconfig(ind)
               changeList = [changeList i]; 
               changedChans(1,i) = configvals(i); 
               changedChans(2,i) = oldconfig(i); 
           end
       end       
    end
    changedChans = configch(changeList);         
    [changedChans, newInds] = setdiff(changedChans, d.scan.loops(1).setchan);
    changeList = changeList(newInds);
    [changedChans, newInds] = setdiff(changedChans, d.scan.loops(2).setchan);
    changeList = changeList(newInds);
    timeInd = strcmpi(changedChans,'Time');
    changeList(timeInd)=[];
    for j = 1:length(changeList)
        if j ==1
            chgStr = [chgStr, sprintf('%s: ', file)];
        end
        chgStr = [chgStr,sprintf('%s: %3.3g to %3.3g. ', configch{changeList(j)}, oldconfig(changeList(j)), configvals(changeList(j)))];
        if j == length(changeList)
            chgStr = [chgStr,sprintf('\n')]; %#ok<*SPRINTFN>
        end
    end        
else
    chgStr =sprintf('%s : No configvals or length changed. \n',file);
end
end