function chans=getMultiGateScan(scan) 
% Give scan, return cell returned with which channels are ramped for multi scna.  
% function chans=getMultiGateScan(scan) 
for i = 1:length(scan.loops(1).trafofn) 
    fn = func2str(scan.loops(1).trafofn(i).fn); 
    chanInds = strfind(fn,'==');  
    chanList{i} = [];
    for j = 1:length(chanInds)
        chanList{i}(end+1) = str2double(fn(chanInds(j)+2)); % This gives the set 
    end    
end
for i = 1:scan.loops(2).npoints
    chanNums{i}=find(cellfun(@(x) ~isempty(find(x==i)),chanList));   
    chans{i}='';
    for j = 1:length(chanNums{i})
    chans{i} = [chans{i} scan.loops(1).setchan{chanNums{i}(j)}, ', ']; 
    end 
    if ~isempty(chans{i})
        chans{i}(end-1:end)=[];
    end
end
end