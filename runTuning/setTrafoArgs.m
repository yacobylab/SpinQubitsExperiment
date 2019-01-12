function scan = setTrafoArgs(scan,chans) 
% Set the trafofn args to take current values of given chans. 
% function scan = setTrafoArgs(scan,chans) 
% FIXME: what is this from? 

for i = 1:length(chans)
    scan.loops(1).trafofn(i).args = smget(chans(i)); 
end
end