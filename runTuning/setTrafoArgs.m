function scan = setTrafoArgs(scan,chans) 
for i = 1:length(chans)
    scan.loops(1).trafofn(i).args = smget(chans(i)); 
end
end