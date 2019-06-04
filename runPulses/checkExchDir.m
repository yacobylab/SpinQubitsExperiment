function checkExchDir
exch1 = linspace(0,-1,10);
global tuneData; 
for i = 1:length(exch1)
    updateExch(struct('exch',[exch1(i),1],'opts','all'));
    %testSep;    
    tuneData.stp.updateGroup; 
    tuneData.stp.run; 
end
end