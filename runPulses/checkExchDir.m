function checkExchDir
global tuneData; 
dict = pdload(tuneData.activeSetName); 
exchStart = dict.exch.val; 
exch1 = linspace(0,-2.5,6);

for i = 1:length(exch1)
    updateExch(struct('exch',[exch1(i),1],'opts','all'));
    %testSep;    
    tuneData.stp.updateGroup; 
    tuneData.stp.run; 
end
updateExch(struct('exch',exchStart,'opts','all'));
end