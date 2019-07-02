function checkExchDir
global tuneData; 
dict = pdload(tuneData.activeSetName); 
exchStart = dict.exch.val; 
exch1 = linspace(0,-3,8);

for i = 1:length(exch1)
    updateExch(struct('exch',[1,exch1(i)],'opts','all'));
    %testSep;    
    tuneData.stp.updateGroup; 
    tuneData.stp.run; 
end
updateExch(struct('exch',exchStart,'opts','all'));
end