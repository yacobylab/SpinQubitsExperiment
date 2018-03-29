function pdsave(name, dict)
% Save a new version of a pulse dictionary.
%function pd=pdsave(namne, dict)
% example: pdsave('right',r), pdsave('left',l) or pdsave('A',a);

global plsdata; global awgdata; 
file=[plsdata.grpdir, 'pd_', name,'_last.mat'];
save(file,'dict');

for i =1:length(awgdata.pulsegroups) 
  %awgdata.pulsegroups(i).changed = true; 
end
% Fixme; there should be codde to find loaded groups that need updating here.   
end