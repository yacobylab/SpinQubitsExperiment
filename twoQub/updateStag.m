function updateStag
% Make stag readout for two qubits.
% Add 0.2 of wait time to offset. 
global tuneData;

name1 = tuneData.alternates(1).activeSetName; 
pd1 = pdload(name1);    
meas1 = pd1.meas; 

meas2.val = [-2,2]; % Fix me to be smart about sepdir. 
meas2.type = 'wait'; 
meas2.time = meas1.time(1)+0.2; 
name2 = tuneData.alternates(2).activeSetName; 
pd2 = pdload(name2);
meas2(2) = pd2.meas; 

meas1(2).type = 'wait'; 
meas1(2).val = [-2,2]; 
meas1(2).time = meas2(2).time+0.2; 
dict = struct('meas',meas1); 
pdsave(sprintf('stag%s',lower(name1(1))),dict);
%fname1 = sprintf('pd_stag%s_last',lower(name1(1))); 
%save(fname1,'dict'); 
%fname2 = sprintf('pd_stag%s_last',lower(name2(1))); 
dict = struct('meas',meas2); 
%save(fname2,'dict');
pdsave(sprintf('stag%s',lower(name2(1))),dict);
end