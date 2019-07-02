%tuneData.stp.search.range = 3000;
tuneData.stp.target = 800; 
tuneData.stp.updateGroup;
%%
profile on; 
scan = fConfSeq(tuneData.stp.plsGrp,{'nloop',2000,'nrep',100, 'datachan',tuneData.dataChan,'opts','ampok'});
scan = measAmp(scan);
side = upper(tuneData.activeSetName(1)); 
scan.loops(1).rng = [0.7 0.6]; 
scan.loops(1).setchan = {'Bz'}; 
%scan.loops(1).ramptime =; 
scan.loops(1).npoints = 50; 

smrun(scan,smnext('funnelR')); 
profile viewer; 