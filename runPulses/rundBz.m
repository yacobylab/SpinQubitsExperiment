function rundBz 
global tuneData; 
nloop = 200; nrep = 20; 
side = upper(tuneData.activeSetName(1));
fname = sprintf('dBz%s',side);
grp = sprintf('dBz_swfb_128_%s',side);
scan = fConfSeq(grp,{'nloop',nloop,'nrep',nrep, 'datachan',tuneData.dataChan,'opts','ampok'});
scan = measAmp(scan);
if ~awgcntrl('ison')
    awgcntrl('on start wait err');
end
smrun(scan,smnext(fname));

sleep('fast'); 
end