function rundBz(opts)
global tuneData; 
if ~exist('opts','var'), opts = ''; end 
nloop = 100; nrep = 20; 
side = upper(tuneData.activeSetName(1));

grp = sprintf('dBz_swfb_128_%s',side);
if isopt(opts,'fb')
    scanOpts = 'swfb ampok'; 
    fname = sprintf('dBz%s',side);
else
    scanOpts = 'ampok';
    fname = sprintf('dBznofb%s',side);
end
if ~awgcntrl('ison'), awgcntrl('on start wait err'); end
scan = fConfSeq(grp,{'nloop',nloop,'nrep',nrep, 'datachan',tuneData.dataChan,'opts',scanOpts});
scan = measAmp(scan);

smrun(scan,smnext(fname));
sleep('fast'); 
end