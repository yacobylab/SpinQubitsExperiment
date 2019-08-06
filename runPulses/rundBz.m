function rundBz(opts,config)
% Using the 'ave' config 
% config: 
%   nloop: default 48 
%   nrep: default 60 
%   ave: set the total number of averages, 48*60.
%   fb: Use feedback
%   long: spacing of 4 ns, 512 ns max  (default spacing 1ns, max 128 ns)
global tuneData;
if ~exist('config','var')
    config = struct;
else
    if iscell(config), config=struct(config{:}); end
end
if ~exist('opts','var'), opts = ''; end
config = def(config,'nloop',48);
if isfield(config,'ave') && ~isempty(config.ave)
    config.nrep = floor(config.ave/config.nloop);
else
    config = def(config,'nrep',60); 
end

if isopt(opts,'both')
    side = 'LR'; 
    datachan = 'both'; 
    grp = sprintf('dBz_stag_%s',side); 
else
    side = upper(tuneData.activeSetName(1));
    datachan = tuneData.dataChan; 
    if isopt(opts,'long')
        maxTime = 512;
    else
        maxTime = 128;
    end
    grp = sprintf('dBz_swfb_%d_%s',maxTime,side);
end
if isopt(opts,'fb')
    scanOpts = 'swfb ampok';
    fname = sprintf('dBz%s',side);
else
    scanOpts = 'ampok';
    fname = sprintf('dBznofb%s',side);
end
if ~awgcntrl('ison'), awgcntrl('on start wait err'); end
scan = fConfSeq(grp,{'nloop',config.nloop,'nrep',config.nrep, 'datachan',datachan,'opts',scanOpts});
scan = measAmp(scan,opts);

smrun(scan,smnext(fname));
sleep('fast');
end