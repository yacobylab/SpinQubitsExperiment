function pg=make1grp(pg,config)
% Right now, this asks that you provide the pulse number, varpar, params,custom
% dictionary
global tuneData; global qdata; 
if iscell(pg), pg=struct(pg{:}); end
if iscell(config), config=struct(config{:}); end

config = def(config,'namePat','grp_%02d_%s'); 
config = def(config,'n',1); 
config = def(config,'side',tuneData.activeSetName); 
config = def(config,'customDict',struct()); 
config = def(config,'opts',''); 
ind = strcmp({qdata.name},config.side);
pg.ctrl='loop pack';
pg.chan = qdata(ind).chan;
pg.name=sprintf(config.namePat,config.n,upper(config.side(1)));
if isopt(config.opts,'stag')
    stagDict = sprintf('stag%s',lower(config.side(1)));
end
if ~isopt(config.opts,'two')
    pg.chan = pg.chan(1:2); 
end
pg.dict={config.customDict,stagDict,config.side};
plsdefgrp(pg);
end