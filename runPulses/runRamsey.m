function [scan,fname]=runRamsey(scanType,config)
% function [scan,fname]=runRamsey(scanType,config)
% scanType: Type of scan to run, default dbz. 
%   dbz
%   redo
%   adPrep
%   dbzCoef
%   single
% config: scan details, see below. 

if ~exist('scanType','var'), scanType = 'dbz'; end
global tuneData; persistent scanOpts;
if ~exist('config','var')
    config = struct;
else
    if iscell(config), config=struct(config{:}); end
end
config = def(config,'opts',''); 
side = upper(tuneData.activeSetName(1));
fname = sprintf('Ramsey%s',side);
dict = pdload(tuneData.activeSetName);
minTime = dict.meas.time(1)+sum(dict.reload.time)+dict.dbzpi.time+0.1;
plsTime = ceil(minTime*4)/4+2;
if ~strcmp(scanType,'redo')
    config = def(config,'npoints',18); npoints = config.npoints;
    config = def(config,'eps',linspace(0.6,3,npoints)); eps = config.eps;
    config = def(config,'evo',1:100); evo=config.evo;
    config = def(config,'plsLength',plsTime); plsLength = config.plsLength;
    config = def(config,'coef',2130); coef=config.coef;
    config = def(config,'lev',0.33); lev = config.lev;
    config = def(config,'opts','');
    awgrm(13,'after'); awgclear('unused'); % Fix me! 
    
    pg.pulses=38;
    pg.varpar = evo';
    pg.chan=[str2double(char(regexp(tuneData.xyChan{1},'\d+','match'))),str2double(char(regexp(tuneData.xyChan{2},'\d+','match')))];
    pg.dict={tuneData.activeSetName};
    pg.ctrl='loop pack';
    %if isopt(config.opts,'traf')
    %    traf = load('Z:/Shannon/Data/imp');
    %end
end
if isopt(config.opts,'trafo')
    %pg.trafofn.func=@kernelTraf; pg.trafofn.args=traf.h2;
    %pg.trafofn.func=@rc_trafofn;         pg.trafofn.args=.5;
    pg.trafofn.func=@skinTraf; pg.trafofn.args=8;
end
switch scanType
    case 'adprep'
        pg.dict={struct('prep',struct('type','@adprep'),'read',struct('type','@adread')),pg.dict};
        for i = 1:length(eps)
            pg.params=[plsLength eps(i) 0]; %Parameters: pulselength, eps, evo
            pg.name = sprintf('Ramsey_%s_%d',side,i);
            plsdefgrp(pg);
            ramsey{i}=pg.name; %#ok<*AGROW>
        end
        awgadd(ramsey);
        awgcntrl('on start wait err');
        scanOpts = '';
    case 'dbz'
        pg.dict={struct('prep',struct('type','@dbzprep'),'read',struct('type','@dbzread')),pg.dict};
        for i = 1:length(eps)
            name{i}= sprintf('Ramsey_L_%d',i);
            pg.params=[plsLength eps(i) 0]; %Parameters: pulselength, eps, evo
            pg.name = name{i};
            plsdefgrp(pg);
            ramsey{i} = pg.name;
        end
        awgadd(ramsey);
        awgcntrl('on start wait err');
        scanOpts = 'swfb';
    case 'dbzCoef'
        pg.dict={struct('prep',struct('type','@dbzprep'),'read',struct('type','@dbzread')),pg.dict};        
        % for % J = J0 e^(-eps/eps0)
        % for now, this is lev=out.decp(2), coef = outdecp(2)
        epsFunc = @(j,j0,eps0) -log(j/j0)*eps0;
        j = linspace(30,450,npoints);
        eps = epsFunc(j,coef,lev);
        for i = 1:length(eps)
            name{i}= sprintf('Ramsey_L_%d',i);
            pg.params=[plsLength eps(i) 0]; %Parameters: pulselength, eps, evo
            pg.name = name{i};
            plsdefgrp(pg);
            ramsey{i} = pg.name;
        end
        awgadd(ramsey);
        awgcntrl('on start wait err');
        scanOpts = 'swfb';
    case 'log'        
        pg.dict={struct('prep',struct('type','@dbzprep'),'read',struct('type','@dbzread')),pg.dict};        
        epsMax = log(max(eps)); epsMin = log(min(eps));
        eps = exp(linspace(epsMax,epsMin,npoints));
        for i = 1:length(eps)
            name{i}= sprintf('Ramsey_L_%d',i);
            pg.params=[plsLength eps(i) 0]; %Parameters: pulselength, eps, evo
            pg.name = name{i};
            plsdefgrp(pg);
            ramsey{i} = pg.name;
        end
        awgadd(ramsey);
        awgcntrl('on start wait err');
        scanOpts = 'swfb';
    case 'single' % Run single ramsey group.
        pg.dict={struct('prep',struct('type','@dbzprep'),'read',struct('type','@dbzread')),pg.dict};
        eps = 1;
        evo = 1:50;
        pg.varpar = evo';
        name= sprintf('Ramsey_%s',side);
        
        pg.params=[plsLength eps 0];
        pg.name = name;
        plsdefgrp(pg);
        awgadd(pg.name);
        awgcntrl('on start wait err raw');
        scanOpts = 'swfb';
        ramsey = {name};
    case 'redo'
        global awgdata;
        plsgrps = {awgdata.pulsegroups.name};
        mask=contains(plsgrps,'ramsey','IgnoreCase',true);
        ramsey = plsgrps(mask);
        scanOpts = 'swfb';
        %scanOpts = '';
end

scan=fConfSeq(ramsey,struct('nloop',60,'nrep',40,'datachan',tuneData.dataChan,'opts',scanOpts));
smrun(scan,smnext(fname)); sleep

end