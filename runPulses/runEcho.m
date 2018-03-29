function runEcho(opts,config)
global tuneData
% Ramsey-echo

if ~exist('config','var')
    config = struct;
else
    if iscell(config), config=struct(config{:}); end
end
if ~exist('opts','var'), opts = ''; end 

dict = pdload(tuneData.activeSetName); 
minTime = dict.meas.time(1)+sum(dict.reload.time)+dict.dbzpi.time+1.5;
plsTime = ceil(minTime*4)/4+0.25;
config = def(config,'npoints',18); npoints = config.npoints;
config = def(config,'nFiles',20); nFiles = config.nFiles;
config = def(config,'eps',linspace(0.9,0.5,nFiles)); eps = config.eps;
config = def(config,'evo',[0:64]'-25); evo=config.evo;
config = def(config,'plsLength',plsTime); plsLength = config.plsLength;
config = def(config,'coef',5415); coef=config.coef;
config = def(config,'lev',0.1744); lev = config.lev;
config = def(config,'pow',-1.375); pow = config.pow; 
config = def(config,'noise',180); noise = config.noise; 
config = def(config,'rng',[60,250]); rng = config.rng; 
if isopt(opts,'fitJ')
    % for % J = J0 e^(-eps/eps0)
    epsFunc = @(j,j0,eps0) -log(j/j0)*eps0;
    j = linspace(rng(1),rng(2),nFiles);
    eps = epsFunc(j,coef,lev);    
    T2func = @(j,n,p) n*j.^p; 
    tEnd = 1.5*T2func(j,noise,pow);
    tStart = linspace(0.07,0.12,nFiles); 
else
    tEnd = linspace(1.2,0.3,nFiles);
    tStart = linspace(0.07,0.12,nFiles); 
end

pg.pulses=22;
pg.ctrl='loop pack';
side = upper(tuneData.activeSetName(1));
namepat='RamseyE_%02d_%s'; 
fname = sprintf('RamseyE%s',side); 
pg.chan=[str2double(char(regexp(tuneData.xyChan{1},'\d+','match'))),str2double(char(regexp(tuneData.xyChan{2},'\d+','match')))];
pg.dict={tuneData.activeSetName};
pg.dict={struct('prep','@dbzprep','read','@dbzread','pi','@dbzpi'),pg.dict};
%Parameters: p(1) = total time p(2) = epsilon (mv) ; p(3)=total evo time ; p(4) = dt in nsec
awgrm(13,'after'); awgclear('unused');
pg.trafofn.func=@skinTraf; pg.trafofn.args=6.5;            

pg.varpar=evo;
if isopt(opts,'redo')
    iStart = config.redo; 
else
    iStart =1;
end
for i = iStart : length(eps)
    echoTimes = linspace(tStart(i),tEnd(i),npoints); 
    jgrpL={};
    for echoTime = echoTimes
        pg.xval = echoTime; 
        pg.params=[plsLength eps(i) echoTime 0];        
        pg.name=sprintf(namepat,length(jgrpL)+1,side);
        pg.trafofn.func = @skineffect_trafofn; pg.trafofn.args = 0;
        plsdefgrp(pg);
        jgrpL=[jgrpL pg.name];
    end    
    awgadd(jgrpL);
    awgcntrl('on start wait err raw');    
    scan=fConfSeq(jgrpL,struct('nloop',60,'nrep',65,'datachan',tuneData.dataChan,'opts','swfb'));
    try 
        smrun(scan,smnext(fname));
    catch 
        fprintf('Trigger issue \n'); 
    end
    sleep
end
end