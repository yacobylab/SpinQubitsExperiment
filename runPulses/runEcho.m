function runEcho(opts,config)
% Set up and run Hahn echo scan. 
% function runEcho(opts,config)
% Tries to be intelligent about the amount of noise in the system. 
% opts: 
%   redo
%   fitJ: Try to measure echo for given range of J instead of eps. Use
%   output from anaExchange for fit fn. 
% See below for config 
global tuneData; global fbdata; 

if ~exist('config','var')
    config = struct;
else
    if iscell(config), config=struct(config{:}); end
end
if ~exist('opts','var'), opts = ''; end

dict = pdload(tuneData.activeSetName);
minTime = dict.meas.time(1)+sum(dict.reload.time)+2*dict.dbzpi.time+.5;
%plsTime = ceil(minTime*4)/4+1;
config = def(config,'npoints',15); npoints = config.npoints; % Number of evo times / scans
config = def(config,'nFiles',12); nFiles = config.nFiles; % Number of scans / eps values. 
config = def(config,'eps',linspace(1.45,0.9,nFiles)); eps = config.eps; % Eps variation
%config = def(config,'evo',[0:256]'-90); evo=config.evo; % Offset time ( from tau/2) 
%config = def(config,'plsLength',plsTime); plsLength = config.plsLength; 
config = def(config,'coef',5415); coef=config.coef; % w/ fitJ
config = def(config,'lev',0.1744); lev = config.lev; % w/ fitJ
config = def(config,'pow',-1.375); pow = config.pow; % w/ fitJ
config = def(config,'noise',180); noise = config.noise; % w/ fitJ % fix what this means... 180
config = def(config,'rng',[60,250]); rng = config.rng; % % w/ fitJ,  range of J to run scans on 
if isopt(opts,'fitJ')
    % for % J = J0 e^(-eps/eps0)
    epsFunc = @(j,j0,eps0) -log(j/j0)*eps0;
    j = linspace(rng(1),rng(2),nFiles);
    eps = epsFunc(j,coef,lev);
    T2func = @(j,n,p) n*j.^p;
    tEnd = 1.5*T2func(j,noise,pow); % Determine appropriate evo time given noise level. 
    tStart = linspace(0.07,0.12,nFiles); % Enough to mostly get past initial decay (fish)
elseif isopt(opts,'log')
    epsMax = log(max(eps)); epsMin = log(min(eps)); 
    eps = exp(linspace(epsMax,epsMin,nFiles));    
    tEnd = 1.6*1./(1:nFiles).^0.55; 
    tStart = 1./(1:nFiles).^0.5 * 0.25;
    %tStart = linspace(0.23,0.1,nFiles);
else
    tEnd = linspace(3^(1/6),0.2^(1/6),12).^6; % what? something wrong here.     
    tStart = linspace(0.12,0.07,nFiles);
end

pg.pulses=22; % Pulse number: check plslist to confirm. FIXME: use a pulse name? 
pg.ctrl='loop pack';
side = upper(tuneData.activeSetName(1));
namepat='RamseyE_%02d_%s';
fname = sprintf('RamseyE%s',side);
pg.chan=[str2double(char(regexp(tuneData.xyChan{1},'\d+','match'))),str2double(char(regexp(tuneData.xyChan{2},'\d+','match')))];
pg.dict={tuneData.activeSetName};
pg.dict={struct('prep','@dbzprep','read','@dbzread','pi','@dbzpi'),pg.dict};
%Parameters: p(1) = total time p(2) = epsilon (mv) ; p(3)=total evo time ; p(4) = dt in nsec
awgrm(fbdata.lastInd,'after'); % Assumes 13 pulses stay on AWG constantly. FIXME and make me smarter. 
awgclear('unused'); % Remove inactive pulses from AWG 
pg.trafofn.func=@skinTraf; pg.trafofn.args=8; % Add trafofn to pulses. 

%pg.varpar=evo;
if isfield(config,'redo') % If something crashes in the middle, restart at the current group. 
    iStart = config.redo;    
else
    iStart =1;
end
if isfield(config,'sing') && ~isempty(config.sing)
    inds = config.sing;
    if isfield(config,'tEnd'), tEnd(inds) = config.tEnd; end
else
    inds = iStart:length(eps); 
end
for i = inds
    echoTimes = linspace(tStart(i),tEnd(i),npoints);
    echoTimes = round(500*echoTimes)/500; 
    jgrpL={};
    evoMax = floor(echoTimes(1)*500)-5;
    %pg.varpar = [-evoMax:evoMax+15]';     
    pg.varpar = (-20:15)'; 
    plsTime = ceil((minTime+tEnd(i))*4)/4;
    for echoTime = echoTimes
        pg.xval = echoTime;
        pg.params=[plsTime eps(i) echoTime 0];
        pg.name=sprintf(namepat,length(jgrpL)+1,side);
        %pg.trafofn.func = @skineffect_trafofn; pg.trafofn.args = 0;
        plsdefgrp(pg);
        jgrpL=[jgrpL pg.name];
    end
    awgadd(jgrpL); % Add groups to AWG 
    awgcntrl('on start wait err raw'); % Turn AWG on 
    scan=fConfSeq(jgrpL,struct('nloop',32,'nrep',120,'datachan',tuneData.dataChan,'opts','swfb'));  % Make scan 
    %try
        smrun(scan,smnext(fname));
    %catch % Generally scans fail due to trigger issue. 
    %    fprintf('Trigger issue \n');
    %end
    sleep
end
end