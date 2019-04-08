function [pars,sdata,figs]=anaExchange(file,config)
% Analyze ramsey and echo noise, performing basic fits.
% function [pars, sdata, figs] = anaExchange(file,config)
% Typical 
% config is a config struct describing how to analyze data. Most of it can
% be automatically generated from scan file, and adding config is used only
% if you need to overwrite the default. 
% The file name is cached in a persistent variable to speed up reanalysis
% of data. 
% Possible fields in config:
%   Often used: 
%   rng: Which evolution times to fit. [1 10] means first 10 xvals.  [5 inf] means 5 onwared.
%   grng: Which groups to fit (e.g. epsilon for ramsey data, tau for echo data).
%       [0.5 inf] means .5 onward.  Uses group xvals
%   frames: reps to fit in the scan. (5:20) means fit reps 5-20.
%
%   Rarely used (for unusual scans or when something's not working automatically):
%   xvals: Represents value that varies across groups. For ramsey, epsilon, for echo, tau (total evolution time)
%          Given as a vector, ie. 0.05:0.05:2. Try to guess from group xvals.
%   channel: A list of channels to fit, (in order data collected in
%       scan). Default is first only
%   t1: Ratio of t1 to tmeas used for scaling data, will be taken from most recent t1 scan
%       if not given. 
%   dbz: a list of which channels have dbz data. defaults to any pulsegroup w/ dbz in the name.
%   side: Which set in tunedata to look at.
%   djde: djde in MHz/mV for computing epsilon noise from J noise. Used for calculating echo noise. 
%   ampCutoff: Cutoff amplitude when fitting data; larger amplitudes indicate bad fit. 
%   fitopts: 
%   offset: How much to offset each data set (depends on data amplitude). Default 0.25
%   
% opts: general set of options that can have the following words
%   types of analysis / plots to make: 
%   amp: Plots amplitude vs. epsilon
%   echocenter: plots center time of echo curve vs evo time
%   echophase: plots phase of echo curve vs evo time
%   freq: Plots freq vs xval (epsilon or tau)
%   nofull: Don't plot the full data set (unaveraged).
%   ramseyT2: plots the ramsey t2 time
%   period: Plots period vs epsilon
%   epsrms or epsRMS: calculate rms epsilon noise.
%   echonoise: computes high frequency noise based on echo decay using the specified value of djde
%   color: show color plot of averaged data in its own figure. 
%
% How to process data: 
%   noscale / linescale. Changes default histogramming behavior. no scale doesn't histogram, linescale does ech line separately.
%   guessxval: tries to guess the xvals (times) from the group def.
%   even /odd fit all even or odd groups only.
%   singlegroup: only one pulsegroup in scan, so average them together.
% How to perform analysis:
%   afitdecay: fit decay only if amplitude is bigger than 5%
%   lines: connect the data with lines (not fitted).
%   fitdecay: fits for the decay envelope
%   power: fits the decay exponent
%
%   gauss
%   linfit / logfit: for fitting amplitude
%   mean
%   color / nocolor
%   nocenter: constrains the center of the decay to be at zero time. Can be used for echo 
%   nodbz: No dbz group included in the scan
%   nofitdbz: do not fit the dbz group
%   nofit: don't perform fit of data.
%   noplot: don't plot data
%   noppt
%   noref: don't remove the dBz measurement when calculating frequency. 
%   ramsey: single evo ramsey
%   ramseyQ
%   fitoffset: allow an offset to be fit with J

%% Load data or cached information
figs=[]; pars = []; plotnum=1; fitdescr=''; cache=0;
persistent lastname; persistent sdata_cache; persistent fileinfo;
if ~exist('config','var'), config=struct(); end
if iscell(config), config=struct(config{:}); end
if ischar(config), config=struct('opts',config); end
config = def(config,'opts','ramsey');

if ~exist('file','var') || isempty(file)
    if ~ischar(lastname), lastname=''; end
    if ~isopt(config.opts,'nofilt')
        file=uigetfile('sm*Ramsey*.mat');
    else
        file=uigetfile('sm*.mat','anaExchange',lastname);
    end
end
if isempty(file), return; end
pars.file = file; 
if strcmp(lastname,file) && ~isempty(sdata_cache) && ~isempty(fileinfo) % if p
    st=dir(file);
    if st.bytes == fileinfo.bytes && st.datenum == fileinfo.datenum, cache=1; end
    fileinfo=st;
end
lastname=file;

if cache % Load the data file
    sdata=sdata_cache;
else
    sdata=load(file);  % Load the scan, data.. This lets us auto-generate some options.
    fileinfo=dir(file);
    sdata_cache=sdata;
end
try
    data=sdata.data; scan=sdata.scan;
catch
    return
end
%% Configure defaults, process scan data
config = def(config,'rng',[]); % Range of points (usually evolution times) to fit, i.e. [10 inf];
config = def(config,'grng',[]); % Range of group xvals (usually eps or tau) to fit, ie [1.1 1.5];
config = def(config,'channel',1); % Range of channels to fit. fixme to be smarter
config = def(config,'frames',1:size(data{config.channel},1)); % Default to fit all frames
config = def(config,'dxlabel','T (ns)');  % Default xlabel for intra-group series
config = def(config,'fitopts','interp');  % Interpolate xvals for fits
config = def(config,'fignum',100); % Base number of figures to output.
config = def(config,'spsize',[3, 3]); % Number of subplots for 'minor' figures.
config = def(config,'dbz', find(cellfun(@(p) ~isempty(p),regexp({scan.data.pulsegroups.name},'[dD][bB][zZ]')))); % Is there a dBz group?
config = def(config,'ampCutoff',1.5); % Cutoff amplitude
ha=[];

if isopt(config.opts,'ramsey') % Sets up default options, xlabel for ramsey.
    if isopt(config.opts, 'singlegroup')
        config.opts = [config.opts ' guessxval fitdecay nocenter gauss color nodbz'];
    else
        config.opts=[config.opts ' guessxval freq amp period color nodbz gauss fitdecay nocenter ramseyQ ramseyT2 epsRMS'];
    end
    config = def(config,'xlabel','\epsilon (mV)');  % Default xlabel for inter-group series
else % Assume echo
    config = def(config,'xlabel','T (\mus)'); % Default xlabel for inter-group series
    config.opts=[config.opts ' guessxval freq amp period color nodbz lines'];
    config = def(config,'offset',0.23);
end
if ~isfield(config,'side') || isempty(config.side) % Use DAQ channel to determine side.
    switch scan.loops(1).getchan{1}
        case 'DAQ2'
            config.side='right';
        case 'DAQ1'
            config.side='left';
        otherwise
            error('Unable to determine side');
    end
end

pars.scantime=getFileTime(file); 
if ~isfield(config,'t1') || isempty(config.t1) % Find t1 before rescaling data.
    [~, config.t1] = att1(config.side,pars.scantime,'before',scan);
end 
exchangeData=setdiff(1:length(scan.data.pulsegroups),config.dbz); % Set of groups other than dBz (i.e.interesting data).
ngrps=length(scan.data.pulsegroups);

% Find set of times within each group, and any other varying parameters in
% case those are xvals. 
if ngrps==1 
    ngrps=scan.data.conf.nrep; % Each rep fitted separately.
    exchangeData=(1:ngrps);
    dt = scan.data.pulsegroups(1).varpar'; % time between points
    xvals = scanRng(scan,1); % Assume that sweep is along first loop of scan.
    config.xvals = xvals;
else 
    for i=1:ngrps
        xvt=scan.data.pulsegroups(i).varpar';
        dxvt=diff(xvt,[],2) ~= 0; % Find part of varpar that varies. 
        [~, ind]=max(sum(dxvt,2)); % Take the column with most points varying? 
        dt(i,:)=xvt(ind,:); %#ok<*AGROW>
        if ismember(i,exchangeData), xv(:,i)=xvt(:); end
    end
end
config = def(config, 'dt', dt); % dts
dt = config.dt; % Hack to keep backward compatibility
if any(size(dt)==1), dt = repmat(dt,size(data{1},2),1); end

% Figure out parameter that varies across groups. If guessxval, looks at
% which pulse params vary. Otherwise, use pulsegroups' varpar. 
if ~isfield(config,'xvals') || isempty(config.xvals)
    if isopt(config.opts,'guessxval') % Guess the group xval from the params
        npars=length(scan.data.pulsegroups(1).params);
        pulseparams = [scan.data.pulsegroups.params];
        pulseparams = pulseparams';
        pulseparams = reshape(pulseparams,npars,length(pulseparams)/npars);
        if size(pulseparams,2) > 1 % Multiple pulse params
            varPar = find(diff(pulseparams(:,exchangeData),[],2) ~= 0); %Find the param that changes
            if isnan(mode(varPar))
                xvals = 1:length(exchangeData);
                warning('Can''t find varying parameter among groups'); 
            else
                xvals = pulseparams(mode(varPar),:)';
            end
        else
            xvals = pulseparams';
        end
    else
        varPar = find(diff(xv,[],2) ~= 0);
        if ~isempty(varPar)
            xvals = xv(mode(varPar),:)';
        else
            fprintf('No xval changes from group to group.  Try guessxval\n');
            xvals = xv(1,:);
        end
    end
else
    xvals=config.xvals;
    if (length(xvals) ~= ngrps), error('The length of xvals (%d) must be the same as the number of groups (%d)\n',length(xvals),length(scan.data.pulsegroups)); end
end
xvals=xvals(:);
if isopt(config.opts,'even') % only ana odd / evengroups
    exchangeData=exchangeData(2:2:end); 
elseif isopt(config.opts,'odd')
    exchangeData=exchangeData(1:2:end);
end

% Scale the data from voltages to binary.
if isopt(config.opts,'linescale') % Scale data line by line.
    t1=ones(config.channel,1).*config.t1; % will break when anaExchange can do more than one channel at once.
    dataAll=anaHistScaleLine(scan,data,t1);
    ylab='P(T)';
    if isfield(config,'offset')
        offset=config.offset;
    else
        offset=1/3;
    end
elseif ~isopt(config.opts,'noscale') % Scale the data by histogramming.    
    dataAll=anaHistScale(scan,data,config.t1);
    ylab='P(T)';
    if isfield(config,'offset')
        offset=config.offset;
    else
        offset=1/2;
    end
else
    dataAll=data;
    ylab='V_{rf} (mV)';
    offset=4e-4;
end

if isfield(scan.data,'setpt') % Use set point for dbz stored in scan
    dbzFreq = scan.data.setpt(1)*2*pi/1000;
else
    dbzFreq = 0.0315 * 2 * pi;
end
%% Fitting and plotting.
for i=1:length(config.channel) 
    data=dataAll{config.channel(i)};
    fignum=config.fignum+100*(i-1);
    % When to create new figures
    newFig = 1; 
    if size(dt,1)>20 % If > 20 groups, make two plots and plot half of data on each. 
        newFig = [newFig floor(size(dt,1)/2)];  
    end
    params=[]; currOff=0; % initialize offset
    
    for j=1:length(exchangeData) % Fit and plot all the nondbz data.
        ind=exchangeData(j);
        if isopt(config.opts, 'singlegroup') % If so, we'll average all the data together.
            currData = squeeze(nanmean(data(config.frames,:)));
        else
            currData=squeeze(nanmean(data(config.frames,ind,:),1))'+currOff;
        end
        currOff=currOff+offset;
        if ~isopt(config.opts,'noplot')
            if j==newFig % Prepare figure
                figure(fignum); fignum = fignum + 1; 
                figs=[figs gcf]; clf; 
                a = gca; 
                ylabel(ylab); 
                xlabel(config.dxlabel); hold on; 
            end                 
            if isopt(config.opts,'lines') % Connect data with lines.
                plot(dt(ind,:),currData,'.-');
            else
                plot(dt(ind,:),currData,'.');                                
            end            
        end
        if ~isopt(config.opts,'nofit')
            if ~isopt(config.opts,'noplot')
                colorOrd = a.ColorOrderIndex;
                if colorOrd>1
                    colorOrd = colorOrd-1;
                else
                    colorOrd = 7;
                end
            end
            a.ColorOrderIndex = colorOrd;                        
            params(j,:)=fitosc(dt(ind,:),currData,config.opts,config.rng,'-');
        end
    end
    if isopt(config.opts,'nofit'), return; end
    if length(config.channel)>1 
        pars.paramsCell{i} = params; pars.dtCell{i}=dt;
    end
    pars.params=params; 
    pars.xvals(i,:)=xvals'; pars.dt=dt; pars.dt1(i,:)=dt(1,:);
    if isopt(config.opts,'echo')
        if length(params) > 3
            if size(params,1) > 10
                meanJ=mean(params(3:end-3,4));
            else
                meanJ=mean(params(:,4));
            end
        else
            meanJ=nan;
        end
        title(sprintf('|J|=%g MHz',1e3*meanJ/(2*pi)));
        fitdescr = [fitdescr sprintf('|J|=%g MHz\n',1e3*meanJ/(2*pi))];
    end % Data has single value of J
    if ~isopt(config.opts,'nodbz')
        for j=1:length(config.dbz)
            dbzdata=squeeze(nanmean(data(config.frames,config.dbz(j),:),1))';
            [plotnum,figs,ha] = nextfig(config,plotnum,fignum,figs,ha);
            plot(dt(1,:),dbzdata,'.');
            xlabel('T (ns)'); ylabel(ylab);
            if ~isopt(config.opts,'nofitdbz')
                fp=fitosc(dt(1,:),dbzdata,['fitdecay nocenter plot ' config.fitopts],[]);
                hold on;
                str=sprintf('T_2^*=%.3g ns, V=%.3f, T=%.3f, phi=%f',1./fp(6),2*sqrt(fp(2)^2+fp(3)^2),2*pi/fp(4),atan2(fp(3),fp(2))-pi);
                title(str);
                pars.dbzt2=1./fp(6);
                fitdescr = [ fitdescr sprintf('dBz_%d: ',config.dbz(j)) str newline ];
                dbzFreq = fp(4);
                pars.omega_dbz=dbzFreq;
            else
                title('dBz reference');
            end
        end
    end %Plot dbz check
    if size(xvals,2)>1, xvals=xvals(exchangeData); end    
    if ~isempty(config.grng)
        ind = find(xvals > config.grng(1) & xvals < config.grng(2));
    else
        ind=1:length(xvals);
    end
    if isopt(config.opts,'amp') % Plot amplitude as a function of eps, fit T2.
        [plotnum,figs,ha] = nextfig(config,plotnum,fignum,figs,ha);
        ampFunc=@(x) 2*(sqrt(x(:, 2).^2 + x(:, 3).^2));
        amp=ampFunc(params);
        if ~isnan(config.ampCutoff) % Only include data with reasonable amplitude.
            mask=amp < config.ampCutoff;
        else
            mask=~isnan(amp);
        end
        plot(xvals(mask),ampFunc(params(mask,:)),'.');
        ylabel(ylab); xlabel(config.xlabel);
        if isopt(config.opts,'echo')
                indAmp = ind & mask;
            if isopt(config.opts,'linfit')            
                fpp = fitwrap('plfit',xvals(indAmp)', ampFunc(params(indAmp,:))', [-1 0], @(p,x) p(1)*x + p(2));
                str = [str sprintf('Amp %g t + %g',fpp(1),fpp(2))];
                fp=[]; fp(1) = fpp(2); fp(2) = (fpp(2)+fpp(1)*mean(xvals(indAmp)))/fpp(1);
            elseif isopt(config.opts,'logfit')                
                fpp = fitwrap('plfit',xvals(indAmp)', log(ampFunc(params(indAmp,:)))', [-1 0], @(p,x) p(1)*x + p(2));
                str = [str sprintf('Amp exp^(%g t + %g)',fpp(1),fpp(2))];
                fp=[]; fp(1) = fpp(2); fp(2) = 1/fpp(1);
            else
                [fp,~,str]=fitdecay(xvals(mask)',ampFunc(params(mask,:))',['plot ' config.opts],config.grng);
                str=['Amp' str];
                str= [str sprintf('Amp=%.3f T_2^{echo}=%.3g Q=%3.1g T=%.3g',fp(1),fp(2),fp(2)/(1e-3*2*pi/meanJ),2*pi/mean(meanJ))];
            end
            fitdescr = [ fitdescr 'Amp: ' str newline ];            
            title(str); 
            pars.maxAmp=max(ampFunc(params(ind,:)));
            pars.T2 = fp(2); pars.Q = fp(2)/(1e-3*2*pi/meanJ);
            pars.J = meanJ;
            pars.T= 2*pi/mean(meanJ); pars.afp = fp; pars.amp = fp(1);
            pars.ampDecay = ampFunc(params(mask,:))';
        else
            title('Amplitude')
        end
    end
    if isopt(config.opts,'period') % Plot period as a function of eps
        [plotnum,figs,ha] = nextfig(config,plotnum,fignum,figs,ha);
        periodFunc=@(params) 2*pi./params(:,4);
        if ~isopt(config.opts,'amp'), mask = 1:size(params,1); end
        mask2 = abs(periodFunc(params))<1e4; % Can assume these are bad fits. 
        maskPer = mask & mask2;
        plot(xvals(maskPer),periodFunc(params(maskPer,:)),'.');
        title('Period'); ylabel('T (ns)'); xlabel(config.xlabel);
    end
    if isopt(config.opts,'freq') % Plot J vs. epsilon
        if ~isopt(config.opts,'noref')
            freqFunc=@(params) abs(sqrt( abs(params(:,4)./(2*pi)).^2-((dbzFreq/((2*pi)))^2)));
        else
            freqFunc=@(params) params(:,4)./(2*pi);
        end
        if isopt(config.opts,'ramsey')
        aliasFreq = 0.5; % in units of GHz
        aliasedData = find(diff(freqFunc(params))<0); % % Find when frequency decreases
        params(aliasedData+1,4)=-params(aliasedData+1,4)+aliasFreq*4*pi;
        end
        [plotnum,figs,ha] = nextfig(config,plotnum,fignum,figs,ha);
        plot(xvals,1e3*real(freqFunc(params)),'.');
        ylabel('J (MHz)'); xlabel(config.xlabel);
        if isopt(config.opts,'ramsey')
            if isopt(config.opts,'fitoffset')
                jFit = fitdecay(xvals',1e3*freqFunc(params)','plot fitoffset',config.grng);
            else
                jFit = fitdecay(xvals',1e3*freqFunc(params)','plot',config.grng);
            end
            str=sprintf('Decay const = %2.2f', jFit(2)); title(str);
            fitdescr = [ fitdescr 'Freq: ' str newline ];
            fitdescr = [fitdescr sprintf('J(eps) = 100 * exp(-(eps-%.4g)/%.3g) + %.3g MHz \n',log(jFit(1)/100)*jFit(2),jFit(2),jFit(3))] ;
            pars.freqfunc=@(eps) 100*exp(-(eps-log(jFit(1)/100)*jFit(2))/jFit(2))+jFit(3);
            pars.jFit = jFit;
            pars.eps=xvals';
        elseif isopt(config.opts,'echo')
            pars.tau=xvals';
        end
        pars.freq = 1e3*freqFunc(params)';
    end
    if isopt(config.opts,'ramseyT2') % Plot T2 vs. epsilon (plot all of these).
        [plotnum,figs,ha] = nextfig(config,plotnum,fignum,figs,ha);
        t2s = abs(1./params(:,6)); pars.t2s = t2s';
        mask = mask & t2s < 1000;
        plot(xvals(mask),t2s(mask),'.-');
        xlabel(config.xlabel); ylabel('T_2^* (ns)');
        [plotnum,figs,ha] = nextfig(config,plotnum,fignum,figs,ha);
        if ~isopt(config.opts,'amp'), mask = 1:size(params,1); end
        plot(1e3*real(freqFunc(params(mask,:))),t2s(mask),'.-');
        xlabel('J (MHz)'); ylabel('T_2^* (ns)');
    end
    if isopt(config.opts,'epsRMS') % Plot djdEps vs T2*, fit slope for low freq. noise
        djde = (1/jFit(2))*(1e3*freqFunc(params)-jFit(3)); %Offset does not contribute to slope
        t2s = abs(1./params(ind,6))';
        noiseSlope = (1./djde(ind))'/t2s; % this is matlab shorthand for least squares slope
        [plotnum,figs,ha] = nextfig(config,plotnum,fignum,figs,ha);
        epsrms=sqrt(2)*noiseSlope/(2*pi); pars.epsrms = epsrms;
        pars.djde = djde';
        % dJ/dE is in MHz/mV=GHz/V, t2* is in ns, so erms in v.        
        plot(1./djde(ind),t2s','.'); hold on; % T2* here is in ns        
        minVal=min(abs(djde(ind)));
        plot([0 1/minVal],[0 1/minVal]./noiseSlope); 
        xlabel('(dJ/d\epsilon)^{-1} (mV/MHz)'); ylabel('T_2^* (ns)');
        title(sprintf('\\epsilon_{RMS} =%g (\\muV)',epsrms*1e6));
        [plotnum,figs,ha] = nextfig(config,plotnum,fignum,figs,ha);
        loglog(djde(ind),t2s','.'); 
        fitdescr = [ fitdescr sprintf('RMS Noise: %g uV\n',epsrms*1e6)];        
        if ~isopt(config.opts,'amp'), mask = 1:size(params,1); end
        [plotnum,figs,ha] = nextfig(config,plotnum,fignum,figs,ha);
        plot(1e3*real(freqFunc(params(mask,:))),real(djde(mask)),'.-');
        xlabel('J'); ylabel('dJ/d\epsilon (MHz/mV)');
    end
    if isopt(config.opts,'echocenter')
        [plotnum,figs,ha] = nextfig(config,plotnum,fignum,figs,ha);
        plot(xvals,params(:,5),'.-');
        xlabel(config.xlabel); ylabel('Echo center (ns)');
    end
    if isopt(config.opts,'echophase')
        [plotnum,figs,ha] = nextfig(config,plotnum,fignum,figs,ha);
        plot(xvals,unwrap(atan2(params(:,3),params(:,2))),'.-');
        xlabel(config.xlabel); ylabel('Echo Phase (radians)');
    end
    if ~isopt(config.opts,'nofull') % Plot the full data set (not averaged).
        [plotnum,figs,ha] = nextfig(config,plotnum,fignum,figs,ha);
        currData=reshape(permute(data,[1 3 2]),size(data,1),size(data,2)*size(data,3));
        imagesc(currData(config.frames,:));
    end
    if isopt(config.opts,'echonoise') && isfield(config,'djde') && ~isempty(config.djde)
        %this only works with power and exponential decay. not 'both'
        djde = config.djde*1e9; %MHz/mV into Hz/V
        if isopt(config.opts,'power')
            beta = fp(4)-1;
            gval=gamma(-1-beta)*sin(pi*beta/2);
        else
            beta = 0;
            gval=pi/2;  % Limit of gval above as beta->0
        end
        pars.Sphi=2*pi/abs(2^-beta*(-2+2^beta)*(1e-6*pars.T2)^(1+beta)*gval);
        pars.Seps2=@(f) pars.Sphi*(2*pi)^(-1-beta)*(djde)^(-2)/f^beta;
        fitdescr = [ fitdescr sprintf('Noise@1Mhz: %g nV (%g nV)(beta=%g)\n',sqrt(pars.Seps(1e6))*1e9,sqrt(pars.Seps2(1e6))*1e9,beta)];
    end
    if isopt(config.opts,'color')
        figure(1122);
        colorData = squeeze(nanmean(data(config.frames,:,:)));
        imagesc(xv(:,1),xvals,colorData,'ButtonDownFcn',@btn);
        colorbar; a = gca; a.YDir = 'Normal'; 
    end
end
figs = unique(figs);
for i = 1:length(figs)
    formatFig(figs(i),'exch full',config.spsize(1),config.spsize(2));
end
if isa(figs,'matlab.ui.Figure'), figs = [figs.Number]; end
if ~isopt(config.opts,'noppt')
    prettyfile=regexprep(file,'(sm_)|(\.mat)','');
    indentdescr= regexprep(fitdescr,'^(.)','\t$1','lineanchors');
    ppt=guidata(pptplot);
    set(ppt.e_file,'String',file);
    set(ppt.e_figures,'String',['[',sprintf('%d ',figs),']']);
    set(ppt.e_title,'String',prettyfile);
    set(ppt.e_body,'String',fitdescr);
    clipboard('copy',sprintf('%s\n%s\n\n',['===' prettyfile],indentdescr));
    set(ppt.exported,'Value',0);
end
fprintf(sprintf('%s\n%s\n\n',['===' prettyfile],indentdescr));
end

function [fitpars,fitfn]=fitosc(x,y,opts,rng,style)
% def'n fitpars: 
% for all: 1: offset 4: freq
% NoDecay: 2: cos coef, 3 sin coef
% Center: 2: cos coef, 3 sin coef, 5: decay center, y(6) 1/t2
% Gen: 2: amp, 3: phase, 5: center, 6: 1/t2
% If style given, plot with given style. 

fitfn.fn = @fioscill; fitfn.args = {1};
cosNoDecay = '@(y,x) y(1) + y(2) * cos(y(4)*x) + y(3) * sin(y(4)*x)';
cosCenter =  '@(y,x) y(1) + (y(2) * cos(y(4) * x) + y(3) * sin(y(4) * x)).*exp(-(x-y(5)).^2 * y(6).^2)';
cosGen =     '@(y,x) y(1) + y(2) * cos(y(4) * x + y(3)).*exp(-(x-y(5)).^2 * y(6).^2)';
if exist('rng','var') && ~isempty(rng)
    pts = x >= rng(1) & x <= rng(2);
    x = x(pts); y = y(pts);
end
if ~isopt(opts,'fitdecay') || (isopt(opts,'afitdecay') && std(y) < 2e-2)
    fitfn.args={2}; % No decay
    fitpars=fitwrap('fine noplot',x,y,fitfn, cosGen, [1 0 1 1 0 0]);
    ig = [fitpars(1), fitpars(2)*cos(fitpars(3)), fitpars(2)*sin(fitpars(3)), fitpars(4:6)];
    fitpars=fitwrap('fine noplot',x,y,ig, cosNoDecay, [1 1 1 1 0 0]);
    fitfn=str2func(cosNoDecay);
elseif ~isopt(opts,'nocenter') % Decay and center
    fitfn.args={2};
    fitpars=fitwrap('fine noplot',x,y,fitfn,cosCenter, [1 1 1 1 0 0]);
    fitpars=fitwrap('fine noplot',x,y,fitpars, cosCenter, [1 1 1 1 1 1]);
    fitfn=str2func(cosCenter);
else  % Decay but no center
    fitfn.args={2};
    fitpars=fitwrap('fine noplot',x,y,fitfn, cosGen, [1 0 1 1 0 0]);
    fitpars = [fitpars(1), fitpars(2)*cos(fitpars(3)), -fitpars(2)*sin(fitpars(3)), fitpars(4:6)];
    fitpars=fitwrap('fine noplot',x,y,fitpars, cosCenter, [1 1 1 1 0 1]);
    fitfn=str2func(cosCenter);
end
if ~isopt(opts,'noplot')
    if isopt(opts,'interp'), x=linspace(x(1,1),x(1,end),1024); end
    if exist('style','var') && ~isempty('style')
        plot(x,fitfn(fitpars,x),style);
    else
        plot(x,fitfn(fitpars,x),'r');
    end
end
end

function [fitpars,fitfn,fitstring]=fitdecay(x,y,opts,rng,style)
% options: gauss, both, power, (none)
% gauss:

if ~isopt(opts,'fitoffset')
    mask=[1 1 0];
else
    mask=[1 1 1];
end
if exist('rng','var') && ~isempty(rng)
    pts = x > rng(1) & x < rng(2);
    x = x(pts); y = y(pts);
end
if isopt(opts,'gauss')
    fitfn=@(p,x) p(1)*exp(-(x/p(2)).^2)+p(3);
    init=[1 max(x)/3 min(y)];
    init(1)=range(y)/(fitfn(init,min(x)-init(3)));
    fmt='Decay: %.3g exp(-(t/%.3g)^2)+%.3g\n';
    fperm=[1 2 3];
elseif isopt(opts,'both')
    fitfn=@(p,x) p(1)*exp(-(x/p(2)).^2-x/p(4))+p(3);
    init=[1 max(x)/3 min(y) max(x)/10];
    init(1)=range(y)/(fitfn(init,min(x))-init(3));
    mask(4)=1;
    fperm=[1 2 3];
    fmt='Decay: %.3g exp(-t/%3g -(t/%3g)^2)+%g\n';
elseif isopt(opts,'power')
    fitfn=@(p,x) p(1)*exp(-(x/p(2)).^p(4))+p(3);
    init=[1 max(x)/3 min(y) 1];
    mask(4)=1;
    init(1)=range(y)/(fitfn(init,min(x))-init(3));
    fperm=[1 2 4 3];
    fmt='Decay: %.3g exp(-(t/%3g)^{%g})+%g\n';
else
    fitfn=@(p,x) p(1)*exp(-x/p(2))+p(3);
    init=[1 max(x)/3 min(y)];
    init(1)=range(y)/(fitfn(init,min(x))-init(3));
    fmt='Decay: %.3g exp(-t/%3g)+%g\n';
    fperm=[1 2 3];
end
if ~isopt(opts,'fitoffset'), init(3)=0; end
fitpars = fitwrap('noplot',x,y,init,fitfn,mask);
fitstring=sprintf(fmt,fitpars(fperm));
if isopt(opts,'plot')    
    hold on;
    if isopt(opts,'interp')
        x=linspace(x(1,1),x(1,end),512);
    end
    if exist('style','var') && ~isempty('style')
        plot(x,fitfn(fitpars,x),style);
    else
        plot(x,fitfn(fitpars,x));
    end
end
end

function [plotnum,figs,ha] = nextfig(config, plotnum, fignum, figs, ha)
% Returns next axis. 
nplot = prod(config.spsize);
figure(1+fignum+floor((plotnum-1)/nplot)); % Check which figure we are on.
if(mod(plotnum-1,prod(config.spsize)) == 0) % If this is the first time we've used the figure.
    clf; ha = tightSubplot(config.spsize,'title');
end
axes(ha(mod(plotnum-1,nplot)+1));
figs=unique([figs gcf]);
plotnum=plotnum+1;
end