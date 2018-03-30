function [pars,sdata,figs ]=anaExchange(file,config)
% Analyze ramsey and echo noise, performing basic fits. 
%function [figs pars sdata] = anaExchange(file,config)
% config is a config struct describing what to do.  File can be blank
% Note that the file name is cached in a persistent variable.
% Possible fields in config:
% ts: the total evolution times as a vector, ie. 0.05:0.05:2.  If ts
%   is nan(default), try to guess from group xvals.
% t1: Default ratio of t1 to tmeas
% dbz: a list of which channels have dbz data. defaults to any pulsegroup w/ dbz in the name.
% rng: points to fit. [1 10] means first 10 xvals.  [5 inf] means 5 onwared.
% grng: groups to fit.  [ .5 inf] means .5 onward.  uses group xvals
% frames: reps to fit in the scan. (5:20) means fit reps 5-20.
% channel: a list of channels to fit.
% nodbz: skip dbz fit
% offset: offset between lines in figure 100
% side: which set in tunedata to look at .
% djde: djde in MHz/mV for computing epsilon noise
% acut: cutoff amplitude
% opts: general set of options that can have the following words
%   afitdecay: fit decay only if amplitude is bigger than 5%
%   amp: Plots amplitude vs. epsilon 
%   epsrms or epsRMS: calculate rms epsilon noise.
%   echonoise: computes high frequency noise based on echo decay using the specified value of djde
%   echocenter: plots center time of echo curve vs evo time
%   echophase: plots phase of echo curve vs evo time
%   even /odd fit all even or odd groups only. 
%   fitdecay: fits for the decay envelope
%   frq: Plots freq vs epsilon
%   gauss
%   guessxval: tries to guess the xvals (times) from the group def.
%   linfit / logfit
%   lines: connect the data with lines (not fitted). 
%   mean
%   nocenter: constrains the center of the decay to be at zero time
%   nocolor: Don't plot the full data set (unaveraged). 
%   nodbz: no dbz group included in the scan
%   nofitdbz: do not fit the dbz group
%   nofit: don't perform fit of data. 
%   noscale / linescale. Changes default histogramming behavior. no scale doesn't histogram, linescale does ech line separately. 
%   offsetoff: 
%   power: fits the decay exponent
%   per: Plots period vs epsilon
%   ramsey: single evo ramsey
%   ramq
%   ramt2: plots the ramsey t2 time
%   rmoutlier
%   singlegroup: only one group in whole thing -- average them together. 

figs=[]; pars = []; plotnum=1; fitdescr=''; cache=0;
persistent lastname; persistent sdata_cache; persistent fileinfo;
if ~exist('file','var') || isempty(file)
    if ~ischar(lastname), lastname=''; end
    file=uigetfile('sm*.mat','ana_echo',lastname);
end
if isempty(file), return; end
if strcmp(lastname,file) && ~isempty(sdata_cache) && ~isempty(fileinfo) % if p
    st=dir(file);
    if st.bytes == fileinfo.bytes && st.datenum == fileinfo.datenum, cache=1; end
    fileinfo=st;
end
lastname=file;
if ~exist('config','var'), config=struct(); end
if iscell(config), config=struct(config{:}); end
if ischar(config), config=struct('opts',config); end
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
scantime=getFileTime(file); pars.scantime = scantime;
config = def(config,'opts','ramsey');   % Random boolean options
config = def(config,'rng',[]); % Range of points to fit, ie [10 inf];
config = def(config,'grng',[]); % Group range points to fit, ie [1.1 1.5];
config = def(config,'channel',1); % Range of channels to fit. fixme to be smarter
config = def(config,'t1',nan); % Default ratio of t1 to tmeas.  fixme to be smarter.
config = def(config,'frames',[]); % Default frames to fit.
config = def(config,'dxlabel','T (ns)');  % Default xlabel for intra-group series
config = def(config,'fitopts','interp');  % Interpolate xvals for fits
config = def(config,'fb',100); % Base number of figures to output.
config = def(config,'spsize',[3, 3]); % Number of subplots for 'minor' figures.
config = def(config,'dbz', find(cellfun(@(p) ~isempty(p),regexp({scan.data.pulsegroups.name},'[dD][bB][zZ]')))); % Is there a dBz group?
config = def(config,'side',[]); % Side of dot examined.
config = def(config, 'ts', nan); % xvals across groups
config = def(config,'acut',2); % Cutoff amplitude
ha=[];
if isopt(config.opts,'ramsey') % Sets up default options, xlabel for ramsey. 
    if isopt(config.opts, 'singlegroup')
        config.opts = [config.opts ' guessxval fitdecay nocenter gauss color'];
    else
        config.opts=[config.opts ' guessxval amp per frq gauss fitdecay nocenter ramq ramt2 epsRMS color nodbz'];
    end
    config = def(config,'xlabel','\epsilon (mV)');  % Default xlabel for inter-group series
else
    config = def(config,'xlabel','T (\mus)'); % Default xlabel for inter-group series
end
if isopt(config.opts,'echo') % Period, amplitude, frequency vs. T
    config.opts=[config.opts ' guessxval freq amp per color nodbz']; 
    config = def(config,'offset',0.15); 
end 
if isempty(config.frames), config.frames=1:size(data{config.channel},1); end % By default use all frames. 
if isempty(config.side) % Use DAQ channel to determine side. 
    switch scan.loops(1).getchan{1}
        case 'DAQ2'
            config.side='right';
        case 'DAQ1'
            config.side='left';
        otherwise
            error('Unable to determine side');
    end
end
if isnan(config.t1), [~, config.t1] = att1(config.side,scantime,'before',scan); end % Find t1 before rescaling data. 
notdbz=setdiff(1:length(scan.data.pulsegroups),config.dbz); % Set of groups other than dBz check (i.e.interesting data).
ngrps=length(scan.data.pulsegroups);
if ngrps==1 % Find the dt's for individual scan lines.
    ngrps=scan.data.conf.nrep; % Each rep fitted separately. 
    notdbz=(1:ngrps);
    dt = scan.data.pulsegroups(1).varpar';    
    ts = scanRng(scan,1); % Assume that sweep is along first loop of scan. 
    config.ts = ts;
else % Guess xvals
    for i=1:ngrps
        xvt=scan.data.pulsegroups(i).varpar';
        dxvt=diff(xvt,[],2) ~= 0;
        [~, ind]=max(sum(dxvt,2));
        dt(i,:)=xvt(ind,:); %#ok<*AGROW>
        if ismember(i,notdbz), xv(:,i)=xvt(:); end
    end
end
config = def(config, 'dt', dt); % dts
dt = config.dt; %hack to keep backward compatibility
if any(size(dt)==1), dt = repmat(dt,size(data{1},2),1); end
if isnan(config.ts)    
    if isopt(config.opts,'guessxval') % Guess the group xval from the params
        npars=length(scan.data.pulsegroups(1).params); 
        pulseparams = [scan.data.pulsegroups.params];        
        pulseparams = pulseparams';
        pulseparams = reshape(pulseparams,npars,length(pulseparams)/npars); 
        if size(pulseparams,2)>1
            r=find(diff(pulseparams(:,notdbz),[],2) ~= 0); %Find the param that changes
            if isnan(mode(r))
                ts=1:length(notdbz);
            else
                ts=pulseparams(mode(r),:)';
            end
        else
            ts = pulseparams'; 
        end
    else
        r = find(diff(xv,[],2) ~= 0);
        if ~isempty(r)
            ts = xv(mode(r),:)';
        else
            fprintf('No xval changes from group to group.  Try guessxval\n');
            ts = xv(1,:);
        end
    end
else
    ts=config.ts;
    if (length(ts) ~= ngrps), error('The length of TS (%d) must be the same as the number of groups (%d)\n',length(ts),length(scan.data.pulsegroups)); end
end
ts=ts(:);
if isopt(config.opts,'even'), notdbz=notdbz(2:2:end); end % only ana odd / evengroups
if isopt(config.opts,'odd'), notdbz=notdbz(1:2:end); end
if isopt(config.opts,'linescale') % scale data line by line.
    t1=ones(config.channel,1).*config.t1; % will break when ana_echo can do more than one channel at once.
    %tMeas = scan.data.pulsegroups(1).readout.readout(3);    
    data_all=anaHistScaleLine(scan,data,t1);
    yl='P(T)';
    if isfield(config,'offset')
        offset=config.offset;
    else
        offset=1/3;
    end
elseif ~isopt(config.opts,'noscale') % Scale the data by histogramming.
    %[~,ratio] =att1(config.side,scantime,'before',scan,data);    
    data_all=anaHistScale(scan,data,config.t1);
    yl='P(T)';
    if isfield(config,'offset')
        offset=config.offset;
    else
        offset=1/2;
    end
else
    data_all=data;
    yl='V_{rf} (mV)';
    offset=4e-4;
end
if isopt(config.opts, 'offsetoff'), offset = 0; end
if isfield(scan.data,'setpt') % Use set point for dbz stored in scan
    dbzFreq = scan.data.setpt(1)*2*pi/1000; 
else
    dbzFreq = 0.0315 * 2 * pi; 
end
for i=1:length(config.channel) % Fitting and plotting.
    data=data_all{config.channel(i)};
    fb=config.fb+100*(i-1);
    if size(dt,1)>20
        newFig = floor(size(dt,1)/2);
    else
        newFig = nan; 
    end
    params=[]; o=0; % initialize offset
    
    for j=1:length(notdbz) % Fit and plot all the nondbz data.
        ind=notdbz(j);
        if isopt(config.opts, 'singlegroup') % If so, we'll average all the data together.
            rdata = squeeze(nanmean(data(config.frames,:)));
        else
            rdata=squeeze(nanmean(data(config.frames,ind,:),1))'+o;
        end
        o=o+offset;
        if ~isopt(config.opts,'noplot')
            if j==1, figure(fb); figs=[figs gcf]; clf; ylabel(yl); xlabel(config.dxlabel); hold on; end %prepare figure
            if j==newFig, figure(2*fb+1); figs=[figs gcf]; clf; ylabel(yl); xlabel(config.dxlabel); hold on; end
            a = gca; 
            if isopt(config.opts,'lines') % Plot lines with data.
                plot(dt(ind,:),rdata,'.-');
            else
                plot(dt(ind,:),rdata,'.');
                colorOrd = a.ColorOrderIndex;
                if colorOrd>1, colorOrd = colorOrd-1;
                else, colorOrd = 7; end
                a.ColorOrderIndex = colorOrd;
            end
            if ~isopt(config.opts,'plot'), config.opts = [config.opts ' plot']; end
        end
        if ~isopt(config.opts,'nofit')
            params(j,:)=fitosc(dt(ind,:),rdata,config.opts,config.rng,'-');
        end
    end
    if isopt(config.opts,'nofit'), return; end
    pars.params=params; pars.ts=ts'; pars.dt=dt; pars.dt1=dt(1,:);
    if length(params) > 3
        if size(params,1) > 10
            jbar=mean(params(3:end-3,4));
        else
            jbar=mean(params(:,4));
        end
    else
        jbar=nan;
    end
    title(sprintf('|J|=%g MHz',1e3*jbar/(2*pi)));
    fitdescr = [fitdescr sprintf('|J|=%g MHz\n',1e3*jbar/(2*pi))];
    if ~isopt(config.opts,'nodbz')
        for j=1:length(config.dbz)
            dbzdata=squeeze(nanmean(data(config.frames,config.dbz(j),:),1))';
            [plotnum,figs,ha] = nextfig(config,plotnum,fb,figs,ha);
            plot(dt(1,:),dbzdata,'b.');
            xlabel('T (ns)'); ylabel(yl);            
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
    end
    if size(ts,2)>1, ts=ts(notdbz); end
    if isopt(config.opts,'amp') % Plot amplitude as a function of eps, fit T2. 
        [plotnum,figs,ha] = nextfig(config,plotnum,fb,figs,ha);
        ampfunc=@(x) 2*(sqrt(x(:, 2).^2 + x(:, 3).^2));
        amp=ampfunc(params);        
        if ~isnan(config.acut) % Only include data with reasonable amplitude. 
            cut = config.acut; 
            mask=amp < cut;
        else
            mask=~isnan(amp);
        end
        if isopt(config.opts, 'rmoutlier'), mask = mask & amp <1.2; end
        plot(ts(mask),ampfunc(params(mask,:)),'b.');
        if isopt(config.opts,'linfit')
            ind = find(ts' > config.grng(1) & ts' < config.grng(2) & mask);
            fpp = fitwrap('plfit',ts(ind)', ampfunc(params(ind,:))', [-1 0], @(p,x) p(1)*x + p(2));
            str = [str sprintf('Amp %g t + %g',fpp(1),fpp(2))];
            fp=[]; fp(1) = fpp(2); fp(2) = (fpp(2)+fpp(1)*mean(ts(ind)))/fpp(1);
        elseif isopt(config.opts,'logfit')
            ind = find(ts' > config.grng(1) & ts' < config.grng(2) & mask(:)');
            fpp = fitwrap('plfit',ts(ind)', log(ampfunc(params(ind,:)))', [-1 0], @(p,x) p(1)*x + p(2));
            str = [str sprintf('Amp exp^(%g t + %g)',fpp(1),fpp(2))];
            fp=[]; fp(1) = fpp(2); fp(2) = 1/fpp(1);
        else
            [fp,~,fstr]=fitdecay(ts(mask)',ampfunc(params(mask,:))',['plot ' config.opts],config.grng);
            str=['Amp' fstr];
            str= [str sprintf('Amp=%.3f T_2^{echo}=%.3g Q=%3.1g T=%.3g',fp(1),fp(2),fp(2)/(1e-3*2*pi/jbar),2*pi/mean(jbar))];           
        end        
        fitdescr = [ fitdescr 'Amp: ' str newline ];
        title(str); ylabel(yl); xlabel(config.xlabel);        
        if ~isempty(config.grng)
            ind = find(ts > config.grng(1) & ts < config.grng(2));
        else
            ind=1:length(ts);
        end
        pars.maxamp=max(ampfunc(params(ind,:)));
        pars.T2 = fp(2); pars.Q = fp(2)/(1e-3*2*pi/jbar); pars.Jbar = jbar;
        pars.T= 2*pi/mean(jbar); pars.afp = fp; pars.amp = fp(1);  
        pars.ampDecay = ampfunc(params(mask,:))'; 
    end   
    if isopt(config.opts,'per') % Plot period as a function of eps
        [plotnum,figs,ha] = nextfig(config,plotnum,fb,figs,ha);
        perfunc=@(params) 2*pi./params(:,4);
        if ~isopt(config.opts,'amp'), mask = 1:size(params,1); end
        mask2 = abs(perfunc(params))<1e4; 
        mask = mask & mask2; 
        plot(ts(mask),perfunc(params(mask,:)),'b.');
        title('Period'); ylabel('T (ns)'); xlabel(config.xlabel);
    end    
    if isopt(config.opts,'frq') % Plot J vs. epsilon 
        if ~isopt(config.opts,'noref')
            freqfunc=@(params) sqrt( abs(params(:,4)./(2*pi)).^2-((dbzFreq/((2*pi)))^2)) .* sign(params(:,4) - dbzFreq);
        else
            freqfunc=@(params) params(:,4)./(2*pi);
        end
        aliasFreq = 0.5; 
        aliasedData = find(diff(freqfunc(params))<0); 
        params(aliasedData+1,4)=-params(aliasedData+1,4)+aliasFreq*4*pi; 
        [plotnum,figs,ha] = nextfig(config,plotnum,fb,figs,ha);
        
        plot(ts,1e3*real(freqfunc(params)),'b.');
        ylabel('J (MHz)'); xlabel(config.xlabel);                
        if isopt(config.opts,'fitoffset')
            jFit = fitdecay(ts',1e3*freqfunc(params)','plot fitoffset',config.grng);
        else
            jFit = fitdecay(ts',1e3*freqfunc(params)','plot',config.grng);                        
        end
        str=sprintf('Decay const = %2.2f', jFit(2)); title(str);
        fitdescr = [ fitdescr 'Freq: ' str newline ];
        fitdescr = [fitdescr sprintf('J(eps) = 100 * exp(-(eps-%.4g)/%.3g) + %.3g Mhz',log(jFit(1)/100)*jFit(2),jFit(2),jFit(3))] ;
        pars.freqfunc=@(eps) 100*exp(-(eps-log(jFit(1)/100)*jFit(2))/jFit(2))+jFit(3);
        pars.jFit = jFit;        
        if isopt(config.opts,'ramsey')
            pars.eps=ts';
        elseif isopt(config.opts,'echo')
            pars.tau=ts';
        end
        pars.freq = 1e3*freqfunc(params)';
    end 
    if isopt(config.opts,'ramt2') % Plot T2 vs. epsilon (plot all of these). 
        [plotnum,figs,ha] = nextfig(config,plotnum,fb,figs,ha);
        t2s = abs(1./params(:,6)); pars.t2s = t2s';
        mask = mask & t2s < 1000; 
        plot(ts(mask),t2s(mask),'.-');
        xlabel(config.xlabel); ylabel('T_2^* (ns)');        
        [plotnum,figs,ha] = nextfig(config,plotnum,fb,figs,ha);
        if ~isopt(config.opts,'amp'), mask = 1:size(params,1); end        
        plot(1e3*real(freqfunc(params(mask,:))),t2s(mask),'.-');
        xlabel('J (MHz)'); ylabel('T_2^* (ns)');
    end    
    if isopt(config.opts,'epsRMS') % Plot djdEps vs T2*, fit slope for low freq. noise
        [plotnum,figs,ha] = nextfig(config,plotnum,fb,figs,ha);        
        djde = (1/jFit(2))*(1e3*freqfunc(params)-jFit(3)); %Offset does not contribute to slope
        noiseslope = (1./djde(ind))'/abs(1./(params(ind,6)))'; % this is matlab shorthand for least squares slope
        plot([0 1/(djde(ind(end)))],[0 1/(djde(ind(end)))]./noiseslope);
        erms=sqrt(2)*noiseslope/(2*pi); pars.erms = erms;
        pars.djde = djde';
        % dJ/dE is in MHz/mV=GHz/V, t2* is in ns, so erms in v.
        plot(1./djde(ind),abs(1./params(ind,6)),'.'); hold on; % T2* here is in ns        
        xlabel('(dJ/deps)^{-1} (mV/MHz)'); ylabel ('T2^* (ns)');
        title(sprintf('RMS noise =%g (\\muV)',erms*1e6));
        fitdescr = [ fitdescr sprintf('RMS Noise: %g uV\n',erms*1e6)];
        [plotnum,figs,ha] = nextfig(config,plotnum,fb,figs,ha);
        if ~isopt(config.opts,'amp'), mask = 1:size(params,1); end
        plot(1e3*real(freqfunc(params(mask,:))),real(djde(mask)),'.-');
        xlabel('J'); ylabel('dJ/d\epsilon (MHz/mV)');
    end
    if isopt(config.opts,'echocenter')
        [plotnum,figs,ha] = nextfig(config,plotnum,fb,figs,ha);
        plot(ts,params(:,5));
        xlabel(config.xlabel); ylabel('Echo center (ns)');
    end
    if isopt(config.opts,'echophase')
        [plotnum,figs,ha] = nextfig(config,plotnum,fb,figs,ha);
        plot(ts,unwrap(atan2(params(:,3),params(:,2))));
        xlabel(config.xlabel); ylabel('Echo Phase (radians)');
    end
    if ~isopt(config.opts,'nocolor') % Plot the full date set (not averaged). 
        [plotnum,figs,ha] = nextfig(config,plotnum,fb,figs,ha);
        rdata=reshape(permute(data,[1 3 2]),size(data,1),size(data,2)*size(data,3));
        imagesc(rdata(config.frames,:));
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
        imagesc(xv(:,1),ts,colorData,'ButtonDownFcn',@btn); 
        colorbar; 
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
    ppt=guidata(pptplot3);
    set(ppt.e_file,'String',file);
    set(ppt.e_figures,'String',['[',sprintf('%d ',figs),']']);
    set(ppt.e_title,'String',prettyfile);
    set(ppt.e_body,'String',fitdescr);
    clipboard('copy',sprintf('%s\n%s\n\n',['===' prettyfile],indentdescr));
    set(ppt.exported,'Value',0);
end
fprintf(sprintf('%s\n%s\n\n',['===' prettyfile],indentdescr));
end

function [fp,ff]=fitosc(x,y,opts,rng,style)
% 1: offset 4: freq 
% NoDecay: 2: cos coef, ,3 sin coef
% Center: 2: cos coef, ,3 sin coef, 5: decay center, y(6) 1/t2
% Gen: 2: amp, 3: phase, 5: center, 6: 1/t2
fig=gcf;
fifn.fn = @fioscill; fifn.args = {1};
cosNoDecay = '@(y, x)y(1)+y(2)*cos(y(4)*x) + y(3) * sin(y(4)*x)';
cosCenter = '@(y, x) y(1) + (y(2) * cos(y(4) * x) + y(3) * sin(y(4) * x)).*exp(-(x-y(5)).^2 * y(6).^2)';
cosGen =    '@(y, x) y(1) + y(2) * cos(y(4) * x + y(3)).*exp(-(x-y(5)).^2 * y(6).^2)';
if exist('rng','var') && ~isempty(rng)
    pts=x>rng(1) & x <rng(2);
    x=x(pts); y=y(pts);
end
if ~isopt(opts,'fitdecay') || (isopt(opts,'afitdecay') && std(y) < 2e-2)
    fifn.args={2}; % No decay
    fp=fitwrap('fine',x,y,fifn, cosGen, [1 0 1 1 0 0]);
    ig = [fp(1), fp(2)*cos(fp(3)), fp(2)*sin(fp(3)), fp(4:6)];
    fp=fitwrap('fine',x,y,ig, cosNoDecay, [1 1 1 1 0 0]);
    ff=str2func(cosNoDecay);
elseif ~isopt(opts,'nocenter') % Decay and center
    fifn.args={2}; 
    fp=fitwrap('fine',x,y,fifn,cosCenter, [1 1 1 1 0 0]);
    fp=fitwrap('fine',x,y,fp, cosCenter, [1 1 1 1 1 1]);
    ff=str2func(cosCenter);
else  % Decay but no center
    fifn.args={2};
    fp=fitwrap('fine',x,y,fifn, cosGen, [1 0 1 1 0 0]);
    fp = [fp(1), fp(2)*cos(fp(3)), -fp(2)*sin(fp(3)), fp(4:6)];
    fp=fitwrap('fine',x,y,fp, cosCenter, [1 1 1 1 0 1]);
    ff=str2func(cosCenter);
end
if isopt(opts,'plot')
    figure(fig); hold on;
    if isopt(opts,'interp'), x=linspace(x(1,1),x(1,end),512); end
    if exist('style','var') && ~isempty('style')
        plot(x,ff(fp,x),style);
    else
        plot(x,ff(fp,x),'r-');
    end
end
end
function [fitpars,fitfn,fitstring]=fitdecay(x,y,opts,rng,style)
% options: gauss, both, power, (none) 
% gauss: 
fig=gcf;
if ~isopt(opts,'fitoffset')
    mask=[1 1 0];
else
    mask=[1 1 1];
end
if exist('rng','var') && ~isempty(rng)
    pts=x>rng(1) & x <rng(2);
    x=x(pts); y=y(pts);
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
fitpars = fitwrap('plinit plfit',x,y,init,fitfn,mask);
fitstring=sprintf(fmt,fitpars(fperm));
if isopt(opts,'plot')
    figure(fig);
    hold on;
    if isopt(opts,'interp')
        x=linspace(x(1,1),x(1,end),512);
    end
    if exist('style','var') && ~isempty('style')
        plot(x,fitfn(fitpars,x),style);
    else
        plot(x,fitfn(fitpars,x),'r-');
    end
end
end
function [plotnum,figs,ha] = nextfig(config, plotnum, fignum, figs, ha)   
   nplot = prod(config.spsize); 
   figure(1+fignum+floor((plotnum-1)/nplot)); % Check which figure we are on. 
   if(mod(plotnum-1,prod(config.spsize)) == 0) % If this is the first time we've used the figure. 
       clf; ha = tightSubplot(config.spsize); 
   end
   axes(ha(mod(plotnum-1,nplot)+1)); 
   figs=unique([figs gcf]); 
   plotnum=plotnum+1;
end