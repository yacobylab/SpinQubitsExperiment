function [out,figs]=mfitEchoFish(file,config)
% Constrained fits of echo data to find values for charge noise.
% function [out,figs]=mfitEchoFish(files,opts)
% opts.grps ([1 inf] default) Which pulsegroups to fit.
% opts.opts ('residuals colorplot plotdbz' default )
%   colorplot 
%   residuals
%   fishless
%   noerrorbars
% opts.tvrange ([-inf,inf] default) Which times to use.
% opts.mfitopts ('plinit plfit plotiterr err' default)
% opts.frames ; Which reps of data to use (e.g. if gradient runs away,
% don't use that)
%  parameters with _i are allowed to vary from curve to curve.
% Parameters kept the same: t2, alpha, a (amplitude), t2*, fish_a (amplitude)
% allowed to vary: freq (rad/s), phase (rad), offset, time offset
% The "fish" is the initial decay in t2* (for short tau), out of phase with
% main oscillation
% model: echo signal: y0_i+a*cos(w_i * x + phi_i) exp(-((x-x0_i)/t2*)^2) exp(-(tau/t2)^alpha)
%        + fish_a * cos(w_i * x + fish_phi_i) exp(-((tau/2 + x)/t2*)^2)
% parameters returned will have this form
%   1   2   3   4     5  6   7      8         9     10
% [t2,alpha,a,fish_a,t2*,w_i,phi_i,fish_phi_i,y0_i,x0_i]
% in out.params, returns with 1st row: J, 2nd: phase, 3rd: fish phase,
% 4th: y offset 5: t offset (for decay)
if ~exist('config','var'), config=struct(); end
if iscell(config), config=struct(config{:}); end

config=def(config,'grps',[1 Inf]);
config=def(config,'opts','residuals colorplot nocenter plotfit'); % nocenter also allowed
config=def(config,'taurange',[-Inf Inf]);
config=def(config,'mfitopts','plinit plfit plotiter err');
config=def(config,'frames',[]);
config=def(config,'fignum',1024); fignum = config.fignum;
config=def(config,'rng',[-Inf,Inf]); rng = config.rng;
figs=[];
nParPerm = 5; nParVar = 5; tCenterVar = 10; phaseVar = 7; freqVar = 6; fishAmp = 4; fishPhase = 8;
% Kind of annoying but I think all of the vars > 6 are actually smaller by
% 1 than in fitfn because tau doesn't count.
if ~exist('file','var') || isempty(file)
    if ~isopt(config.opts,'nofilt')
        file=uigetfile('sm*Ramsey*.mat');
    else
        file=uigetfile('sm*.mat');
    end    
end
out.file = file; 
%% Load the data.
% Get scaled data and pulsegroup info data
s=procPlsData(file,struct('grps',config.grps,'opts','noplot noppt','frames',config.frames));
if ~isfield(s,'data'), return; end
if ~isopt(config.opts,'quiet')
    fprintf('Processing file %s \n',file(1:end-4))
end
% s.tv = tau value
goodTau=find(s.tv > config.taurange(1) & s.tv < config.taurange(2));
s.tv=s.tv(goodTau); s.grps=s.grps(goodTau);

dataAll=squeeze(nanmean(s.data{1}(:,s.grps,:),1));
if any(any(isnan(s.data{1}))), fprintf('NAN data: variances will be wrong\n'); end
% denominator here gets variance of mean, only if no NANs in data
nFrames = length(find(all(~isnan(s.data{1}),[2,3]))); 
if nFrames <= 1, return; end 
xs=squeeze(nanstd(s.data{1}(:,s.grps,:),1)).^2/nFrames;

if isopt(config.opts,'colorplot') && isopt(config.opts,'residuals')
    figure(fignum); clf;
    ha = tightSubplot([1,3],'nolabely'); 
    axInd = 2;
elseif isopt(config.opts,'colorplot') 
    figure(fignum); clf;
    ha = gca; 
elseif isopt(config.opts,'residuals')
    ha = tightSubplot([1,2]); 
    axInd = 1; 
end

if isopt(config.opts,'colorplot') % Make a 2D color plot of the data and variances    
    figs=[figs fignum]; fignum=fignum+1;
    [mx,my]=meshgrid(s.xv{config.grps(1)},s.tv);
    %subplot(1,2,1);
    axes(ha(1)); 
    pcolor(mx,my,dataAll); shading flat; colorbar;
    xlabel('Time (ns)'); ylabel('\tau (\mus)');
    title('Singlet Probability');    
    %subplot(1,2,2);
    %pcolor(mx,my,xs); shading flat; colorbar;
    %xlabel('Time (ns)'); ylabel('\tau (ns)');
    %title('Variance'); 
end
%% Stuff the data structure and initial guess
if isopt(config.opts,'fishless')
% Guess   T2 = max tau, alpha, amp = std of first row of data, fish amp, T2* = 30 ns
    initial=[max(s.tv), 1.3, 2*nanstd(dataAll(1,:)), 0, 10];
else
    initial=[max(s.tv), 1.3, 5*nanstd(dataAll(1,:)), 2*nanstd(dataAll(1,:)), 10];
end
parInit=[];
for j=1:length(s.grps) % Collect data to go to mfitwrap
    % Check this: have we not already filtered for groups?
    goodGrps=find(s.xv{s.grps(j)} > rng(1) & s.xv{s.grps(j)} < rng(2));
    data(j).x=s.xv{s.grps(j)}(goodGrps);
    data(j).y=dataAll(j,goodGrps); %#ok<*AGROW>
    data(j).vary=xs(j,goodGrps);
    
    offsetInd = (j-1)*nParVar; currInd = offsetInd + nParPerm;
    %               Offset,amplitude, freq, phase,             t offset,t2*           tau/t2   alpha  fish_amp
    fitFn = '@(p,x) p(%d)+p(3)*cos(p(%d)*x+p(%d)) .* exp(-abs((x-p(%d))/p(5)).^2 - abs(%f/p(1))^p(2)) + p(4) * cos(p(%d)*x+p(%d)) .* exp(-((1e3*%f/2 + x)/p(5)).^2)'; 
    % Write in the fitfn with the current indices. 
    model(j).fn=str2func(sprintf(fitFn, currInd+4,currInd+1,currInd+2,currInd+5,s.tv(j),currInd+1,currInd+3,s.tv(j)));
    parInitOld = parInit;
    parInit=fioscill(data(j).x,data(j).y,2); % initial guess for freq, phase.
    if (nanstd(data(j).y) < 2.8*nanmean(sqrt(data(j).vary))) && j > 1 % propagate forward frequency guess when signal is very small
        %parInit(4)=initial(end-4);
        parInit(2:4) = parInitOld(2:4); 
        %fprintf('Propagating omega forward on group %d due to small signal \n',j);
    end
    %see above,       freq        phase,           amp offset,   time offset
    initial=[initial, parInit(4), parInit(3), 0.1, mean(data(j).y), 2];
end

for i=freqVar:nParVar:length(initial) % Median filter the initial frequency guess.
    initial(i)=median([initial(max(freqVar,i-nParVar)), initial(i), initial(min(end,i+nParVar))]);
end
out.s=s; out.mdata=data; out.mmod=model; out.initial=initial; 
out.p=initial;

% Start wby performing mfitwrap of each data set individually. 
% Allow only echo center and fish phase, phase to be fit. 
t2s=[];
for j = 1:length(s.grps)
    mask=initial*0;
    mask(phaseVar+nParVar*(j-1))=1; %phase
    %mask(5) = 1; % Try fitting t2*     
    if isopt(config.opts,'nocenter')
        mask(tCenterVar+nParVar*(j-1))=1; % time center
    end
    if isopt(config.opts,'debug') 
        config.mfitopts = [config.mfitopts,' plinit plfit optimplot'];
    end    
    if ~isopt(config.opts,'fishless')
        mask(8+nParVar*(j-1))=1; % fish phase
        %m1 = mask; m1(4) = 1; 
         if j == 1
             d1.x = data(1).x(1:floor(end/3)); d1.y = data(1).y(1:floor(end/3)); 
             d1.vary = data(1).vary(1:floor(end/3)); 
             mod1.fn = @(p,x) p(1) + p(2).*cos(p(4).*x+p(3)).*exp(-((1e3*s.tv(1)/2 + x)/p(5)).^2);
             beta0 = [0.55, 0.15, out.p(freqVar),3*pi/2,10];
             %pars=mfitwrap(d1,mod1,beta0,config.mfitopts,ones(1,length(beta0))); %out.p is the initial guess.                                
             pars = fitwrap('plinit plfit',d1.x,d1.y,beta0,mod1.fn); 
         end
    end
    out.p=mfitwrap(data(j),model(j),out.p,config.mfitopts,mask); %out.p is the initial guess.
    %t2s = [t2s out.p(5)];
end
%out.p(5) = nanmean(t2s); 
mask=ones(1,length(initial));
mask(2)=0; % alpha
% If all tau >> t2*, won't need to fit initial decay
if isopt(config.opts,'fishless')
    mask(fishPhase:nParVar:end)=0; % fish phase
    mask(fishAmp)=0; % fish amp
end
% Don't fit the center point of the echo (not allowing for pulse errors)
if ~isopt(config.opts,'nocenter') 
    mask(tCenterVar:nParVar:end)=0;
end
% First perform a fit without fitting alpha 
[out.p, out.chisq,out.cov]=mfitwrap(data,model,out.p,config.mfitopts,mask);
mask(2)=1;
[out.p, out.chisq, out.cov,out.exitflag,out.output]=mfitwrap(data,model,out.p,config.mfitopts,mask);
%% Put data in out struct 
out.t2 = out.p(1); out.alpha=out.p(2); out.amp = out.p(3); out.fishAmp = out.p(4); 
out.t2s = out.p(5); out.tMax = max(s.tv);
out.params = reshape(out.p(nParPerm+1:end),nParVar,(length(out.p)-nParPerm)/nParVar);
stdData = sqrt(diag(out.cov)); 
out.t2Err = stdData(1); out.alphaErr = stdData(2); out.ampErr = stdData(3); 
out.fishAmpErr=stdData(4); out.t2sErr = stdData(5); 
stdParams = reshape(stdData(nParPerm+1:end),nParVar,(length(stdData)-nParPerm)/nParVar);
out.jErr = stdParams(1,:)*1e3/2/pi; 
out.phaseDataErr = stdParams(2,:); out.ampOffsetErr = stdParams(4,:); 
out.ampOffset = out.params(4,:); out.tOffset = out.params(5,:); 
out.tOffsetErr = stdParams(5,:); 
freqData = out.params(1,:); 
out.J = freqData*1e3/2/pi; % In MHz. 
out.phaseData=out.params(2,:); out.fishPhase = out.params(3,:); 
out.phaseDataErr = stdParams(5,:); 
%% Make plots of fitted data
if isopt(config.opts,'residuals') 
    offset=0.25; % Try making this smart     
    resHt=.18; 
    for j=1:length(out.mdata)
        dataAll=out.mdata(j).x;
        y=out.mdata(j).y;
        % Plot data                
        if isopt(config.opts,'noerrorbars')
            plot(ha(axInd),dataAll,y+offset*(j-1),'.');
        else
            errorbar(ha(axInd),out.mdata(j).x,out.mdata(j).y+offset*(j-1), sqrt(out.mdata(j).vary),'CapSize',2);
        end        
        hold(ha(axInd),'on');    
        colorOrd = ha(axInd).ColorOrderIndex;
        if colorOrd>1 
            colorOrd = colorOrd-1;
        else
            colorOrd = 7;
        end
        ha(axInd).ColorOrderIndex = colorOrd;
        xVec=linspace(min(out.mdata(j).x),max(out.mdata(j).x),512);
        yVec=out.mmod(j).fn(out.p,xVec);
        plot(ha(axInd),xVec,yVec+offset*(j-1));        
        hold(ha(axInd+1),'on'); 
        % Plot residual 
        yVec=out.mmod(j).fn(out.p,dataAll);        
        plot(ha(axInd+1),dataAll,offset*resHt*(j-1)+dataAll*0,'k'); % zero error line        
        errorbar(ha(axInd+1),dataAll,y-yVec+offset*resHt*(j-1), sqrt(out.mdata(j).vary),'CapSize',1);        
    end    
    out.phaseData = unwrap(out.phaseData,pi/2);     
    out.fishPhase = unwrap(out.phaseData,pi/2);     
    if isopt(config.opts,'plotfit')
        figs = [figs fignum];
        equivTime = out.phaseData./freqData;
        ampFn = out.ampOffset + out.amp.*exp(-s.tv/out.t2);
        %ampFn = @(p,x) p(1) + p(2).*exp(-x/out.t2);
        % if equivTime > max t, probably fit failed
        equivTime(abs(equivTime)>max(data(1).x)) = nan;         
        plot(ha(axInd),-equivTime,(0:length(out.phaseData)-1)*offset+ampFn,'k.-');
        title(ha(axInd),'Fits + Phase')
        title(ha(axInd+1),sprintf('Errors, \\chi^2=%3f',out.chisq));
        fitAxis(ha(axInd));
        fitAxis(ha(axInd+1));
        xlabel(ha(axInd),'dt (ns)')
        xlabel(ha(axInd+1),'dt (ns)')
        
        figure(fignum); clf;
        ha = tightSubplot([2,2]);
        
        errorbar(ha(1),out.J,out.jErr,'.-'); ylabel(ha(1),'J (MHz)');
        jMean=mlePar(out.J,out.jErr); 
        title(ha(1),sprintf('J = %3.3f MHz',jMean)); 
        axes(ha(2));
        %plot(unwrap(phaseData,pi/4)/pi); % FIXME: Write a version that works.
        errorbar(out.phaseData/pi,out.phaseDataErr/pi,'.-'); hold on;
        % plot(out.fishPhase/pi);
        ylabel('Phase (rad/\pi)');
        axes(ha(3));
        errorbar(out.ampOffset,out.ampOffsetErr,'.-');
        ylabel('Amplitude offset');
        axes(ha(4));
        errorbar(out.tOffset,out.tOffsetErr,'.-');
        ylabel('t_{offset} (ns)');
        formatFig(fignum,'exch full',2,2);
    end
end
end