function [out,figs]=mfitEchoFish(file,config)
% Constrained fits of echo data to find values for charge noise.
% function [out,figs]=mfitEchoFish(files,opts)
% opts.grps ([1 inf] default) Which pulsegroups to fit.
% opts.opts ('residuals colorplot plotdbz' default )
%   colorplot 
%   residuals
%   fishless
%   noerrorbars
% opts.taurange ([-inf,inf] default) Which times to use.
% opts.mfitopts ('plinit plfit plotiterr err' default)
% opts.frames ; Which reps of data to use (e.g. if gradient runs away,
% don't use that)
%  parameters with _i are allowed to vary from curve to curve.
% Parameers kept the same: t2, alpha, a (amplitude), t2*, fish_a (amplitude)
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
config=def(config,'grps',[1 Inf]);
config=def(config,'opts','residuals colorplot nocenter'); % nocenter also allowed
config=def(config,'taurange',[-Inf Inf]);
config=def(config,'mfitopts','plinit plfit plotiter err');
config=def(config,'frames',[]);
config=def(config,'fignum',1024); fignum = config.fignum;
config=def(config,'xrng',[-Inf,Inf]); xrng = config.xrng;
figs=[];
nParPerm = 5; nParVar = 5; tCenterVar = 10; phaseVar = 7; freqVar = 6; fishAmp = 4; fishPhase = 8;
% Kind of annoying but I think all of the vars > 6 are actually smaller by
% 1 than in fitfn because tau doesn't count.
if ~exist('file','var') || isempty(file), file = getFiles; end
%% Load the data.
% Get scaled data and pulsegroup info data
s=procPlsData(file,struct('grps',config.grps,'opts','noplot noppt','frames',config.frames));
if ~isopt(config.opts,'quiet')
    fprintf('Processing file %s \n',file(1:end-4))
end
% s.tv = tau value
goodTau=find(s.tv > config.taurange(1) & s.tv < config.taurange(2));
s.tv=s.tv(goodTau); s.grps=s.grps(goodTau);

x=squeeze(nanmean(s.data{1}(:,s.grps,:),1));
if any(any(isnan(s.data{1}))), fprintf('NAN data: variances will be wrong\n'); end
% denominator here gets variance of mean, only if no NANs in data
xs=squeeze(nanstd(s.data{1}(:,s.grps,:),1)).^2/(size(s.data{1},1));
figure(fignum); clf;

if isopt(config.opts,'colorplot') % Make a 2D color plot of the data and variances
    ha = tightSubplot([1,3]); 
    if ~isopt(config.opts,'onefig')
        figs=[figs fignum]; fignum=fignum+1;
    end    
    [mx,my]=meshgrid(s.xv{config.grps(1)},s.tv);
    %subplot(1,2,1);
    axes(ha(1)); 
    pcolor(mx,my,x); shading flat; colorbar;
    xlabel('Time (ns)'); ylabel('\tau (\mu s)');
    title('Singlet Probability');
    axInd = 2; 
    %subplot(1,2,2);
    %pcolor(mx,my,xs); shading flat; colorbar;
    %xlabel('Time (ns)'); ylabel('\tau (ns)');
    %title('Variance');
else
   ha = tightSubplot([1,2]); axInd = 1;  
end
%% Stuff the data structure and initial guess
% Guess T2 = max time, alpha = 1.7, amp = standard deviation of data , T2* = 30 ns
% FIXME: could make T2* smarter.
if isopt(config.opts,'fishless')
    initial=[max(s.tv), 1.3, 2*nanstd(x(1,:)), 0, 10];
else
    initial=[max(s.tv), 1.3, 5*nanstd(x(1,:)), 2*nanstd(x(1,:)), 10];
end
for j=1:length(s.grps) % Collect data to go to mfitwrap
    % Check this: have we not already filtered for groups?
    goodGrps=find(s.xv{s.grps(j)} > xrng(1) & s.xv{s.grps(j)} < xrng(2));
    data(j).x=s.xv{s.grps(j)}(goodGrps);
    data(j).y=x(j,goodGrps); %#ok<*AGROW>
    data(j).vary=xs(j,goodGrps);
    
    offsetInd = (j-1)*nParVar; currInd = offsetInd + nParPerm;
    %               Offset,amplitude, freq, phase,             t offset,t2*           tau/t2   alpha  fish_amp
    fitFn = '@(p,x) p(%d)+p(3)*cos(p(%d)*x+p(%d)) .* exp(-abs((x-p(%d))/p(5)).^2 - abs(%f/p(1))^p(2)) + p(4) * cos(p(%d)*x+p(%d)) .* exp(-((1e3*%f/2 + x)/p(5)).^2)'; 
    % Write in the fitfn with the current indices. 
    model(j).fn=str2func(sprintf(fitFn, currInd+4,currInd+1,currInd+2,currInd+5,s.tv(j),currInd+1,currInd+3,s.tv(j)));
    parInit=fioscill(data(j).x,data(j).y,1); % initial guess for freq, phase.
    if (nanstd(data(j).y) < nanmean(sqrt(data(j).vary))) && j > 1 % propagate forward frequency guess when signal is very small
        parInit(4)=initial(end-4);
        fprintf('Propagating omega forward on group %d due to small signal \n',j);
    end
    %see above,      freq      phase,                    offset,   time offset
    initial=[initial, parInit(4), atan2(parInit(2),parInit(3)), 0.1, mean(data(j).y), 2];
end

for i=freqVar:nParVar:length(initial) % Median filter the initial frequency guess.
    initial(i)=median([initial(max(freqVar,i-nParVar)), initial(i), initial(min(end,i+nParVar))]);
end
out.s=s; out.mdata=data; out.mmod=model; out.initial=initial; out.p=initial;

% Start with an initial mfitwrap of each data set individually. 
for j = 1:length(s.grps)
    mask=initial*0;
    mask(phaseVar+nParVar*(j-1))=1; %phase
    if ~isopt(config.opts,'fishless')
        mask(8+nParVar*(j-1))=1; % fish phase
    end
    if isopt(config.opts,'nocenter')
        mask(tCenterVar:nParVar:end)=1; % time center
    end
    out.p=mfitwrap(data(j),model(j),out.p,config.mfitopts,mask); %out.p is the initial guess.
end

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

% Put data in out struct 
out.t2 = out.p(1); out.alpha=out.p(2); out.amp = out.p(3); out.fishAmp = out.p(4); 
out.t2s = out.p(5); out.tMax = max(s.tv);
out.params = reshape(out.p(nParPerm+1:end),nParVar,(length(out.p)-nParPerm)/nParVar);

%% Make plots of fitted data

if isopt(config.opts,'residuals') 
    offset=0.25; % Try making this smart     
    resHt=.18; 
    for j=1:length(out.mdata)
        x=out.mdata(j).x;
        y=out.mdata(j).y;
        % Plot data                
        if isopt(config.opts,'noerrorbars')
            plot(ha(axInd),x,y+offset*(j-1),'.');
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
        yVec=out.mmod(j).fn(out.p,x);        
        plot(ha(axInd+1),x,offset*resHt*(j-1)+x*0,'k'); % zero error line        
        errorbar(ha(axInd+1),x,y-yVec+offset*resHt*(j-1), sqrt(out.mdata(j).vary),'CapSize',1);        
    end
    phaseData=out.p(phaseVar:nParVar:end); freqData = out.p(freqVar:nParVar:end);
    phaseData = unwrap(phaseData,pi/2); 
    %phaseData = mod(phaseData,pi);
    equivTime = phaseData./freqData;
    plot(ha(axInd),equivTime,(1:length(phaseData))*offset+.15,'k.-');    
    title(ha(axInd),'Fits + Phase')
    title(ha(axInd+1),sprintf('Errors, \\chi^2=%3f',out.chisq));
    fitAxis(ha(axInd));
    fitAxis(ha(axInd+1)); 
    xlabel(ha(axInd),'dt (ns)')
    xlabel(ha(axInd+1),'dt (ns)')
    
    figure(fignum); clf;
    ha = tightSubplot([2,2]);
    axes(ha(1));
    plot(1e3*out.params(1,:)/2/pi); ylabel('J (MHz)');
    axes(ha(2));
    %plot(unwrap(phaseData,pi/4)/pi); % FIXME: Write a version that works. 
    plot(phaseData/pi); 
    ylabel('Phase (rad/\pi)');
    axes(ha(3));
    plot(out.params(4,:)); ylabel('Amplitude offset');
    axes(ha(4));
    plot(out.params(5,:)); ylabel('t_{offset} (ns)');
    formatFig(fignum,'exch full',2,2);
end
end