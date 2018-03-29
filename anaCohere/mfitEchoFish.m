function [out,figs]=mfitEchoFish2(files,config)
% function [out,figs]=mfitEchoFish(files,opts)
% opts.grps ([1 inf] default)
% opts.opts ('residuals colorplot plotdbz' default )
% opts.taurange ([-inf,inf] default)
% opts.mfitopts ('plinit plfit plotiterr err' default)
% opts.frames ; usual
% Simultaneous constrained fit on echo data.
%  parameters with _i are allowed to vary from curve to curve.
% kept the same: t2, alpha, a, t2*, fish_a
% allowed to vary: freq (rad/s), phase (rad), offset, time offset
% The "fish" is another out of phase oscillation, decaying at t2* (for
% short tau)
% model: echo signal: y0_i+a*cos(w_i * x + phi_i) exp(-((x-x0_i)/t2*)^2) exp(-(tau/t2)^alpha)
%         + fish_a * cos(w_i * x + fish_phi_i) exp(-((tau/2 + x)/t2*)^2)
% parameters returned will have this form
%   1   2   3   4     5  6   7      8   9      10   11
% [t2,alpha,a,fish_a,t2*,w_i,phi_i,fish_phi_i,y0_i,x0_i]
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
if ~exist('files','var') || isempty(files), files = get_files; end
%% Load the data.
s=procPlsData(files,struct('grps',config.grps,'opts','noplot noppt','frames',config.frames));
good=find(s.tv > config.taurange(1) & s.tv < config.taurange(2));
s.tv=s.tv(good); s.grps=s.grps(good);
x=squeeze(nanmean(s.data{1}(:,s.grps,:),1));
if any(any(isnan(s.data{1}))), fprintf('NAN data: variances will be wrong\n'); end
xs=squeeze(nanstd(s.data{1}(:,s.grps,:),1)).^2/(size(s.data{1},1));% denominator here gets variance of mean, only if no NANs in data
if isopt(config.opts,'colorplot') % Make a 2D color plot of the data and variances
    figure(1024); clf;
    %figs=[figs fignum]; fignum=fignum+1;
    [mx,my]=meshgrid(s.xv{config.grps(1)},s.tv);
    subplot(1,2,1);
    pcolor(mx,my,x); shading flat; colorbar;
    subplot(1,2,2);
    pcolor(mx,my,xs); shading flat; colorbar;
end
%% Stuff the data structure and initial guess
if isopt(config.opts,'fishless') % Guess T2 = max time, alpha = 1.7, amp = std. , T2* = 30 ns    
    initial=[max(s.tv) 1.3 2*nanstd(x(1,:)) 0 10];
else
    initial=[max(s.tv) 1.3 5*nanstd(x(1,:)) 2*nanstd(x(1,:)) 10];
end 
for j=1:length(s.grps) % Collect data to go to mfitwrap2    
    good=find(s.xv{s.grps(j)} > xrng(1) & s.xv{s.grps(j)} < xrng(2));
    data(j).x=s.xv{s.grps(j)}(good);    
    data(j).y=x(j,good); %#ok<*AGROW> 
    data(j).vary=xs(j,good);
    
    off = (j-1)*nParVar; bs = off + nParPerm; 
    model(j).fn=str2func(sprintf('@(p,x) p(%d)+p(3)*cos(p(%d)*x+p(%d)) .* exp(-abs((x-p(%d))/p(5)).^2 - abs(%f/p(1))^p(2)) + p(4) * cos(p(%d)*x+p(%d)) .* exp(-((1e3*%f/2 + x)/p(5)).^2)',bs+4,bs+1,bs+2,bs+5,s.tv(j),bs+1,bs+3,s.tv(j)));    
    pInit=fioscill(data(j).x,data(j).y,1); % initial guess for freq, phase. 
    if (nanstd(data(j).y) < nanmean(sqrt(data(j).vary))) && j > 1 % propagate forward frequency guess when signal is very small
        pInit(4)=initial(end-4);
        fprintf('Propagating omega forward on group %d\n',j);
    end
    %see above,      freq      phase,                    offset,   time offset
    initial=[initial, pInit(4), atan2(pInit(2),pInit(3)), 0.1, mean(data(j).y), 2];
end
for i=freqVar:nParVar:length(initial) % Median filter the initial frequency guess.
    initial(i)=median([initial(max(freqVar,i-nParVar)), initial(i), initial(min(end,i+nParVar))]);
end
out.s=s; out.mdata=data; out.mmod=model; out.initial=initial; out.p=initial;
for j=1:length(s.grps)
    mask=initial*0;
    mask(phaseVar+nParVar*(j-1))=1; %phase 
    if ~isopt(config.opts,'fishless')
        mask(8+nParVar*(j-1))=1; % fish phase
    end
    if isopt(config.opts,'nocenter')
        mask(tCenterVar:nParVar:end)=1; % time center
    end
    out.p=mfitwrapFast2(data(j),model(j),out.p,config.mfitopts,mask); %out.p is the initial guess. 
end

mask=ones(1,length(initial));
mask(2)=0; % alpha
if isopt(config.opts,'fishless')
    mask(fishPhase:nParVar:end)=0; % fish phase 
    mask(fishAmp)=0; % fish amp
end
if ~isopt(config.opts,'nocenter')
    mask(tCenterVar:nParVar:end)=0;
end
[out.p, out.chisq,out.cov]=mfitwrapFast2(data,model,out.p,config.mfitopts,mask);
mask(2)=1; 
[out.p, out.chisq, out.cov,out.exitflag,out.output]=mfitwrapFast2(data,model,out.p,config.mfitopts,mask);

out.t2 = out.p(1); out.alpha=out.p(2); out.amp = out.p(3); out.t2s = out.p(5); out.fishAmp = out.p(4); out.tMax = max(s.tv); 
out.params = reshape(out.p(nParPerm+1:end),nParVar,(length(out.p)-nParPerm)/nParVar); 
if isopt(config.opts,'residuals') % Make a pretty plot of fit qualities    
    offset=0.3;
    figure(fignum); clf; hold on;
    figs=[figs];
    for j=1:length(out.mdata)
        x=out.mdata(j).x;
        y=out.mdata(j).y;
        subplot(1,2,1); hold on; a=gca;
        if isopt(config.opts,'noerrorbars')
            plot(x,y+offset*(j-1),'.');
        else
            errorbar(out.mdata(j).x,out.mdata(j).y+offset*(j-1), sqrt(out.mdata(j).vary),'CapSize',2);
        end
        colorOrd = a.ColorOrderIndex; 
        if colorOrd>1, colorOrd = colorOrd-1; 
        else, colorOrd = 7; end
        a.ColorOrderIndex = colorOrd; 
        xVec=linspace(min(out.mdata(j).x),max(out.mdata(j).x),512);
        yVec=out.mmod(j).fn(out.p,xVec);
        plot(xVec,yVec+offset*(j-1));
                
        subplot(1,2,2); hold on;        
        yVec=out.mmod(j).fn(out.p,x);
        em=.25;
        plot(x,offset*em*(j-1)+x*0,'k');% zero error line
        errorbar(x,y-yVec+offset*em*(j-1), sqrt(out.mdata(j).vary),'CapSize',1);
        title(sprintf('Errors, \\chi^2=%3f',out.chisq));        
    end
    subplot(1,2,1);
    phaseData=out.p(7:nParVar:end); freqData = out.p(6:nParVar:end); 
    phaseData = mod(phaseData,pi); 
    equivTime = phaseData./freqData;
    plot(equivTime,(1:length(phaseData))*offset+.15,'k.-');
    title('Fits + Phase')
    fixLog; 
    subplot(1,2,2); fixLog;
    fignum = fignum+1; 
    figure(fignum); clf; 
    subplot(2,2,1);     
    plot(out.params(1,:)/2/pi); ylabel('J (GHz)'); 
    subplot(2,2,2); 
    plot(phaseData); ylabel('Phase (rad)'); 
    subplot(2,2,3); 
    plot(out.params(4,:)); ylabel('Offset (V)'); 
    subplot(2,2,4); 
    plot(out.params(5,:)); ylabel('t offset (ns)'); 
end
end