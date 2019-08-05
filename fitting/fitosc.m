function [fitpars,fitfn,x,res,mse,err]=fitosc(x,y,opts,rng,style)
% Fit decaying sinusoid. 
% function [fitpars,fitfn]=fitosc(x,y,opts,rng,style)
% opts: 
%   nodecay
%   afitdecay: 
%   center: 
% def'n fitpars: 
% for all: 1: offset 4: freq
% NoDecay: 2: cos coef, 3 sin coef
% Coefs: 2: cos coef, 3 sin coef, 5: decay center, y(6) 1/t2
% Phase: 2: amp, 3: phase, 5: center, 6: 1/t2
% rng: x-rng for data to fit. 
% If style given, plot with given style. 

if ~exist('opts','var'), opts = ''; end
beta0.fn = @fioscill; beta0.args = {1};
cosNoDecay = '@(y,x) y(1) + y(2) * cos(y(4)*x) + y(3) * sin(y(4)*x)';
cosCoefs =   '@(y,x) y(1) + (y(2) * cos(y(4) * x) + y(3) * sin(y(4) * x)).*exp(-(x-y(5)).^2 * y(6).^2)';
cosPhase =   '@(y,x) y(1) + y(2) * cos(y(4) * x + y(3)).*exp(-(x-y(5)).^2 * y(6).^2)';
if exist('rng','var') && ~isempty(rng)
    pts = x >= rng(1) & x <= rng(2);
    x = x(pts); y = y(pts);
end

beta0.args={2};
if isopt(opts,'nodecay') || (isopt(opts,'afitdecay') && std(y) < 2e-2) % No decay, centered
    fitpars=fitwrap('fine noplot',x,y,beta0, cosPhase, [1 0 1 1 0 0]); % Don't fit amplitude, include unfit decay 
    ig = [fitpars(1), fitpars(2)*cos(fitpars(3)), fitpars(2)*sin(fitpars(3)), fitpars(4:6)];
    [fitpars,res,~,~,mse,err]=fitwrap('fine noplot',x,y,ig, cosNoDecay, [1 1 1 1 0 0]);
    fitfn=str2func(cosNoDecay);
elseif isopt(opts,'center') && ~isopt(opts,'phase') % Decay and center    
    fitpars=fitwrap('fine noplot',x,y,beta0,cosCoefs, [1 1 1 1 0 0]); % Don't fit center or decay
    [fitpars,res,~,~,mse,err]=fitwrap('fine noplot',x,y,fitpars, cosCoefs, [1 1 1 1 0 1]);
    fitfn=str2func(cosCoefs);
elseif isopt(opts,'phase')
    fitpars=fitwrap('fine noplot',x,y,beta0,cosPhase, [1 0 1 1 0 0]); % Don't fit amp or decay    
    [fitpars,res,~,~,mse,err]=fitwrap('fine noplot',x,y,fitpars, cosPhase, [1 1 1 1 0 1]);
    fitfn=str2func(cosPhase);
else  % Decay but no center    
    fitpars=fitwrap('fine noplot',x,y,beta0,cosPhase, [1 0 1 1 0 0]); % Don't fit amp or decay
    fitpars = [fitpars(1), fitpars(2)*cos(fitpars(3)), -fitpars(2)*sin(fitpars(3)), fitpars(4:6)];
    [fitpars,res,~,~,mse,err]=fitwrap('fine noplot',x,y,fitpars, cosCoefs, [1 1 1 1 1 1]);
    fitfn=str2func(cosCoefs);
end
err = squeeze(err(:,:,1));
if isopt(opts,'pldata')    
    plot(x,y,'.-');    
end
if isopt(opts,'samecolor')
    a =gca;
    colorOrd = a.ColorOrderIndex;
    if colorOrd>1
        colorOrd = colorOrd-1;
    else
        colorOrd = 7;
    end
    a.ColorOrderIndex = colorOrd;
end
if ~isopt(opts,'noplot')
    if isopt(opts,'interp'), x=linspace(x(1,1),x(1,end),1024); end
    if exist('style','var') && ~isempty('style')
        plot(x,fitfn(fitpars,x),style);
    elseif isopt(opts,'red')
        plot(x,fitfn(fitpars,x),'r');
    else
        plot(x,fitfn(fitpars,x));
    end
end
end