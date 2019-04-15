function [fitpars,fitfn,fitstring]=fitdecay(x,y,opts,rng,style)
% Fit different types of exponential decay 
% options: gauss, both, power, fitoffset, plot
%   types of fits:
%   gauss: Fit to gaussian centered at 0
%   both: Fit to gaussian multiplied by exp. decay
%   power: fit to arbitrary power of decay
%   none of those given: 
%   fitoffset: allow the exponential not to decay to 0 
%   plot: plot the data on current axis. 
%   interp: use interpolated data for fit fn instead of x data. 
% rng: rng of x data to fit 
% style: style to plot data in 

if ~isopt(opts,'fitoffset')
    mask=[1 1 0];
else
    mask=[1 1 1];
end
if exist('rng','var') && ~isempty(rng) % Fit data within given x-rng. 
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