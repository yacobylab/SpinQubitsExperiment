function [fp,fitfn] = fioscill(x, y, ctrl)
% Initial guesses for parameters to fit oscillating data
% function fp = fioscill(x, y, ctrl)
% for both:
% fp(1): y offset, fp(4): frequency fp(5): shift (always 0), decay
% prefactor (very basic calc, 1/2 of inverse mean time)
% ctrl = 1: returns fp(2): cos coef, fp(3): s coef
% ctrl = 2: returns fp(2): global coef, fp(3): phase
% returns pars: offset, amplitude, phase, freq, shift, decay prefac = 1/mean(abs(x));
% Return fitfn as well. 

if ~exist('ctrl','var'), ctrl = 2; end
dx = mean(diff(x));
nx = length(x);
window=1;
% hack to get phase shift, factors determined experimentally
ft = fft((y-mean(y)).*window) .* exp(1i * x(1) * (0:nx-1) * 2*pi /(dx*nx))/nx;
ft = ft(1:round(end/2));
ftAbs = abs(ft); 
[~, freqInd] = max(ftAbs(2:end)); % Find frequency
xi = max(1,freqInd-2) : min(length(ft)-1,freqInd+2); % Average points around peak. 
meanxi = sum(ftAbs(xi+1).*xi)/sum(ftAbs(xi+1));
fp(4) = 2 * pi * meanxi / ((nx+1) * dx); % frequency
f = ft(freqInd+1);
fp(1) = mean(y);
fp(5) = 0;
fp(6) = 1 / (2 * mean(abs(x)));
switch ctrl
    case 1        
        fp(2) = 3*real(f);
        fp(3) = 3*imag(f);
        fitfn = '@(p,x) p(1) + (p(2).*cos(p(4) .* x) + p(3).*sin(p(4) .* x)).*exp(-(x-p(5)).^2.*p(6).^2)'; 
    case 2
        fp(3) = angle(f);
        fp(2) = 3*abs(f);
        fitfn = '@(p,x) p(1) + p(2).*cos(p(4) .* x + p(3)).*exp(-(x-p(5)).^2.*p(6).^2)'; 
end