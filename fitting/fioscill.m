function fp = fioscill(x, y, ctrl)
% Initial guesses for parameters to fit oscillating data
% function fp = fioscill(x, y, ctrl)
% for both:
% fp(1): y offset, fp(4): frequency fp(5): shift (always 0), decay
% prefactor (very basic calc, 1/2 of inverse mean time)
% ctrl = 1: returns fp(2): sin coef, fp(3): cos coef
% ctrl = 2: returns fp(2): global coef, fp(3): phase
% returns pars: offset, amplitude, phase, freq, shift, decay prefac = 1/mean(abs(x));

dx = mean(diff(x));
nx = length(x);
window=1;
% hack to get phase shift, factors determined experimentally
ft = fft((y-mean(y)).*window) .* exp(1i * x(1) * (0:nx-1) * 2*pi /(dx*nx))/nx;
ft = ft(1:round(end/2));
[~, freqInd] = max(abs(ft(2:end))); % Find frequency
xi = max(1,freqInd-2) : min(length(ft)-1,freqInd+2);
meanxi = sum(abs(ft(xi+1)).*xi)/sum(abs(ft(xi+1)));
fp(4) = 2 * pi * meanxi / ((nx+1) * dx);
f = ft(freqInd+1);
fp(1) = mean(y);
fp(5) = 0;
fp(6) = 1 / (2 * mean(abs(x)));
switch ctrl
    case 1
        fp(3) = -3*real(f);
        fp(2) = -3*imag(f);
    case 2
        fp(3) = angle(f);
        fp(2) = -3*abs(f);
end