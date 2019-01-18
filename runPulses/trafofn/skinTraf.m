function wf=skinTraf(wf, chan,atten, clock)
% This function takes a waveform sampled at 1ns/point, and precompensates
% it for skin effect damping.  Atten specifies the amount of damping in dB at 1 Ghz.

if ~exist('clock','var') || isempty(clock), clock=1e9; end
if atten ==0, return; end
alpha=atoa(atten,2*pi*clock);
downsamp=1;

% Double the sampling rate of the data
y=wf(chan,:);
for i=1:downsamp
  y=([ y(1) y(floor(1+.5:0.5:end+.5)) ] + [ y(floor(1:0.5:end)) y(end) ])/2;
end

nyquist=(2^downsamp)*pi*clock;  %1 Ghz nyquist frequency on the resampled data
k= ([ 0:(length(y)/2-1) length(y)/2 -((length(y)/2-1):-1:1) ])*nyquist/(length(y)/2);

if(atten > 0)
  wf=real(ifft(fft(y)./skin(k,alpha,clock),'symmetric'));
else
  wf=real(ifft(fft(y).*skin(k,alpha,clock),'symmetric'));
end
for i=1:downsamp
  wf = (wf(1:2:end)+wf(2:2:end))/2;
end
end

function j=skin(k,alpha,clock)
% Skin effect transfer function.  k=omega
%j=exp(-sqrt((1i*k+alpha).*(1i*k))).*exp(1i*k);  % Second term mostly cancels phase  
  kmax=2*pi*clock;
  j=exp(-(1+1i)*sqrt(k*alpha)) .* exp(1i*(k/kmax)*sqrt(kmax*alpha));
end

function alpha=atoa(atten,w)
  alpha = ((log(10)*atten/(20))^2)/w;    
end