function wf=rc_trafofn(wf, channel, tau)
% This function takes a waveform sampled at 1ns / point, and precompensates
% it for an RC rise time.  The argument is the rise time in ns.
%allows for different arguments on different channels

if size(wf,1)>1, wf=wf(channel,:); end
if length(tau)>1, tau=tau(channel); end
if tau == 0 || abs(tau) < 1e-3, return; end
downsamp=1;
y=wf;
for i=1:downsamp % smooths by twos.  
  y=([ y(1) y(floor(1+.5:0.5:end+.5)) ] + [ y(floor(1:0.5:end)) y(end) ])/2;
end
nyquist=(2^downsamp)*pi*1e9;  % 1 Ghz nyquist frequency on the resampled data
k= ([ 0:(length(y)/2-1) length(y)/2 -((length(y)/2-1):-1:1) ])*nyquist/(length(y)/2);
if tau > 0
  wf=real(ifft(fft(y)./rc(k,tau),'symmetric'));
else
  wf=real(ifft(fft(y).*rc(k,tau),'symmetric'));
end
for i=1:downsamp
  wf = (wf(1:2:end)+wf(2:2:end))/2;
end
end

function j=rc(k,tau)
% Skin effect transfer function. k=omega
   tau=abs(tau)*1e-9;
   j=exp(1i*tau*k/1.464)./(1i*k*tau+1);
end