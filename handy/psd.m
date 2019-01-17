function [sigf,freq] = psd(sig,Fs) 
% compute power spectral density. 
%function [sigf,freq] = psd(sig,Fs) 
% sig is linear signal 
% Fs is sampling frequency 
if ~exist('Fs','var') 
    Fs = 1; 
end
if round(length(sig)/2)~=length(sig)/2 % odd number of points 
    sig(end)=[];
end
if all(isnan(sig))
elseif any(isnan(sig))    
    sig = sig(1:find(isnan(sig),1)-1); 
end
N= length(sig); 
inds = 1:N/2+1; 
xdft = fft(sig); 
xdft = xdft(inds); 
sigf = (1 / (Fs*N))*abs(xdft).^2;
sigf(2:end-1) = 2*sigf(2:end-1);
freq = 0:Fs/N:Fs/2; 

end


