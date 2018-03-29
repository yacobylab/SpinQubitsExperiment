startupsm

global awgdata; global plsdata;
global scandata; global tuneData; 
global smdata; global fbdata; 

%fprintf('Reset the RF freqs \n'); 


% % 2017-06-02: added to fix problem of dimension mismatch
% 
%         % Subscripted assignment dimension mismatch.
%         %Error in awgcntrl (line 89)
%         %             val(end+1) = fscanf(awgdata(a).awg,'%f');
%         %Error in awgcntrl (line 70)
%         %           if any(any(awgcntrl('israw')))
% 
% 
% awgdatasave = awgdata;
% clear awgdata
% 
% awgdata(1)= awgdatasave;
% global awgdata