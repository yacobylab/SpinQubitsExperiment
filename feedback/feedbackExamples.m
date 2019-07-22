%% Check if pumping works. This will pump, then measure gradient repeatedly 
characPump;  % singlet directoin 
characPump('',struct('opts','tl')); % triplet direction 
%% Measure, set gradient 
getgradient % 
setgradient 
%%
fbdata.params(1).setpoint = 30; % Set value for gradient 
%% Parameters to change in fb 
% Change to 32 or 64 if you are using charac pump and want to see better
% what is happening. 
fbdata.nloopfb = 16;
%%
rundBz('fb')
%% If stp/tl fits work, this will remake dictionary elements for feedback to be centered
tuneData.tl.updateGroup('target'); 
tuneData.stp.updateGroup('target'); 
%%
tuneData.stp.run
tuneData.tl.run
%% If STP/TL centered, update feedback groups
feedbackTest('all'); 
%% For centering, first update the groups so that STP/TL centered (see above) 
awgcntrl('on start raw err'); % start by turning on raw mode 
%% See how stp and tl have moved 
tuneData.stp.run
tuneData.tl.run
%% calculate how this moves the measurement point 
tuneData.tmp.run
tuneData.center('noconfirm')
%%
tuneData.dotCenter