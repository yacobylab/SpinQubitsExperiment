function closeVal = closeParallel(chans, lockin, zeroVal)
% Closes off parallel channels (for double-bonded ohmics). Give list of
% channels, chans, that will fully close all connected channels (see measureQPCauto
% for examples). 
% lockin is channel measuring current.
% Return value at which current shuts off. 
% Wait for current to be below zeroVal, default 
%function closeVal = closeParallel(chans, lockin, zeroval)
global smdata; 
keepGoing = 1; newkeepgo = 1;
if ~exist('zeroVal','var'), zeroVal = 0; end
initVal = -.05;
smset(chans,initVal); currVal = initVal; 
stopVal = smdata.channels(chl(chans{1})).rangeramp(1); % minimum value  
while keepGoing && currVal >= stopVal
    sminc(chans,-0.05); currVal = currVal -0.05; % Increment by 50 mV steps. 
    pause(1.5); % for equilibration
    val1 = abs(cell2mat(smget(lockin))); 
    val2 = abs(cell2mat(smget(lockin))); % average 2 values 
    fprintf('Current is %3.3f nA \n',1e9/2*(val1+val2)); 
    oldkeepgo = newkeepgo;
    newkeepgo = val1 > zeroVal || val2 > zeroVal;
    if ~newkeepgo && ~oldkeepgo
        keepGoing =0;
    end
end
closeVal = cell2mat(smget(chans{1})); 
fprintf('Close val is %2.2f \n',closeVal);
end