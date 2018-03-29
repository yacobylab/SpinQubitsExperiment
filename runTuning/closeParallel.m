function closeVal = closeParallel(chans, lockin, zeroval)
% closes off parallel channel using list of gates chans, lockin channel. 
%function closeVal = closeParallel(chans, lockin, zeroval)
global smdata; 
keepgo = 1; newkeepgo = 1;
% tauChan = findChans(lockin,18);
% tau = smget(tauChan);
% smset(tauChan,0.3);
initval = -.05;
smset(chans,initval); currVal = initval; 
stopVal = smdata.channels(chl(chans{1})).rangeramp(1); 
while keepgo && currVal >= stopVal
    sminc(chans,-0.05); currVal = currVal -0.05; 
    pause(1.5);
    val1 = abs(cell2mat(smget(lockin)));
    val2 = abs(cell2mat(smget(lockin)));
    fprintf('Current is %3.3f nA \n',1e9/2*(val1+val2)); 
    oldkeepgo = newkeepgo;
    newkeepgo = val1 > zeroval || val2 > zeroval;
    if ~newkeepgo && ~oldkeepgo
        keepgo =0;
    end
end
closeVal = cell2mat(smget(chans{1})); 
fprintf('Close val is %2.2f \n',closeVal);
%smset(tauChan,cell2mat(tau)); 
end