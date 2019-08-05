function data=autoscan(scanname,opts,config)
% Run tuning scans. Checks that scanning within ramprates, sets other channels correctly, sets sensor scan ranges to correspond to CB scan, etc.
% function autoscan(scanname,config)
% If a config or opt is given, need to also give the option 'run' if you want the scan to be run.
% Lockin scans across dot: CB, NT, qpc, hyst, leads, twoD, offT
% Lockin scans across sensor: SDLock, sensLock, sensorLock, twoDsens, junc
% DAQ scans: sens, sensor, RF, phase, SD, Calib, sensGate
% last (run previous scan again)
% configs:
%    sensorGate: change the trafofn sensor gate on all relevant scans
%    sensGet: change lockin getchan on all sensing lockin scans.
%    CBGet: change all lead scans' getchan.
%    init: start scanning: run scan on a/b gates, set opposite side channels to neg values, clear trafofn.
%    trafofn: load trafofn from different scan.
%    nudge : if pos, move trafofn 10 mv in pos direction, v.v. for neg.
%    copy: copies range of scan in config to scanname. maybe also want to copy trafo?
%    sweep: runs a sweep changing the value of a channel through some range for npoints. also give a config of ssetchan, srng, snpoints.
%    sensorVal: print out the optimal value for sensor gate given certain value for a gate.
% opts:
%    gofast: change ramprate to maximum.
%    gca: zoom in on scan, then call to change range.
%    fast: only 8 points in y, max ramprate.
%    rot: rotate 90 deg (check that this works!)
%    flip: first loop direction (not sure if this works)
%    hyst: check possible hyst by running 3 scans w/ rot and flip of current CB scan.
%    print: do print of scan plus scantime
%    trafo: run sensorTrafofn on scan, make 2nd setchan of scan the sensor setchan (change this to make it more general?)
%           trafa: do auto version of trafo
%    one/two: start collecting both Lockin channels or just one.
%    startSens: copy params from CB scan, set up trafo, set active scan to sensLock
%    startDAQ: copy params from sensLock scan, set set active scan to sens
%    load: load the range, trafofn, npoints from an old scan.
%    ldrng: load rng from old scan
%    ldtraf: load trafofn from old scan
%    lockin: configure lockins in usual way for charge sensing.
%    set: Set values of 2D scan to where you click.
%    csens: copy sensor sensGate range to sd, csd: copy sd range to sensor.
%    close: close off parallel channel, save value to scandata.
%    checkhy: run a set of hyst scans at different ramprates.
%    b: set b gates to clicked val, set CB scan to scan a gates. [largely defunct for current tuning]
%    up/ down : run a scan going up in rng or down. (Like frame up and down in AFM)
% can also edit many vals of scans this way, by specifying in config:
%    rng/npoints/setchan (with loop number)
%    rn11, rn12, rn21, rn22: change one value of rng. first number is loop, 2nd number is which rng val .
%    pointSpace (with loop number) : change npoints to get given spacing in mV
%    between points
%    ramprate : if using lock in, change npoints to get certain ramprate. (set ramptime to min for given int. time first)
%               if using DAQ, change ramptime.
%    ramptime : for loop 1
%    getchan (for loop 2)
%    move (with loop number) : shifts scan by some amt.
%    in opts, can also give finer1, finer2, coarser1, coarser2. This doubles / halves the number of points in given loop.
%    center: change the center of the scan. if you give an empty bracket, will have you click on screen.
%    setchan1/2 : cell or string.
%    diff1/2: magnitude of rng of scan, with same center as before.
% for running scans, some special features. Can turn them off by setting scandata.autoramp = 0
% if you give an opt, won't run scan unless you also give 'run' as an opt.
%   sens: charge sensing scan with a gates.
%   sensor : sensor scan of nearby a gate vs. current sensor Gate. Sets inner a gate to its mean value in activeScan. Sets rng/setchan of 1st
% loop to 2nd loop of active scan.
% phase: sweep phase of phase shifter on readout. Set phase to value where DAQ max.
% SD: sensor dot scan. Set a gates to their mean value for curr. activeScan
% RF: Sweep RF frequency vs. phase. On resonance, signal should go to 0 (if impedance matched). Sets a gates to their mean vals in active
%     scan. Sets freq to value for line with minimum average.
% Calib: run 1d scan of sensor top gate.
% sensGate:
% Lockin scans:
% sensLock: charge sensing scan with a gates
% SDLock: sensor dot scan. set a gates to their mean vals in activeScan.
% sensorLock: Run set inner a gate to mean val in activeScan. Set rng/setchan
% of 1st loop to 2nd loop of activeScan.
% CB: coulomb blockade wall-wall scan of a gates or a and b gates.
% NT: Nose and Tail swept against all A and B gates, to get a sense of how to set nose tail. probably phasing this out.
% qpc: After setting nose and tail, sweep each a/b gate separately, then a/b of each side together.
% hyst: Run buffered qpc sweep (tail + other gate) in each direction at a specific ramprate.
% junc: very small scan of a gates, usually of a signal junction.
% twoD: Run a 2D scan. Run with keepval to keep nonscanned channels at
% initial values (instead of setting to 0).
% offT
% twoDsens: can run with option keepval to set channels back to their initial values at the end of
% each scan, not to 0.

%% Configuring scan list.
global scandata; global smdata;
minRamp = 0.0625;
if strcmp(scanname,'last') || isempty(scanname)
    scanname = scandata.last;
end
if ~exist('opts','var'), opts = ''; end
if isstruct(opts) || iscell(opts)
    config = opts;
    opts = '';
end
if ~exist('config','var')
    config = ''; currConfigs = {''};
else
    if iscell(config), config=struct(config{:}); end
    currConfigs = fieldnames(config);
end

% check that any used configs , opts, scans actually exist
daqscans = {'sens','sensor','phase','SD','RF','Calib','sensGate','juncd','oneD','RFgate'};
lockscans = {'SDLock','sensLock','CB','sensorLock','NT','qpc','hyst','junc','twoD','offT','twoDsens','leads'};
optsList = {'init', 'gca','fast','rot','flip','hyst','print','trafo','trafa','bgates','one','two','close','startSens','load','sweep','lockin','copy','checkhy','run','ldtraf','ldrng','keepvals','set','csd','csens',...
    'finer1','finer2','coarser1','coarser2','fn','crs','sensorVal','up','down','startDAQ','gofast','tuneSens','trafx','slp'}; %longPrep???
configList = {'switch','sensorGate','sensGet','CBget','init','nudge','trafofn','copy','move1','move2','rng1','rng2','rn11','rn22','rn12','rn21','npoints1','npoints2','pointSpace1','pointSpace2','ramprate',...
    'setchan1','setchan2','ramptime','snpoints','srng','ssetchan','center','diff1','diff2','sensorVal'};
scanList = [daqscans(:); lockscans(:); {'last'}];
currOpts = strsplit(opts);

if ~any(strcmpi(scanname,scanList)) && ~isempty(scanname)
    error('%s is not a valid scan',scanname);
end
for i = 1:length(currOpts)
    if ~any(strcmpi(currOpts{i},optsList)) && ~isempty(currOpts{i})
        error('%s is not a valid option',currOpts{i});
    end
end
for i = 1:length(currConfigs)
    if ~any(strcmpi(currConfigs{i},configList)) && ~isempty(currConfigs{i})
        error('%s is not a valid config',currConfigs{i});
    end
end
if ~iscell(scandata.(scanname).loops(1).setchan)
    scandata.(scanname).loops(1).setchan = {scandata.(scanname).loops(1).setchan};
end
scandata.(scanname).data.rangeramp = smdata.channels(chl(scandata.(scanname).loops(1).setchan{1})).rangeramp(3);
taunum= 18; s=scandata.name(1); s = upper(s);
if isfield(scandata.(scanname).loops,'setchan') && ischar(scandata.(scanname).loops(1).setchan)
    scandata.(scanname).loops(1).setchan = {scandata.(scanname).loops(1).setchan}; 
end
if isfield(scandata.(scanname).loops,'setchan') && ischar(scandata.(scanname).loops(2).setchan)
    scandata.(scanname).loops(2).setchan = {scandata.(scanname).loops(2).setchan}; 
end
if isfield(scandata.(scanname).loops,'getchan') && ischar(scandata.(scanname).loops(2).getchan)
    scandata.(scanname).loops(2).getchan = {scandata.(scanname).loops(2).getchan}; 
end
if isfield(scandata.(scanname).loops,'getchan') && ischar(scandata.(scanname).loops(1).getchan)
    scandata.(scanname).loops(1).getchan = {scandata.(scanname).loops(1).getchan}; 
end
%% Run and edit scans
if (isempty(config) && isempty(opts)) || isopt(opts,'run')
    scandata.last = scanname;
    % Check ramptime. If it is too small compared to averaging time, set to 1.5x tau, then round it to an acceptable time for the lock in.
    if any(strcmp(scanname,lockscans)) 
        oldtime = scandata.(scanname).loops(1).ramptime;
         % Lock in must be exp of 2. Times will typically be 62.5 ms (16 Hz), 0.0156 ms (64 Hz)
         scandata.(scanname).loops(1).ramptime = -2^(round(log2(abs(oldtime))));
         if oldtime ~=scandata.(scanname).loops(1).ramptime
             fprintf('Reset ramptime to %3.3f \n', scandata.(scanname).loops(1).ramptime);
         end
         try
             tauchan = findChans(scandata.(scanname).loops(2).getchan{1}, taunum); % Averaging time for lock in
             tau = cell2mat(smget(tauchan));             
             if abs(scandata.(scanname).loops(1).ramptime)< 1.5 * tau
                 oldTime = abs(scandata.(scanname).loops(1).ramptime);
                 scandata.(scanname).loops(1).ramptime= - 1.5 * tau;
                 scandata.(scanname).loops(1).ramptime = -2^(round(log2(abs(scandata.(scanname).loops(1).ramptime))));
                 warning('Ramping too fast for averaging time. Slowing down from %3.3f to %3.3f',oldTime,abs(scandata.(scanname).loops(1).ramptime));
             end
         catch
             warning('Not using a real Lock in channel or Tau doesn''t exist');
         end         
    end
    if ~isempty(scanname) && isfield(scandata,scanname) % Check ramprate to see if it's within correct range.
        [isgood,rng,rate] = smramprate(scandata.(scanname));
        if rate <= 1/2*scandata.(scanname).data.rangeramp && abs(scandata.(scanname).loops(1).ramptime) >= 2*minRamp && any(strcmp(scanname,lockscans))
            fprintf('Going slow. Consider reducing ramptime \n');
        end
        if ~isgood
            rat = abs(diff(rng) / diff(scandata.(scanname).loops(1).rng));
            newramptime = scandata.(scanname).loops(1).ramptime / (rat-0.05);
            if any(strcmp(scanname,lockscans))
                scandata.(scanname).loops(1).ramptime = -2^(ceil(log2(abs(newramptime))));
            else
                scandata.(scanname).loops(1).ramptime = -round(abs(newramptime),3)-1e-3;
            end
            warning('Ramping too fast. Increasing ramptime to %3.3f', scandata.(scanname).loops(1).ramptime);
        end
    end
    switch scanname
        case 'RF'
            if scandata.autoramp
                smset(scandata.(scandata.activeScan).loops(1).setchan{1},mean(scandata.(scandata.activeScan).loops(1).rng));
                smset(scandata.(scandata.activeScan).loops(2).setchan{1},mean(scandata.(scandata.activeScan).loops(2).rng));
            end
            oldVal = cell2mat(smget(scandata.RF.loops(2).setchan)); oldPhaseVal = cell2mat(smget(scandata.RF.loops(1).setchan));
            data=smrun(scandata.RF,smnext(sprintf('RF%s',s)));
            meanVal = max(abs(data{1}'));
            [~,ind] = min(abs(meanVal)); % you should really make sure it's small enough -- add a check.
            freqVals = scanRng(scandata.RF,2);
            figure(18); plot(freqVals/1e6,meanVal); xlabel('Frequency (MHz)');
            ylabel('Mean abs(Phase)');
            smset(scandata.RF.loops(2).setchan,freqVals(ind));
            fprintf('Setting %s to %4.4g MHz from %4.4g MHz \n',scandata.RF.loops(2).setchan{1},freqVals(ind),oldVal);
            
            phaseData = data{1}(ind,:);
            [~,indP]=max(abs(phaseData));
            phaseVals = linspace(scandata.(scanname).loops(1).rng(1),scandata.(scanname).loops(1).rng(2),scandata.(scanname).loops(1).npoints);
            smset(scandata.(scanname).loops(1).setchan,phaseVals(indP));
            fprintf('Setting %s to %3.3g from %3.3g V \n',scandata.(scanname).loops(1).setchan{1},phaseVals(indP),oldPhaseVal);
            fprintf('Mean val changes: %3.3f to %3.3f \n',data{1}(ind),data{1}(indP))
        case 'RFgate'
            if scandata.autoramp
                smset(scandata.(scandata.activeScan).loops(1).setchan{1},mean(scandata.(scandata.activeScan).loops(1).rng));
                smset(scandata.(scandata.activeScan).loops(2).setchan{1},mean(scandata.(scandata.activeScan).loops(2).rng));
            end
            oldVal = cell2mat(smget(scandata.RFgate.loops(2).setchan)); 
            oldPhaseVal = cell2mat(smget(scandata.RFgate.loops(1).setchan));
            data=smrun(scandata.RFgate,smnext(sprintf('RFgate%s',s)));
            meanVal = max(abs(data{1}'));
            [~,ind] = min(abs(meanVal)); % you should really make sure it's small enough -- add a check.
            freqVals = scanRng(scandata.RFgate,2);
            figure(18); plot(freqVals/1e6,meanVal); xlabel('Gate Value (V)');
            ylabel('Mean abs(Phase)');
            smset(scandata.RFgate.loops(2).setchan,freqVals(ind));
            fprintf('Setting %s to %4.4g V from %4.4g V \n',scandata.RFgate.loops(2).setchan{1},freqVals(ind),oldVal);
            
            phaseData = data{1}(ind,:);
            [~,indP]=max(abs(phaseData));
            phaseVals = linspace(scandata.(scanname).loops(1).rng(1),scandata.(scanname).loops(1).rng(2),scandata.(scanname).loops(1).npoints);
            smset(scandata.(scanname).loops(1).setchan,phaseVals(indP));
            fprintf('Setting %s to %3.3g from %3.3g V \n',scandata.(scanname).loops(1).setchan{1},phaseVals(indP),oldPhaseVal);
            fprintf('Mean val changes: %3.3f to %3.3f \n',data{1}(ind),data{1}(indP))
        case 'phase'
            setChan = scandata.(scanname).loops(1).setchan{1}; 
            oldVal = cell2mat(smget(setChan));
            data=smrun(scandata.(scanname));
            phaseData = data{1};
            [maxVal,ind]=max(abs(phaseData));
            phaseVals = scanRng(scandata.(scanname),1);
            [~,oldInd] = min(abs(phaseVals - oldVal));
            oldMax = phaseData(oldInd);            
            smset(setChan,phaseVals(ind));
            fprintf('Setting %s to %3.3g V from %3.3g V, from %3.3f to %3.3f\n',setChan,phaseVals(ind),oldVal,oldMax,maxVal);
        case 'sensGate'
            oldVal = cell2mat(smget(scandata.(scanname).loops(2).setchan));
            data=smrun(scandata.(scanname),smnext(sprintf('sensGate%s',s)));
            meanVal = nanmean(abs(data{1}'));
            [~,ind] = min(abs(meanVal));
            gateVals = linspace(scandata.(scanname).loops(2).rng(1),scandata.(scanname).loops(2).rng(2),scandata.(scanname).loops(2).npoints);
            %figure(18); plot(gatevals,mnval);
            smset(scandata.(scanname).loops(2).setchan,gateVals(ind));
            fprintf('Setting %s to %3.3f from %3.3f \n',scandata.(scanname).loops(2).setchan,gateVals(ind),oldVal);
        case 'Calib'
            smrun(scandata.(scanname));
        case 'sens'
            data=smrun(scandata.sens,smnext(sprintf('sens%s',s)));
        case 'SDLock'
            if scandata.autoramp
                smset(scandata.(scandata.activeScan).loops(1).setchan{1},mean(scandata.(scandata.activeScan).loops(1).rng));
                smset(scandata.(scandata.activeScan).loops(2).setchan{1},mean(scandata.(scandata.activeScan).loops(2).rng));
            end
            data=smrun(scandata.SDLock,smnext(sprintf('SD%sLock',s)));
        case 'SD'
            if scandata.autoramp
                smset(scandata.(scandata.activeScan).loops(1).setchan{1},mean(scandata.(scandata.activeScan).loops(1).rng));
                smset(scandata.(scandata.activeScan).loops(2).setchan{1},mean(scandata.(scandata.activeScan).loops(2).rng));
            end
            data=smrun(scandata.SD,smnext(sprintf('SD%s',s)));
            diffData = diff(data{1},1,2);
            [maxCol,ind2]=max(abs(diffData));
            [~,ind1] = max(abs(maxCol));
            ind2 = ind2(ind1);
            loop1Vals = scanRng(scandata.(scanname),1);
            loop2Vals = scanRng(scandata.(scanname),2);
            smset(scandata.(scanname).loops(1).setchan{1},loop1Vals(ind1));
            smset(scandata.(scanname).loops(2).setchan{1},loop2Vals(ind2));
        case 'sensLock'
            data=smrun(scandata.sensLock,smnext(sprintf('sensLock%s',s)));
        case 'junc'
            data=smrun(scandata.junc,smnext(sprintf('junc%s',s)));
        case 'juncd'
            if scandata.autoramp
                scandata.juncd.loops(2).trafofn = scandata.sens.loops(2).trafofn;
                scandata.juncd.loops(2).setchan{2} = scandata.sens.loops(2).setchan{2}; 
            else
                scandata.juncd.loops(2).setchan = scandata.juncd.loops(2).setchan(1); 
            end
            data=smrun(scandata.juncd,smnext(sprintf('junc%s',s)));
        case 'oneD'
            smset(scandata.sensor.loops(2).setchan,scandata.config.sensVals(1));
            smset(scandata.sensor.loops(1).setchan,scandata.config.sensVals(2));
            scandata.(scanname).loops(1).rng = scandata.sens.loops(1).rng;
            scandata.(scanname).loops(1).setchan = scandata.sens.loops(1).setchan(1);
            scandata.(scanname).loops(1).npoints = scandata.sens.loops(1).npoints;
            scandata.(scanname).loops(1).ramptime = scandata.sens.loops(1).ramptime;
            data=smrun(scandata.oneD);
            xvals = scanRng(scandata.oneD);
            dataDiff = smooth(diff(data{1}),10);
            maxSlp = abs(max(dataDiff)./(xvals(2)-xvals(1)))*4/5;
            scandata.config.outSlp = maxSlp;
        case 'CB'
            data=smrun(scandata.CB,smnext(sprintf('Wall%s',s)));
        case 'leads'
            data=smrun(scandata.leads,smnext(sprintf('Leads%s',s)));
        case 'twoD'
            smset(scandata.config.closeGates,scandata.config.closeVal)
            if strcmp(scandata.name,'right')
                chanList = {'4a','3a','4b','4a','N34'};
            else
                chanList = {'1a','2a','1b','2b','N12'};
            end
            scandata.(scanname).loops(1).rng = scandata.hyst.loops(1).rng; % or do you want to make this more negative
            scandata.(scanname).loops(2).rng = scandata.hyst.loops(1).rng;
            scandata.(scanname).loops(2).setchan = scandata.config.tgate;
            autoscan(scanname,'fs');
            for i = 1:length(chanList)
                scandata.(scanname).loops(1).setchan = chanList(i);
                data=smrun(scandata.(scanname),smnext(sprintf('twoD%s_%s_%s',s,scandata.(scanname).loops(1).setchan{1},scandata.(scanname).loops(2).setchan{1})));
                smset(setdiff({smdata.channels(1:18).name},scandata.config.closeGates),0);
            end
        case 'offT'
            smset(scandata.config.closeGates,scandata.config.closeVal)
            if strcmp(scandata.name,'right')
                offList = {'4a','4b'; '4b','N34'; '3a','3b'; '3b','N34'};
            else
                offList = {'1a','1b'; '1b','N12'; '2a','2b'; '2b','N12'};
            end
            scandata.(scanname).loops(1).rng = scandata.hyst.loops(1).rng; % you will have to check if this is right
            scandata.(scanname).loops(2).rng = scandata.hyst.loops(1).rng;
            
            for i = 1:size(offList,1)
                smset(scandata.config.tgate,scandata.config.tval);
                scandata.(scanname).loops(1).setchan = offList(i,1);
                scandata.(scanname).loops(2).setchan = offList(i,2);
                data=smrun(scandata.(scanname),smnext(sprintf('offT%s_%s_%s',s,scandata.(scanname).loops(1).setchan{1},scandata.(scanname).loops(2).setchan{1})));
                smset(setdiff({smdata.channels(1:18).name},scandata.config.closeGates),0);
            end
        case 'twoDsens'
            if strcmp(scandata.name,'right')
                chanList = {'SD4top','4b';'SD4mid','4a';'SD4bot','4a'};
            else
                chanList = {'SD1top','1b';'SD1mid','1a';'SD1bot','1a'};
            end
            if isopt(opts,'keepvals')
                chanvals = cell2mat(smget(1:18));
            else
                chanvals = 0;
            end
            scandata.(scanname).loops(1).rng = [0 -.4]; % or do you want to make this more negative?
            scandata.(scanname).loops(2).rng = [0 -.4];
            for i = 1:size(chanList,1)
                scandata.(scanname).loops(2).setchan = chanList(i,1);
                scandata.(scanname).loops(1).setchan = chanList(i,2);
                data=smrun(scandata.(scanname),smnext(sprintf('twoD%s_%s_%s',s,scandata.(scanname).loops(1).setchan{1},scandata.(scanname).loops(2).setchan{1})));
                smset(1:18,chanvals);
            end
        case 'sensor'
            if scandata.autoramp
                smset(scandata.(scandata.activeScan).loops(1).setchan{1},mean(scandata.(scandata.activeScan).loops(1).rng));
            end
            if ~isempty(scandata.activeScan) && scandata.autoramp
                scandata.sensor.loops(1).rng = scandata.(scandata.activeScan).loops(2).rng;
                %scandata.sensor.loops(1).rng(1) = scandata.sensor.loops(1).rng(1)+ 0.015;
                scandata.sensor.loops(1).rng = fliplr(scandata.sensor.loops(1).rng);
                scandata.sensor.loops(1).setchan = scandata.(scandata.activeScan).loops(2).setchan(1);
            end
            fname = smnext(sprintf('sensor%s',s));
            data=smrun(scandata.sensor,fname);
            if isopt(opts,'slp')
                [~,out] = plotChrg('maxSlp invis',{['sm_' fname]});
                scandata.config.inSlp = out.sens;
                scandata.config.sensVals = [out.yval,out.xval];
            end
        case 'sensorLock'
            if scandata.autoramp
                smset(scandata.(scandata.activeScan).loops(1).setchan{1},mean(scandata.(scandata.activeScan).loops(1).rng));
            end
            if ~isempty(scandata.activeScan) && scandata.autoramp
                scandata.sensorLock.loops(1).rng = scandata.(scandata.activeScan).loops(2).rng;
                scandata.sensorLock.loops(1).rng(1) = scandata.sensorLock.loops(1).rng(1)+ 0.015;
                scandata.sensorLock.loops(1).setchan = scandata.(scandata.activeScan).loops(2).setchan(1);
            end
            
            % lockChan = findChans(scandata.(scanname).loops(2).getchan{1},6);
            % vout = cell2mat(smget(lockChan));
            %div = scandata.config.div(1) / (scandatscandata.sensorLock.loops(2)a.config.div(1) + scandata.config.div(2));
            %scandata.(scanname).data.div = div; scandata.(scanname).data.vout = vout;
            data=smrun(scandata.sensorLock,smnext(sprintf('sensorLock%s',s)));
        case 'qpc'
            smset(scandata.config.tgate,scandata.config.tval); smset(scandata.config.ngate,scandata.config.nval);
            data = smrun(scandata.qpc,smnext(sprintf('oneGate%s',s)));
            gateInner = data{1}(5,:); gateOuter = data{1}(6,:);
            gateVals = linspace(scandata.qpc.loops(1).rng(1), scandata.qpc.loops(1).rng(2), scandata.qpc.loops(1).npoints);
            inZero = find(gateInner <1e-10);
            if length(inZero)>=5
                inV = gateVals(inZero(5)) - 0.035;
            else
                inV = min(gateVals);
            end
            outZero = find(gateOuter <1e-10);
            if length(outZero)>=5
                outV = gateVals(outZero(5)) - 0.035;
            else
                outV = min(gateVals);
            end
            if outV < -.6
                scandata.config.outRng = [outV outV+0.6];
            else
                scandata.config.outRng = [outV 0];
            end
            if inV < -.6
                scandata.config.inRng = [inV inV+0.6];
            else
                scandata.config.inRng = [inV 0];
            end
        case 'NT'
            smset(scandata.config.closeGates,scandata.config.closeVal)
            data = smrun(scandata.NT,smnext(sprintf('NT%s',s)));
            tnval = ginput(1);
            tval = tnval(2); %nval = tval;
            scandata.config.tval = tval; %scandata.config.nval = nval;
        case 'hyst'
            smset(scandata.config.closeGates,scandata.config.closeVal)
            scandata.hyst.loops(1).rng = [0 scandata.config.hystVal];
            if isopt(opts,'checkhy')
                rampRates = [300 240 120 60 30 10]*1e-3;
                for i = 1:length(rampRates)
                    smset('LockinTauA',0.01);
                    autoscan('hyst','',struct('ramprate',rampRates(i)));
                    autoscan('hyst');
                end
                %                 smset('LockinTauA',0.3);
                %                 autoscan('hyst','',struct('ramprate',1e-3));
                %                 scan = scandata.hyst; scan.loops(3).npoints = 1;
                %                 scan.loops(2).procfn(2).dim = scan.loops(1).npoints;
                %                 smrun(scan,smnext(sprintf('hyst%sPos',s)));
                %                 scan.loops(1).rng = fliplr(scan.loops(1).rng);
                %                 data = smrun(scan,smnext(sprintf('hyst%sNeg',s)));
                %                 smset('LockinTauA',0.03);
            else
                scandata.hyst.loops(2).procfn(2).dim = scandata.hyst.loops(1).npoints;
                scandata.hyst.loops(1).rng = sort(scandata.hyst.loops(1).rng);
                smrun(scandata.hyst,smnext(sprintf('hyst%sPos',s)));
                scandata.hyst.loops(1).rng = fliplr(scandata.hyst.loops(1).rng);
                data = smrun(scandata.hyst,smnext(sprintf('hyst%sNeg',s)));
            end
        case 'vanderpauw'
            % The best choice when there is no switch is to use OL1 to OS4top for current, OL2 to OS4bot for voltage.
            % Using only sensors, can do: OS1top to OS4bot current and voltage OS1bot to OS4top
            smaSwitch('','ground')
            smset('LockinTauA',0.3); smset('LockinTauB',0.3);
            smaSRS830('LockinB',struct('ref','ext'));
            smset(scandata.config.closeGates,scandata.config.closeVal);
            switchChans = smdata.inst(inl('Switch')).data.outChans;
            permsA = [1 3 2 4; 3 1 4 2; 2 4 1 3; 4 2 3 1];
            permsB = [3 4 2 1; 4 3 1 2; 2 1 3 4; 1 2 4 3];
            for i=1:size(permsA,1)
                smaSwitch({switchChans(permsA(i,1)),'VA'; switchChans(permsA(i,2)),'IA'; switchChans(permsA(i,3)),'VB';switchChans(permsA(i,4)),'IB'});
                pause(3);
                imn(i) = cell2mat(smget('LockinA')); vop(i) = cell2mat(smget('LockinB'));
                imny(i) = cell2mat(smget('LockinY')); vopy(i) = cell2mat(smget('LockinBY'));
                smaSwitch({switchChans(permsB(i,1)),'VA'; switchChans(permsB(i,2)),'IA'; switchChans(permsB(i,3)),'VB';switchChans(permsB(i,4)),'IB'});
                pause(3);
                ino(i) = cell2mat(smget('LockinA')); vpm(i) = cell2mat(smget('LockinB'));
                inoy(i) = cell2mat(smget('LockinY')); vpmy(i) = cell2mat(smget('LockinBY'));
            end
            rmnop = mean(vop./imn); rnopm = mean(vpm./ino);
            func = @(f) abs(cosh (( rmnop/rnopm -1)/(rmnop/rnopm+1)*log(2)./f)-1/2*exp(log(2)./f));
            f=fminunc(func,0.5);
            rho = pi / log(2) * (rmnop + rnopm)/2 * f;
            smset('LockinTauA',0.03); smset('LockinTauB',0.03);
            fprintf('Sheet resistance is %3.3f Ohms/Square \n',rho);
            scandata.config.rho = rho;
    end
end % Running scan

if exist('config','var') && ~isempty(config) && ~isempty(scanname) % Configuring scan params
    if isfield(config,'rng1')
        scandata.(scanname).loops(1).rng = config.rng1;
    end
    if isfield(config,'rn11')
        scandata.(scanname).loops(1).rng(1) = config.rn11;
    end
    if isfield(config,'rn12')
        scandata.(scanname).loops(1).rng(2) = config.rn12;
    end
    if isfield(config,'rn21')
        scandata.(scanname).loops(2).rng(1) = config.rn21;
    end
    if isfield(config,'rn22')
        scandata.(scanname).loops(2).rng(2) = config.rn22;
    end
    if isfield(config,'ramptime')
        scandata.(scanname).loops(1).ramptime = -abs(config.ramptime);
    end
    if isfield(config,'npoints1')
        scandata.(scanname).loops(1).npoints = config.npoints1;
    end
    if isfield(config,'rng2')
        scandata.(scanname).loops(2).rng = config.rng2;
    end
    if isfield(config,'move1')
        scandata.(scanname).loops(1).rng = scandata.(scanname).loops(1).rng + config.move1;
    end
    if isfield(config,'move2')
        scandata.(scanname).loops(2).rng = scandata.(scanname).loops(2).rng + config.move2;
    end
    if isfield(config,'ramptime2')
        scandata.(scanname).loops(2).ramptime = config.ramptime2;
    end
    if isfield(config,'npoints2')
        scandata.(scanname).loops(2).npoints = config.npoints2;
    end
    if isfield(config,'pointSpace1')
        scandata.(scanname).loops(1).npoints = ceil(abs(diff(scandata.(scanname).loops(1).rng)) / config.pointSpace1);
    end
    if isfield(config,'pointSpace2')
        scandata.(scanname).loops(2).npoints = ceil(abs(diff(scandata.(scanname).loops(2).rng)) / config.pointSpace2);
    end
    if isfield(config,'ramprate')
        if any(strcmp(scanname,lockscans)) % For the lock in, probably need to change npoints
            tauchan = findChans(scandata.(scanname).loops(2).getchan{1}, taunum);
            tau = cell2mat(smget(tauchan));
            scandata.(scanname).loops(1).ramptime= - 1.5 * tau;
            scandata.(scanname).loops(1).ramptime = -2^(round(log2(abs(scandata.(scanname).loops(1).ramptime))));
            scandata.(scanname).loops(1).npoints = ceil(abs(diff(scandata.(scanname).loops(1).rng)/scandata.(scanname).loops(1).ramptime/config.ramprate));
            fprintf('New ramptime %3.3f, npoints %d \n',scandata.(scanname).loops(1).ramptime, scandata.(scanname).loops(1).npoints);
        else % For daq, change ramptime
            scandata.(scanname).loops(1).ramptime = -abs(diff(scandata.(scanname).loops(1).rng)/scandata.(scanname).loops(1).npoints/config.ramprate);
        end
    end
    if isfield(config,'getchan')
        scandata.(scanname).loops(2).getchan = {config.getchan};
    end
    if isfield(config,'setchan1')
        if iscell(config.setchan1)
            scandata.(scanname).loops(1).setchan = config.setchan1;
        else
            scandata.(scanname).loops(1).setchan = {config.setchan1};
        end
    end
    if isfield(config,'setchan2')
        if iscell(config.setchan2)
            scandata.(scanname).loops(2).setchan = config.setchan2;
        else
            scandata.(scanname).loops(2).setchan = {config.setchan2};
        end
    end
    if isfield(config,'diff1')
        scandata.(scanname).loops(1).rng = mean(scandata.(scanname).loops(1).rng)+config.diff1/2*[-1 1];
    end
    if isfield(config,'diff2')
        scandata.(scanname).loops(2).rng = mean(scandata.(scanname).loops(2).rng)+config.diff2/2*[-1 1];
    end
    if isfield(config,'center') % Center scan at point clicked or arg given.
        if isempty(config.center)
            figure(1000);
            center = ginput(1);
        else
            center = config.center;
        end
        if length(center)==2
            rng1 = diff(scandata.(scanname).loops(1).rng);
            rng2 = diff(scandata.(scanname).loops(2).rng);
            scandata.(scanname).loops(1).rng = center(1) + [rng1/2,-rng1/2];
            scandata.(scanname).loops(2).rng = center(2) + [rng2/2,-rng2/2];
        else
            rng1 = diff(scandata.(scanname).loops(1).rng);
            scandata.(scanname).loops(1).rng = center(1) + [rng1/2,-rng1/2];
        end
    end
end
%% Configuring scans, changing instr. 

% start scanning: run scan on a/b gates, set opposite side channels to neg values, clear trafofn.
if isopt(opts,'init') 
    if scandata.config.switchOn
        smaSwitch({'OL1','VA'; 'OL2','IA'; 'OS1top','open';'OS1bot','open'})
    end
    if strcmp(scandata.name,'right')
        swapscan('right');
        fprintf('Initializing right side. Closing all left hand gates. \n')
        %         scandata.CB.loops(1).setchan = {'3a','3b'};
        scandata.CB.loops(1).setchan = {'3a'};
        %         scandata.CB.loops(2).setchan = {'4a','4b'};
        scandata.CB.loops(1).setchan = {'4a'};
        smset({'SD4top','SD4bot'},0)
    else
        swapscan('left')
        fprintf('Initializing left side. Closing all right hand gates. \n')
        %         scandata.CB.loops(1).setchan = {'2a','2b'};
        scandata.CB.loops(1).setchan = {'2a'};
        scandata.CB.loops(2).setchan = {'1a'};
        %         scandata.CB.loops(2).setchan = {'1a','1b'};
        smset({'SD1top','SD1mid','SD1bot'},0)
    end
    scandata.CB.loops(1).npoints = 128;
    scandata.CB.loops(2).npoints = 32;
    smset(scandata.config.closeGates,scandata.config.closeVal);
    %smset( scandata.config.ngate,scandata.config.nval);
    smset(scandata.config.tgate,scandata.config.tval);
    fprintf('Current getchan is %s. \n',scandata.CB.loops(2).getchan{1});
    fprintf('Setting averaging time to 0.03 s \n')
    tauchan = findChans(scandata.(scanname).loops(2).getchan{1}, taunum);
    smset(tauchan,0.03);
    fprintf('Clearing trafofn, setting settle time, ramping a and b gates \n')
    scandata.CB.loops(2).trafofn =[];
    scandata.CB.loops(1).ramptime = -0.0625;
    if isfield(scandata.config,'inRng') && ~isempty(scandata.config.inRng)
        scandata.CB.loops(1).rng = scandata.config.inRng;
        scandata.CB.loops(2).rng = scandata.config.outRng;
    else
        scandata.CB.loops(1).rng = [-.2 -.9];
        scandata.CB.loops(2).rng = [-.2 -.9];
    end
    scandata.CB.loops(1).settle = 1;
    autoscan('CB','fs');
    scandata.activeScan = 'CB';
    autoscan('CB');
end
if isopt(opts,'startSens') % copy params from CB scan, set up trafo, set active scan to sensLock
    fprintf('Copying CB range, setting point spacing to 1 mV, setting up trafofn, changing active scan \n');
    autoscan('sensLock',struct('copy','CB'));
    autoscan('sensLock',struct('pointSpace1',1e-3,'ramptime',-0.0625));
    autoscan('sensorLock',struct('setchan2','SD1mid','rng2',[-.4 -.8]));
    scandata.activeScan = 'sensLock';
end
if isopt(opts,'startDAQ') %  copy params from sensLock scan, set set active scan to sens
    scandata.sens.loops(1).rng = scandata.sensLock.loops(1).rng;
    scandata.sens.loops(2).rng = scandata.sensLock.loops(2).rng;
    scandata.sensor.loops(2).rng = scandata.sensorLock.loops(2).rng;
    scandata.activeScan = 'sens';
end

if isopt(opts,'two') % start collecting both Lockin channels 
    g1 = scandata.(scanname).loops(2).getchan{1};
    g2 = setdiff(scandata.config.getchans,g1);
    scandata.(scanname).loops(2).getchan(2) = g2;
    scandata.(scanname).disp(3:4) = scandata.(scanname).disp(1:2);
    scandata.(scanname).disp(3).channel = 2;
    scandata.(scanname).disp(4).channel = 2;
    fprintf('Adding other getchan \n');
end  
if isopt(opts,'one') % Change to collect from one lock in channel 
    scandata.(scanname).loops(2).getchan = scandata.(scanname).loops(2).getchan(1);
    scandata.(scanname).disp(3:4) = [];
    fprintf('Removing second getchan \n');
end
%% Tuning up sensor
if isopt(opts,'tuneSens') % Run sensor scan, then fit trafofn, then run sensing scan. 
    if strcmp(scandata.activeScan,'sens')
        autoscan('sensor');
        autoscan('sens','trafa');
        autoscan('sens');
    else
        autoscan('sensorLock');
        autoscan('sensLock','trafa');
        autoscan('sensLock');
    end
end
if isopt(opts,'ldtraf') % load trafofn from old scan
    file = uigetfile('');
    d = load(file);
    scandata.(scanname).loops(2).trafofn = d.scan.loops(2).trafofn;
    fprintf('Loaded trafofn from %s \n',file(1:end-4));
end
if isopt(opts,'trafo') || isopt(opts,'trafa') || isopt(opts,'trafx')
    if ~(strcmp(scanname,'sens') || strcmp(scanname,'sensLock'))
        warning('Adding trafofn to weird scan. Are you sure you want to do this?');
    end
    if strcmp(scanname,'sens') && ~isopt(opts,'trafx')
        scandata.(scanname).loops(2).setchan{2} = scandata.sensor.loops(2).setchan{1};
        fprintf('Setting trafofn chan to %s \n', scandata.(scanname).loops(2).setchan{2});
    elseif strcmp(scanname,'sensLock') && ~isopt(opts,'trafx')
        scandata.(scanname).loops(2).setchan{2} = scandata.sensorLock.loops(2).setchan{1};
        fprintf('Setting trafofn chan to %s \n', scandata.(scanname).loops(2).setchan{2});
    end
    if isopt(opts,'trafo')
        scandata.(scanname) = sensorTrafofn(scandata.(scanname)');
    elseif isopt(opts,'trafa')
        scandata.(scanname) = sensorTrafofn(scandata.(scanname),'','auto');
    else
        scandata.(scanname).loops(1).setchan{2} = scandata.sensor.loops(2).setchan{1};
        scandata.(scanname) = sensorTrafofn(scandata.(scanname),'','autox');
    end
end
if isfield(config,'sensorVal') % print out the optimal value for sensor gate given certain value for a gate.
    if isfield(scandata.sens.loops,'trafofn') && length(scandata.sens.loops(2).trafofn) > 1
        data=scandata.sens.loops(2).trafofn(2).fn(config.sensorVal,0,scandata.sens.loops(2).trafofn(2).args{:});
    else
        data=scandata.sens.loops(1).trafofn(2).fn(config.sensorVal,0,scandata.sens.loops(1).trafofn(2).args{:});
    end
end
%% Changing get or set chans in multiple scans.

if isfield(config,'CBget') % Change all lead scans' getchan.
    fprintf('Changing getchan for CB, NT, qpc \n')
    scandata.CB.loops(2).getchan = {config.CBget};
    scandata.NT.loops(2).getchan = {config.CBget};
    scandata.qpc.loops(2).getchan = {config.CBget};
end
if isfield(config,'sensGet') % Change sensor scans' getchan 
    fprintf('Changing getchan for SD, sensor, sens \n')
    scandata.SDLock.loops(2).getchan = {config.sensGet};
    scandata.sensLock.loops(2).getchan = {config.sensGet};
    scandata.sensorLock.loops(2).getchan = {config.sensGet};
end

% Add/change the sensing gate to all sensing scans. 
% With config = rm, remove sensing gate
if isfield(config,'sensorGate') 
    if strcmp(config.sensorGate,'rm')
        scandata.sensLock.loops(2).setchan = scandata.sensLock.loops(2).setchan{1};
        fprintf('Removing sensor from sens scan \n');
    else
        scandata.sensLock.loops(2).setchan{2} = config.sensorGate;
        scandata.sens.loops(2).setchan{2} = config.sensorGate;
        scandata.juncd.loops(2).setchan{2} = config.sensorGate;
        scandata.SD.loops(2).setchan{1} = config.sensorGate;
        scandata.sensor.loops(2).setchan{1} = config.sensorGate;
        scandata.sensorLock.loops(2).setchan{1} = config.sensorGate;
        fprintf('Changing setchan in sens, SD, sensor scans to %s \n', config.sensorGate);
    end
end
%% Set values in scan and copying values from other scans. 
if isopt(opts, 'bgates') % If you start by scanning A and B gates together, select B values and then remove them from the scan. 
    pt = ginput(1);
    if strcmp(scandata.name,'right')
        smset('3b',pt(1));
        smset('4b',pt(2));
        scandata.CB.loops(1).setchan = {'3a'};
        scandata.CB.loops(2).setchan = {'4a'};
        fprintf('Set 3b to %3.3f and 4b to %3.3f \n',pt(1),pt(2));
    else
        smset('2b',pt(1));
        smset('1b',pt(2));
        scandata.CB.loops(1).setchan = {'2a'};
        scandata.CB.loops(2).setchan = {'1a'};
        fprintf('Set 2b to %3.3f and 1b to %3.3f \n',pt(1),pt(2));
    end
end
if isopt(opts,'set') % Set scanned gates to clicked value 
    % Edited to make it easier to click on charge scan data. 
    %figure(1000); 
    a = gca;
    xchan = a.XLabel.String;
    ychan = a.YLabel.String;
    if contains(xchan,',')
        com=strfind(xchan,',');
        xchan = xchan(1:com(1)-1); 
    end
    if contains(ychan,',')
        com=strfind(ychan,',');
        ychan = ychan(1:com(1)-1); 
    end
    newvals = ginput(1);
    smset(xchan,newvals(1));
    smset(ychan,newvals(2));
    data = newvals; 
end
if isopt(opts,'ldrng') % Load and copy range from file to current scan. 
    file = uigetfile('');
    d = load(file);
    scandata.(scanname).loops(1).rng = d.scan.loops(1).rng;
    scandata.(scanname).loops(2).rng = d.scan.loops(2).rng;
    fprintf('Copying ranges from file %s \n', file(1:end-4));
end
if isopt(opts,'load') % Load and copy range, trafofn, npoints from file to current scan
    file = uigetfile('');
    d = load(file);
    scandata.(scanname).loops(1).rng = d.scan.loops(1).rng;
    scandata.(scanname).loops(2).trafofn = d.scan.loops(2).trafofn;
    scandata.(scanname).loops(2).rng = d.scan.loops(2).rng;
    scandata.(scanname).loops(1).npoints = d.scan.loops(1).npoints;
    scandata.(scanname).loops(2).npoints = d.scan.loops(2).npoints;
    fprintf('Copying ranges, points, trafofn from file %s \n', file(1:end-4));
end
if isfield(config,'copy') % Copy range and trafofn from another scan. 
    scandata.(scanname).loops(2).rng = scandata.(config.copy).loops(2).rng;
    scandata.(scanname).loops(1).rng = scandata.(config.copy).loops(1).rng;
    scandata.(scanname).loops(1).trafofn = scandata.(config.copy).loops(1).trafofn;
    scandata.(scanname).loops(1).setchan = scandata.(config.copy).loops(1).setchan;
    fprintf('Copying ranges, trafofn, and loop1 setchan from %s to %s \n', config.copy, scanname);
end
if isopt(opts,'csens') % Copy sensor sensGate range to SD scan, 
    if strcmp(scandata.activeScan,'sensLock')
        scandata.SDLock.loops(2).rng = scandata.sensorLock.loops(2).rng;
    elseif strcmp(scandata.activeScan,'sens')
        scandata.SD.loops(2).rng = scandata.sensor.loops(2).rng;
    end
end
if isopt(opts,'csd') % copy SD scan range to sensor scan.
    if strcmp(scandata.activeScan,'sensLock')
        scandata.sensorLock.loops(2).rng = scandata.SDLock.loops(2).rng;
    elseif strcmp(scandata.activeScan,'sens')
        scandata.sensor.loops(2).rng = scandata.SD.loops(2).rng;
    end
end
%% Intelligently adjusting scans
if isopt(opts,'hyst')
    autoscan('CB','flip'); data=autoscan('CB');
    good = checkNan(data);  if ~good, return; end
    data = autoscan('CB','rot');
    good = checkNan(data);  if ~good, return; end
    autoscan('CB','flip');     data = autoscan('CB','rot');
end % Run 3 scans with all possible flips.
if isopt(opts,'gca') % Zoom in on scan. Increase ramprate to max.
    scandata.(scanname) = smscanpar(scandata.(scanname),'gca');
end
if isopt(opts,'flip') % Flip scan in horizontal dir.
    scandata.(scanname)=smflip(scandata.(scanname),1);
end
if isopt(opts,'fast') % Run mini scan that should go faster
    scantmp = scandata.(scanname);
    if any(strcmp(scanname,lockscans))
        autoscan(scanname,'',struct('ramptime',-0.0625));
    end
    autoscan(scanname,'',struct('ramprate',smdata.channels(chl(scandata.(scanname).loops(1).setchan{1})).rangeramp(3),'npoints2',8));
    autoscan(scanname)
    scandata.(scanname)=scantmp;
    return;
end
if isopt(opts,'crs') % Coarse scan, 12 mV spacing in X direction.
    if any(strcmp(scanname,lockscans))
        autoscan(scanname,'',struct('ramptime',-0.0625));
    end
    setChan = chl(scandata.(scanname).loops(1).setchan{1});
    autoscan(scanname,'',struct('ramprate',smdata.channels(setChan).rangeramp(3),'pointSpacing1',12e-3));
    autoscan(scanname)
end
if isopt(opts,'fn') % Fine scan, 1 mV spacing in X and 3 mV in Y.
    autoscan(scanname,'',struct('pointSpacing1',1e-3,'pointSpacing2',3e-3));
    autoscan(scanname)
end
if isopt(opts,'rot') % Rotate scan 90 degrees and run.
    scantmp = scandata.(scanname);
    autoscan(scanname,'',struct('ramprate',70e-3,'rng1',scantmp.loops(2).rng,'rng2',scantmp.loops(1).rng));
    scandata.(scanname).loops(1).setchan = scantmp.loops(2).setchan;
    scandata.(scanname).loops(2).setchan = scantmp.loops(1).setchan;
    scandata.(scanname).data.opts = 'rotate';
    data = autoscan(scanname);
    scandata.(scanname)=scantmp;
    return;
end
if isopt(opts,'up') % Frame up
    rng = scandata.(scanname).loops(2).rng;
    sortRng = sort(rng);
    scandata.(scanname).loops(2).rng = sortRng;
    autoscan(scanname);
end
if isopt(opts,'down') % Frame down 
    rng = scandata.(scanname).loops(2).rng;
    sortRng = sort(rng);
    scandata.(scanname).loops(2).rng = [sortRng(2) sortRng(1)];
    autoscan(scanname);
end
if isopt(opts,'finer1') % Double npoints 
    scandata.(scanname).loops(1).npoints = 2* scandata.(scanname).loops(1).npoints;
elseif isopt(opts,'coarser1') % Halve npoints
    scandata.(scanname).loops(1).npoints = 0.5* scandata.(scanname).loops(1).npoints;
end
if isopt(opts,'finer2') % Double npoints
    scandata.(scanname).loops(2).npoints = 2* scandata.(scanname).loops(2).npoints;
elseif isopt(opts,'coarser2') % Halve npoints
    scandata.(scanname).loops(2).npoints = 0.5* scandata.(scanname).loops(2).npoints;
end
if isfield(config,'nudge')
    %FIX ME
    traf = scandata.(scanname).loops(2).trafofn(2);
    if contains(func2str(traf.fn),'interp1')
        if strcmp(config.nudge,'pos')
            traf.args{2} = traf.args{2}+0.005;
        else
            traf.args{2} = traf.args{2}-0.005;
        end
    else
        if strcmp(config.nudge,'pos')
            %scandata.(scanname).loops(2).trafofn(2).args{1}(3)=scandata.(scanname).loops(2).trafofn(2).args{1}(3)+0.01;
            for i = 1:length(scandata.(scanname).loops(2).trafofn(2).args)
                scandata.(scanname).loops(2).trafofn(2).args{i}(1) = scandata.(scanname).loops(2).trafofn(2).args{i}(1)+0.005;
            end
        else
            %scandata.(scanname).loops(2).trafofn(2).args{2}=scandata.(scanname).loops(2).trafofn(2).args{2}-0.01;
            scandata.(scanname).loops(2).trafofn(2).args{1}(3)=scandata.(scanname).loops(2).trafofn(2).args{1}(3)-0.005;
        end
    end
    scandata.(scanname).loops(2).trafofn(2) = traf; 
end
if isopt(opts,'gofast') % change ramprate to maximum.
    autoscan(scanname,struct('ramprate',smdata.channels(chl(scandata.(scanname).loops(1).setchan{1})).rangeramp(3)));
end
%% Auxiliary (inst control, print info, etc. )
if isopt(opts,'sweep')
    setvals = linspace(config.srng(1),config.srng(2),config.snpoints);
    for i = 1:config.snpoints
        smset(config.ssetchan, setvals(i));
        autoscan(scanname);
    end
end
if isopt(opts,'close') % close off parallel channel, save value to scandata.
    closeGates = [scandata.sets(1).config.closeGates, scandata.sets(2).config.closeGates];
    lockin = findChans(scandata.CB.loops(2).getchan,1);
    scandata.config.closeVal = closeParallel(closeGates, lockin, 1e-10);
    scandata.sets(1).config.closeVal = scandata.config.closeVal;
    scandata.sets(2).config.closeVal = scandata.config.closeVal;
end
if isopt(opts, 'lockin')
    smaSRS830('LockinB',struct('opts','def')); smaSRS830('LockinB',struct('meas','a','line','all'));
    smset('LockinTauB',0.03);
    smset('LockExcB',0.010);
    smset('LockSensB',.1);
    smrange({'LockinB','QPCbufB'},'div',-1e8);
    fprintf('Initialize lock in for sensing. Setting Tau, Exc, Sens, divider, filtering \n');
end
if isfield(config,'switch') % If Jim's switch is connected, use this to change it's value.
    switch config.switch
        case 'CB'
            smaSwitch({'OL1','VA'; 'OL2','IA'; 'OS1top','open';'OS1bot','open'})
        case 'sens'
            smaSwitch({'OL1','open'; 'OL2','open'; 'OS1top','VB';'OS1bot','IB'})
        case 'both'
            smaSwitch({'OL1','VA'; 'OL2','IA'; 'OS1top','VB';'OS1bot','IB'})
        case 'ground'
            smaSwitch({'OL1','ground'; 'OL2','ground'; 'OS1top','ground';'OS1bot','ground'})
    end
end
if isopt(opts,'print') %Print info about scan
    niceprint(scandata.(scanname));     
    return
end
sleep
end

function good= checkNan(data)
if any(isnan(data{1}))
    action = lower(questdlg('go to next scan?','question','yes','no','yes'));
    switch action
        case 'no'
            fprintf('Quitting')
            good = 0;
        case 'yes'
            good = 1;
        otherwise
            good = 1;
    end
else
    good = 1;
end
end