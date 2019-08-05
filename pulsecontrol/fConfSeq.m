function scanSeq = fConfSeq(plsGrpList, conf)
% Configure scan with pulsegroups. 
% function scanSeq = fConfSeq(plsgrp, conf)
% allows multi length pulses. 
% Generate a scan using pulse groups plsgrp.
% conf can have any or none of:
%  nloop: default 200 number of ties to repeat experiment without interuption
%  oversamp: default 1 (number of measured points per pulse. main time to get more is t1 scan, when we want every point. 
%  nrep:    default 50.
%        number of times to repeat experiment (with fb between reps). (total is nrep*nloop)
%  datachan: default {'DAQ1','DAQ2'}; % which channel DAQ used (can be both)
%  auxchan: default {'Time'}; % by default we measure time at each rep. 
%  fastmode: default 1 % multiple pulses on single line, 2 is similar, except we keep 
%           first loop which does all the pulses in one point in first loop.
%  setmask: default oversamp < 2
%  snglshot: default 2:   0 - no snglshot  1 - smasnglshot   2 - built-in snglshot
%  hwsampler: hardware sample rate, for DAQ 40 Mhz default, nan to leave alone.
%  npulse: default nan usually defined by number of pulses in group. 
%          If given, will repeat pulse npulse times. 
%  sampler: default nan
%  extclk: default based on DAQ inst. 
%      valid options are: pol, fb, nodisp, nosave, ampok, nocheck,raw, nosave
%            raw: save raw (non-averaged) data as well.
%            pol: add polarization code
%             fb: add feedback control code
%         nodisp: disable display
%          ampok: skip check that amplifier is off on awg
%        nocheck: skip sanity checks on groups
%         nosave: save only at end
%           swfb: enable software feedback
%       ramseyfb: enable software ramsey feedback
%          qsave: save qdata with the scan
%% Set up defaults
global awgdata; global smdata; global plsdata; global fbdata; 
if ~exist('conf','var'), conf=struct(); end
if iscell(conf), conf=struct(conf{:}); end
daqChan=smchanlookup('DAQ1');
daqInst=smdata.channels(daqChan).instchan(1);
conf = def(conf,'opts','');  opts=conf.opts;
conf = def(conf,'nloop',200); nloop = conf.nloop;
conf = def(conf,'oversamp',1);  oversamp = conf.oversamp;
conf = def(conf,'nrep',50);  nrep = conf.nrep;
conf = def(conf,'sampler',nan);  sampler = conf.sampler;
conf = def(conf,'npulse',nan);  npulse=conf.npulse;
conf = def(conf,'datachan',{'DAQ1','DAQ2'}); 
if strcmp(conf.datachan,'both'), conf.datachan = {'DAQ1','DAQ2'}; end
datachan = conf.datachan;
conf = def(conf,'auxchan',{'Time'}); auxchan = conf.auxchan;
conf = def(conf,'fastmode',1);  fastmode=conf.fastmode;
conf = def(conf,'setmask', oversamp < 2 ); setmask=conf.setmask;
conf = def(conf,'snglshot',2); snglshot = conf.snglshot;
conf = def(conf,'hwsampler',40e6);  hwsampler = conf.hwsampler;
conf = def(conf,'extclock',smdata.inst(daqInst).data.extclk); extclk=conf.extclock;

if ischar(datachan), datachan = {datachan}; end
if ischar(auxchan), auxchan = {auxchan}; end
if isnan(oversamp) || oversamp==0, oversamp = 1; end
if feature('IsDebugMode'), warning('debugger is on, fool!'); end
scanSeq.cleanupfn=[]; scanSeq.data.conf=conf;  % Handy documentation.
%% check AWG groups up to date, configure DAQ and readout. 
plsGrpList = awggrpind(plsGrpList); % List of groups to run. 
if plsGrpList(1) > 0 % Check if pulse info needs to be updated, set up sampler, npulse.
    if isnan(npulse)
        npulse = awgdata(1).pulsegroups(plsGrpList(1)).npulse(1); % number of pulses in that group. 
    elseif npulse < 0
        npulse = awgdata(1).pulsegroups(plsGrpList(1)).npulse(1) * abs(npulse); % So if negative number, we repeat each pulse npulse x.
    end
    if ~isopt(opts,'nocheck') % Check that npulse correct and group is updated. 
        for i = 1:length(plsGrpList)
            if npulse ~= awgdata(1).pulsegroups(plsGrpList(i)).npulse(1)
                warning('fConfSeq:PulseNum','Pulse number mismatch; length of group %s (%d) is not %d\n',...
                    awgdata(1).pulsegroups(plsGrpList(i)).name,awgdata(1).pulsegroups(plsGrpList(i)).npulse(1),npulse);
            end
            if awgdata(1).pulsegroups(plsGrpList(i)).changed
                error('fConfSeq:StalePulse','Group needs uploading\n');
            end
        end
    end
    if isnan(sampler)
        plsLength = awgdata(1).pulsegroups(plsGrpList(1)).zerolen;
        maxplsTime = abs(max(abs(plsLength(:, 1)))) * max(1, awgdata(1).pulsegroups(plsGrpList(1)).nrep(1));
        sampler = awgdata(1).clk/maxplsTime; % inverse pulse time.
        for i=2:length(plsGrpList)
            plsLength = awgdata(1).pulsegroups(plsGrpList(i)).zerolen;
            plsTimeMean = abs(mean(plsLength(:, 1))) * max(1, awgdata(1).pulsegroups(plsGrpList(i)).nrep(1));
            samplerMean = awgdata(1).clk/plsTimeMean;
            if sampler ~= samplerMean
                warning('fConfSeq:PulseLen','Pulse length mismatch; sample rate for group %s is %g, not %g\n',awgdata(1).pulsegroups(plsGrpList(i)).name,samplerMean,sampler);
            end
        end
    end
    if (isfield(awgdata(1).pulsegroups(plsGrpList(1)), 'jump') && ~isempty(awgdata(1).pulsegroups(plsGrpList(1)).jump)) || ...
            awgdata(1).pulsegroups(plsGrpList(1)).nrep(1) == 0 || ...  % single pulse repeated indefinitely -> same logic
            awgdata(1).pulsegroups(plsGrpList(1)).nrep(1) == 1  % single pulse repeated indefinitely -> same logic
        npulse = npulse*nloop;
        fastmode = 1;
    else
        nloop = 0;
    end
    seqind = [awgdata(1).pulsegroups(plsGrpList).seqind];
    scanSeq.data.pulsegroups = awgdata(1).pulsegroups(plsGrpList);
else % Single pulse, preceded by trigger.
    seqind = awgseqind(abs(plsGrpList))-1; % subtract pulse to jump to preceding trigp.
    if isnan(npulse)
        npulse=1;
    end
    if isnan(sampler)
        pd = awgdata(1).pulsedata(abs(plsGrpList(1)));
        sampler = pd.clk/(pd.pulsetab(1, end) * pd.tbase);
    end
    if nloop > 0 % treat same way as jump sequence
        npulse = npulse*nloop;
        fastmode = 1;
    end    
end

% Start by creating first loop containing npoints = number of pulses in each line,with ramptime = pulse time.
scanSeq.loops(1).npoints = abs(npulse*oversamp);
scanSeq.loops(1).rng = [];
scanSeq.configfn.fn = @smabufconfig2; 
scanSeq.configfn.args = {'arm pls', 1}; % This configures only one channel. Seems to work.
scanSeq.loops(1).ramptime = 1/(oversamp * sampler); %  pulse time / number of points per pulse. 
scanSeq.loops(2).getchan = [datachan, auxchan];
scanSeq.loops(2).setchan = [];
% Set pulseline at the beginning of each line. 
scanSeq.loops(2).prefn(2).fn = '@(x, seqind, loop) smset(''PulseLine'', seqind(mod(floor(x(loop)-1), end)+1),'''',''quiet'')'; 
scanSeq.loops(2).prefn(2).args = {seqind, 2};

if length(plsGrpList) == 1  % Second loop is repeating single group
    scanSeq.loops(2).npoints = nrep;
else % Second loop goes through grp list, third loop does reps. 
    scanSeq.loops(2).npoints = length(seqind); 
    scanSeq.loops(3).npoints = nrep;
    scanSeq.loops(3).setchan = 'count';
end
if extclk % Use RFfreq3 as clock. Otherwise, set the samprate before configuring DAQ 
    if ~isnan(hwsampler)
        smset('RFfreq3',hwsampler);
        scanSeq.consts(1).setchan='RFfreq3';
        scanSeq.consts(1).val=hwsampler;
        scanSeq.consts(2).setchan='samprate';
        scanSeq.consts(2).val=hwsampler;
    end
    smdata.inst(daqInst).data.extclk=1;
else
    if ~isnan(hwsampler)
        smset('samprate',hwsampler);  % Yuck.
        scanSeq.consts(1).setchan='samprate';
        scanSeq.consts(1).val=hwsampler;
    end
end
scanSeq = scanSeq.configfn.fn(scanSeq, scanSeq.configfn.args{:}); % This will configure desk basedo 
if scanSeq.loops(1).npoints ~= abs(npulse*oversamp) % checks for rounding error of machine (daq)
    error('Invalid record length.');        
end
scanSeq.disp=[];
for i=1:length(datachan) % Set up disp
    scanSeq.disp(end+1).loop = 2;
    scanSeq.disp(end).channel = i;
    scanSeq.disp(end).dim = 1;
    scanSeq.disp(end).mod = 2;
    scanSeq.disp(end+1).loop = 2;
    scanSeq.disp(end).channel = i;
    scanSeq.disp(end).dim = 2;
    scanSeq.disp(end).mod = 2;
end 
if npulse*oversamp/nloop == 1 % make disp work in degenerate case
    scanSeq.disp(2) = []; 
end
switch fastmode % if 1, delete the first loop.
    case 1 % remove inner loop for fast acq.
        scanSeq.configfn.args{3} = [scanSeq.loops(1).npoints, 1/abs(scanSeq.loops(1).ramptime)];        
        if length(datachan)>1
            scanSeq.configfn.args{1} = 'fast arm pls chans';
        else
            scanSeq.configfn.args{1} = 'fast arm pls';
        end
        scanSeq.loops(1) = [];
        scanSeq.loops(1).prefn(2).args{2} = 1; % loop var for setting pulseline. 
        scanSeq.loops(end).setchan = 'count';
        if smdata.inst(daqInst).data.nBuffers == 0 % we use async. reads now; no polling.  interrupt driven.
            scanSeq.loops(1).waittime = -1;
        end
        [scanSeq.disp.loop] = deal(1);
    case 2 % shorten inner loop instead.
        scanSeq.configfn.args{3} = [scanSeq.loops(1).npoints, 1/abs(scanSeq.loops(1).ramptime)];
        scanSeq.configfn.args{1} = 'fast arm pls';
        scanSeq.configfn.args{4} = 2;
        scanSeq.loops(1).ramptime = scanSeq.loops(1).ramptime * scanSeq.loops(1).npoints;
        scanSeq.loops(1).npoints = 1;
        scanSeq.loops(end).setchan = 'count';
        scanSeq.saveloop = [length(scanSeq.loops) 10];
end
if nloop % Create procfn for reshaping and averaging data. 
    npulse = npulse/nloop; % restore previous value
    scanSeq.loops(1).procfn(1).fn(1).fn = @reshape;
    scanSeq.loops(1).procfn(1).fn.args = {[npulse * oversamp, nloop]};
    scanSeq.loops(1).procfn(1).fn(2).fn = @nanmean;
    scanSeq.loops(1).procfn(1).fn(2).args = {2};
    scanSeq.loops(1).procfn(1).dim = npulse * oversamp;
    scanSeq.loops(1).procfn(2:length(datachan)) = scanSeq.loops(1).procfn(1);
    scanSeq.loops(1).procfn(length(datachan)+(1:length(auxchan))).fn = [];
end
%% Set up the mask and procfn. 
if setmask && oversamp == 1 % Set up the mask. 
    samprate = cell2mat(smget('samprate'));    
    for i = 1:length(plsGrpList)
        rdout(i) = awgdata(1).pulsegroups(plsGrpList(i)).readout;
        masktmp = [];
        for j = 1:size(rdout(i).readout,3) %this is across pulses in the group. 
            if isfield(rdout(i),'plens')
                nreps = max(1, awgdata(1).pulsegroups(plsGrpList(1)).nrep(1));
                lineTime = abs(rdout(i).plens(j)) * nreps/awgdata(1).clk; 
            else
                lineTime = 1/sampler;
            end
            maskVec = false(2,round(samprate*lineTime)); 
            for k = 1:2 % channels                                 
                chan=find(rdout(i).readout(:,1) == k);
                if isempty(chan)
                    chan = find(rdout(i).readout(:,1) == 0);
                end
                if isempty(chan)
                    chan=1;
                end
                % The proliferation of rounds below is needed to avoid machine rounding error.
                %maskInd are the indices in terms of sampling intervals on the DAQ where the card should acquire.
                maskInd=round(round(samprate * cumsum(rdout(i).readout(chan,2:3,j)))*round(plsdata.tbase*1e-3)*1e-6); 
                maskVec(k,maskInd(1):maskInd(2))=true;
            end
            masktmp = [masktmp, repmat(maskVec,1,rdout(i).reps(j))]; % reps is always one?
            if j==1 % simple error checking
                firstRdout=diff(maskInd);
            else
                if diff(maskInd)~=firstRdout
                    error('Mask size not consistent');
                end
            end
        end
        masks{i} = (masktmp>0); %#ok<*AGROW> %need to make it a logical
    end
    scanSeq.loops(1).prefn(3:end+1)=scanSeq.loops(1).prefn(2:end);
    cntrlfn_str= func2str(smdata.inst(daqInst).cntrlfn);
    scanSeq.loops(1).prefn(2).fn=['@(x,masks,loop) ',cntrlfn_str,'([smchaninst(''' datachan{1} '''),6],masks{mod(floor(x(loop)-1),end)+1})'];
    scanSeq.loops(1).prefn(2).args = {masks,1};
    if length(unique([rdout.plens]))>1 % pulses changing lengths!        
        scanSeq.loops(1).prefn(3:end+1)=scanSeq.loops(1).prefn(2:end);
        % next line will set the numPls in smdata.inst(14).data so the daq driver knows what to do
        scanSeq.loops(1).prefn(2).fn = sprintf('@(x,t) %s([%d,%d,%d],t)',cntrlfn_str,daqInst,7,1); % 7 is this field
        tmp = vertcat(awgdata(1).pulsegroups(plsGrpList).npulse);
        scanSeq.loops(1).prefn(2).args = {tmp(:,1)};
    end
    scanSeq.cleanupfn.fn = @smaconfigwrap;
    scanSeq.cleanupfn.args = {smdata.inst(daqInst).cntrlfn, [smchaninst(datachan{1}), 6], []};
end
scanSeq.cleanupfn(end+1).fn = @smaconfigwrap; %now set the n_pls_in_grp flag to empty
scanSeq.cleanupfn(end).args = {smdata.inst(daqInst).cntrlfn, [daqInst, 7, 1], []}; 
if isopt(opts,'raw') % Set up extra channel that doesn't average data. 
    cpfn.fn = []; %copy function
    cpfn.args = {};
    scanSeq.loops(1).procfn(1).fn = [cpfn, scanSeq.loops(1).procfn(1).fn];
    scanSeq.loops(1).procfn(1).fn(1).inchan = 1:length(datachan);
    scanSeq.loops(1).procfn(1).fn(1).outchan = length(scanSeq.loops(1).procfn) + (1:length(datachan));
    for k=scanSeq.loops(1).procfn(1).fn(1).outchan
        scanSeq.loops(1).procfn(k).fn=[];
        scanSeq.loops(1).procfn(k).dim=nloop*npulse * oversamp;
    end
end
% time convention: time specifies end of sampling interval, i.e. time = (1:n)*dt,% consistent with loadpulse.m. Want samples ending
scanSeq.data.snglshot = snglshot;  % Later processing wants this.
switch snglshot
    case 1 % use single shot framework to gather histogram and T1 info.
        clear pp;
        ndc = length(datachan);
        [pp.datadef(1:ndc).type] = deal('mom'); % std readout
        [pp.datadef(ndc + (1:ndc)).type] = deal('ghist'); %histograms
        for i = 1:ndc
            pp.datadef(i).par = {(1:npulse) + (i-1)*npulse};
            chan = smchaninst(datachan{i});
            pp.datadef(i+ndc).par = {linspace(-6e-3, 16e-3, 45) + fbdata.refval(chan(2)), (1:npulse) + (i-1)*npulse};
        end
        pp.datadef(2*ndc+1).type = 'ave'; % T1 measurement
        % this makes some of the above code obsolete, e.g. setting sample count in configfn and  processing functions. Likely requires fastmode.
        scanSeq = smaSnglShot2(scanSeq, pp, 'raw dec', [], nloop);
    case 2 % Configure histogram using fbdata. 
        scanSeq.loops(1).procfn(1).fn(2:end+1) = scanSeq.loops(1).procfn(1).fn(1:end);
        scanSeq.loops(1).procfn(1).fn(1).fn = [];
        scanSeq.loops(1).procfn(1).fn(1).args = {};
        scanSeq.loops(1).procfn(1).fn(1).inchan = 1:length(datachan);
        scanSeq.loops(1).procfn(1).fn(1).outchan = length(scanSeq.loops(1).procfn) + (1:length(datachan));
        if oversamp > 1
            for i = 1:length(datachan)
                npts=1000;
                scanSeq.loops(1).procfn(end+1).dim = [npts oversamp];
                scanSeq.loops(1).procfn(end).fn(end+1).fn = @reshape;
                scanSeq.loops(1).procfn(end).fn(end).args = {[oversamp nloop]};
                scanSeq.loops(1).procfn(end).fn(end+1).fn = @transpose;
                scanSeq.loops(1).procfn(end).fn(end).args = {};
                scanSeq.loops(1).procfn(end).fn(end+1).fn = @histc;
                chan = smchaninst(datachan{i});
                %scanSeq.loops(1).procfn(end).fn(end).args = {linspace(-20e-3, 20e-3, npts) + fbdata.refval(chan(2))};
                scanSeq.loops(1).procfn(end).fn(end).args = {linspace(-120e-3, 120e-3, npts)};
            end
        else
            for i = 1:length(datachan)
                npts=1000;
                scanSeq.loops(1).procfn(end+1).dim = npts;
                scanSeq.loops(1).procfn(end).fn(1).fn = @histc;
                chan = smchaninst(datachan{i});
                %scanSeq.loops(1).procfn(end).fn(1).args = {linspace(-20e-3, 20e-3, npts) + fbdata.refval(chan(2))};
                scanSeq.loops(1).procfn(end).fn(end).args = {linspace(-120e-3, 120e-3, npts)};                
            end
        end
end
%% Checks, feedback, cleanupfns. 
scanSeq.loops(1).stream = 1; 
if ~isopt(opts,'nosave')
    scanSeq.saveloop = [length(scanSeq.loops), 10];
else
    scanSeq.saveloop = [-1, 1];
end
if ~isopt(opts,'ampok') %check for channels in amp mode
    israw = awgcntrl('israw');
    for i = 1:length(israw)
        if israw(i) ~=1, fprintf('channel %d is in amp mode\n', i); end
    end
end
if ~contains(opts,'offok') %check if AWGs on
    ison = awgcntrl('ison');
    if any(~ison)
        disp('AWG is off.  <a href="matlab:awgcntrl(''on start wait err'');">Turn it on?</a> <a href="matlab:disp(''exiting...'');dbquit">Exit?</a> ')
        keyboard
    end
end
if isopt(opts,'pol')
    scanSeq = fScanPol(scanSeq);
elseif isopt(opts,'fb')
    [~,scanSeq] = fScanPol(scanSeq);
    fbdata.x=plsinfo('xval',plsGrpList(1));
    fbdata.pulses=1:length(fbdata.x);
    scanSeq.configfn(end+1).fn=@smaconfigwrap;
    scanSeq.configfn(end).args{1}=@feedback_control;
    figure(1034);
end
if isopt(opts,'swfb'), scanSeq=swfbConfig(scanSeq); end
if isopt(opts,'ramseyfb'), scanSeq=swfbConfig2(scanSeq); end
if isopt(opts,'multifb'), scanSeq=swfbConfigMulti_v2(scanSeq); end
if isopt(opts,'qsave')
    scanSeq.configfn(2:end+1)=scanSeq.configfn;
    scanSeq.configfn(1).fn=@saveQubitParams;
    scanSeq.configfn(1).args={};
end
if isopt(opts,'predbz') || isopt(opts,'logdbz')
    a= struct('datachan',datachan);
    cc=struct('args',a,'name','pre_dbz','n_rep',1);
    if isfield(scanSeq,'configfn') && ~isempty(scanSeq.configfn) %want this to run before configuring daq card
        scanSeq.configfn(2:end+1) = scanSeq.configfn;
        scanSeq.configfn(1).fn = @sma_dBz_save;
        scanSeq.configfn(1).args = {cc};
    else %should really never happen
        scanSeq.configfn.fn = @sma_dBz_save;
        scanSeq.configfn.args = {cc};
    end
end
if isopt(opts,'postdbz') || isopt(opts,'logdbz')
    a= struct('datachan',datachan);
    cc=struct('args',a,'name','post_dbz','n_rep',1);
    scanSeq.cleanupfn(end+1).fn = @sma_dBz_save;
    scanSeq.cleanupfn(end).args = {cc};
end
if isopt(opts,'FPGA') || isopt(opts,'fpga')
    scanSeq.cleanupfn(end+1).fn = @get_FPGA_freqs;
    scanSeq.cleanupfn(end).args = {};
    if isopt(opts,'no_params')
        scanSeq.configfn(end+1).fn = @store_FPGA_params;
        scanSeq.configfn(end).args = {};
    end
    scanSeq.loops(1).postfn(end+1).fn = @fpga_control_output;
    scanSeq.loops(1).postfn(end).args = {};
end
if isopt(opts,'nodisp'), scanSeq.disp=[]; end
end