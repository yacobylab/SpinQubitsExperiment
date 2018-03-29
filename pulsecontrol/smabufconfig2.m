function scan = smabufconfig2(scan, ctrl, getchanInd, config, loop)
% scan = smabufconfig2(scan, cntrl, getchanInd,config,loop
% Configure buffered acquisition for fastest loop using drivers. Usually used as configfn. 
% Typically includes triggering, arming, and configuring. Flow is: 
% At beginning of scan, when smabufconfig called, call cntrlfn with op 5 to configure the readout. 
% At beginning of each outer loop point (loop with getchan), instrument armed. 
% On first point of each inner loop point, instrument triggered along w/ setchans that ramp. 
% ctrl:   trig : use smatrigfn for triggering
%         arm : use smatrigfn to arm insts in loops(2).prefn(1) using arg 4 of ctrl function
%         fast: Acquire buffered date in 1st loop. Hence, don't use rate and time of first loop for timing. Config gives [npts, rate, nrec(optional)]
%         end: when used with arm, arm in a new prefn in the readout loop.  Otherwise uses the first one. 
%         pls: 
% getchanInd: indices to loops(2).getchan that do buffered readout (and must be armed, triggered).
% config: indices of setchans in inner loop to trigger, unless 'fast' (see above).
% loop: loop to perform buffered readout on. Default readout loop is 2, unless 'fast' given, in which case loop is 1.
% Possible extensions (not implemented):
% - configure decimation (see smarampconfig for code)

global smdata;
if ~exist('ctrl','var')
    ctrl = '';
end
if isopt(ctrl, 'fast') % Set which loop is used for readout, which ramps to trigger.
    if ~exist('loop','var') || isempty(loop)
        loop = 1;
    end
else
    if ~exist('loop','var') || isempty(loop)
        loop = 2;
    end
    if loop == 1
        error('Need to use fast control if you want readout in first loop')
    end
    setchans = smchaninst(scan.loops(loop-1).setchan);
    if exist('config','var') && ~isempty(config)
        setchans = setchans(config, :);
    end
end
getchans = smchaninst(scan.loops(loop).getchan); % Only select getchans with index getchanInd.
if exist('getchanInd','var') && ~isempty(getchanInd) && getchanInd ~= 0
    getchans = getchans(getchanInd, :);
end
if isopt(ctrl, 'fast') % configure readout inst.
    for i = 1:size(getchans, 1)
        args = num2cell(config);
        if isopt(ctrl,'pls')
            [config(1), config(2)] = smdata.inst(getchans(i, 1)).cntrlfn([getchans(i, :), 5], args{:},scan.data.pulsegroups(1).npulse(1),'pls');
        else
            [config(1), config(2)] = smdata.inst(getchans(i, 1)).cntrlfn([getchans(i, :), 5], args{:});
        end
    end
else
    for i = 1:size(getchans, 1)
        if isopt(ctrl,'pls')
            [scan.loops(loop-1).npoints, rate] = smdata.inst(getchans(i, 1)).cntrlfn([getchans(i, :), 5], scan.loops(loop-1).npoints, 1/abs(scan.loops(loop-1).ramptime),scan.data.pulsegroups(1).npulse(1),'pls');
        else
            [scan.loops(loop-1).npoints, rate] = smdata.inst(getchans(i, 1)).cntrlfn([getchans(i, :), 5], scan.loops(loop-1).npoints, 1/abs(scan.loops(loop-1).ramptime));
        end
        scan.loops(loop-1).ramptime = sign(scan.loops(loop-1).ramptime)/abs(rate);
    end
    if strfind(ctrl, 'trig') % add trigger to scan
        scan.loops(loop-1).trigfn.fn = @smatrigfn;
        scan.loops(loop-1).trigfn.args = {[setchans; getchans]};
    end
end
if strfind(ctrl, 'arm') %arm inst before each readout. 
    if strfind(ctrl,'end') % bad hack, need to fix this
        scan.loops(loop).prefn(end).fn = @smatrigfn;
        scan.loops(loop).prefn(end).args = {getchans, 4};
    else
        scan.loops(loop).prefn(1).fn = @smatrigfn;
        scan.loops(loop).prefn(1).args = {getchans, 4};
    end
end
end