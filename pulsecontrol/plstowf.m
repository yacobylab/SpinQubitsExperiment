function pulse = plstowf(pulse, dict)
%Convert any (valid) pulse format to wf format. (int, tab, wf, elem)
% pulse = plstowf(pulse, dict)
% If pulse is an integer, it is taken as a database index. 
% Otherwise, pulse is a struct with fields:
% format: 'tab', 'wf', 'elem
% data: struct with following fields
%   format = 'tab':
% 	pulsetab
% 	marktab
% 	pulsefn
% 	readout
%   format = 'elem': val, time (array)
%   format = 'wf': wf, marker, readout
%
% name 
% taurc (optional)
% xval
%
% MARKTAB:
%  each column gives a marker pulse
%  [ start_time; ch1_mk1_wid ; ch1_mk2_wid ; ch_2_mk1_wid ...]
%  ie,
%  [ 2 ; 0 ; 1 ; 0 ; 1 ] fires markers ch1mk1 and ch2mk2 for 1 us starting at 2us.

global plsdata; global awgdata;
dt=1e-11; 
pulse = plsdefault(pulse);
if strcmp(pulse.format, 'wf') %alreay in final format
    return
elseif strcmp(pulse.format, 'elem')
    if exist('dict','var') && ~isempty(dict)
        pulse = pdapply(dict,pulse);
    end
    pulse = plstotab(pulse);
elseif ~strcmp(pulse.format, 'tab')
    error('Invalid format %s.', pulse.format);
end
pulseInfo = pulse.data;
pulse.data = [];
if ~isfield(pulseInfo, 'marktab'),     pulseInfo.marktab = []; end
if ~isfield(pulseInfo, 'pulsefn') 
    pulseInfo.pulsefn = [];
elseif ~isempty(pulseInfo.pulsefn) && ~isfield(pulseInfo.pulsefn, 'args')
    [pulseInfo.pulsefn.args] = deal(cell(2, 0));
end
if ~isfield(pulseInfo, 'readout'),    pulseInfo.readout = []; end
clk = unique([awgdata.clk]);
for c=1:length(clk)
    pulsetab = pulseInfo.pulsetab;
    nchan = size(pulsetab, 1)-1;
    npoints = round(max(pulsetab(1, :)) * plsdata.tbase * clk(c)/1e9);
    data = zeros(nchan, npoints+1);
    time = linspace(pulsetab(1, 1), pulsetab(1, end), npoints+1);
    if pulse.taurc == Inf
        avg = zeros(nchan, 1);
    else
        avg = 0.5 * sum((pulsetab(2:end, 2:end) + pulsetab(2:end, 1:end-1)).*repmat((pulsetab(1, 2:end) - pulsetab(1, 1:end-1)), nchan, 1), 2)./(pulsetab(1, end) - pulsetab(1, 1));
        pulsetab(2:end, :) = pulsetab(2:end, :) - repmat(avg, 1, size(pulsetab, 2));
    end
    for i = 2:size(pulsetab, 2) % writes the pulse into data using lines to connect the corners defined in pulsetab
        timeMask = time >= pulsetab(1, i-1)-dt & time <= pulsetab(1, i)+dt;
        for j = 1:nchan
            data(j, timeMask) = ((-pulsetab(j+1, i-1) + pulsetab(j+1,i)) * time(timeMask) + pulsetab(j+1,i-1) * pulsetab(1, i) - pulsetab(j+1,i) * pulsetab(1, i-1))./(pulsetab(1, i) -  pulsetab(1, i-1));
        end
    end
    for i = 1:length(pulseInfo.pulsefn)
        timeMask = time > pulseInfo.pulsefn(i).t(1) & time <= pulseInfo.pulsefn(i).t(2);
        for j = 1:nchan
            data(j, timeMask) = pulseInfo.pulsefn(i).fn{j}(time(timeMask)-pulseInfo.pulsefn(i).t(1), pulseInfo.pulsefn(i).args{j, :}) - avg(j);
        end
    end   % lets pulses be defined with functions (eg. sin, cos) instead of just lines
    data(:, end) = [];
    if any(isfinite(pulse.taurc)) % calculates input voltage based on output voltage (different bc of bias T.
        if length(pulse.taurc) == 1
            vc = cumsum(data, 2) * (pulsetab(1, end) - pulsetab(1, 1))/(npoints * pulse.taurc);
        else
            vc = cumsum(data, 2) * (pulsetab(1, end) - pulsetab(1, 1))./repmat(npoints * pulse.taurc', 1, npoints) ;
        end
    else
        vc=0;
    end
    marker = zeros(nchan, npoints, 'uint8');
    marktab = pulseInfo.marktab; % extend marktab to be right dimensions
    marktab(end+1:2*nchan+1,:) = 0;
    for i = 1:size(pulseInfo.marktab, 2)
        calcOnlyOnce=time(1:end-1) >= marktab(1, i) - dt;
        for j = 1:nchan
            for k = 1:2
                timeMask = calcOnlyOnce & time(1:end-1) < marktab(1, i) + marktab(j*2+k-1, i)-2e-11;
                marker(j, timeMask) = bitor(marker(j, timeMask), k);
            end
        end
    end
    pulse.data(c).marker = marker;
    pulse.data(c).wf = data + vc;
    pulse.data(c).readout = pulseInfo.readout;
    if isfield(pulseInfo,'elem')
        pulse.data(c).elem=pulseInfo.elem;
    end
    pulse.data(c).clk = clk(c);
    pulse.data(c).pulsetab = pulsetab;
    pulse.data(c).marktab = marktab;
end
pulse.format = 'wf';
end