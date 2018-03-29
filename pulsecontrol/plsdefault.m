function pulse = plsdefault(pulse)
% set defaults and guess format if not given
% Return pulse in database if pulse is a number 
% pulse = plsdefault(pulse)

global plsdata;
if ~isstruct(pulse) % look up pulse in plsdata 
    if pulse > length(plsdata.pulses)
       error('Requested pulse %d, but only %d are defined.  Did you plssync?',pulse,length(plsdata.pulses)); 
    end
    pulse = plsdata.pulses(pulse);        
    return;
end
if ~isfield(pulse, 'xval'),      [pulse.xval] = deal([]); end
if ~isfield(pulse, 'taurc'),     [pulse.taurc] = deal(Inf); end
if ~isfield(pulse, 'name'),      [pulse.name] = deal(''); end
if ~isfield(pulse, 'pardef'),    [pulse.pardef] = deal([]); end
if ~isfield(pulse, 'trafofn'),    [pulse.trafofn] = deal([]); end
if ~isfield(pulse(1), 'format') || isempty(pulse(1).format) % only implemented for single pulse
    if isreal(pulse.data)
        pulse.format = 'ind';
    elseif isfield(pulse.data, 'type')
        pulse.format = 'elem';
    elseif isfield(pulse.data, 'wf')
        pulse.format = 'wf';
    elseif isfield(pulse.data, 'pulsetab') || isfield(pulse.data, 'marktab') || isfield(pulse.data, 'pulsefn')
        pulse.format = 'tab';
    else
        error('Invalid format.\n')
    end
end 
end