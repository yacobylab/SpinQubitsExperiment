function plsnum = plsreg(pulse, plsnum)
% Adds pulse to plsdata.pulses.   
% plsnum = plsreg(pulse, plsnum)
% The return value plsnum is the pulse index for the pulse.
% plsnum defaults to adding to the end.

global plsdata;
pulse = plsdefault(pulse);

if nargin < 2
    plsnum = length(plsdata.pulses) + 1;
end

% check format?
if plsnum < length(plsdata.pulses) && ~isequal(plsdata.pulses(plsnum),pulse)
   fprintf('pulse %i already exists.',plsnum);
    q=(input('Overwrite [y/n]','s'));
    if ~strcmpi(q,'y')
        return;
    end
end

if isfield(plsdata,'pulses') && ~isempty(plsdata.pulses)
    plsdata.pulses(plsnum) = orderfields(pulse, plsdata.pulses);
else
    plsdata.pulses(plsnum) = pulse;
end

