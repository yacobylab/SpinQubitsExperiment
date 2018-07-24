function rinds = plslist(rng, name)
% Print list of pulses
% inds = plslist(rng, name)
% rng: pulse indices to be shown. Default all if not given or empty.
% single number means last n pulses.
% name: show only pulses with this name.
% inds: indices of pulses displayed

global plsdata

if nargin < 1 || isempty(rng)
    rng = 1:length(plsdata.pulses);
end

if rng(1) <= -1
    rng = length(plsdata.pulses)+rng:length(plsdata.pulses);
end

pulses = plsdata.pulses(rng);

if nargin >= 2
    inds = strcmp(name, {pulses.name});
    rng = rng(inds);
    pulses = pulses(inds);
end;

if nargout >= 1
    rinds = rng;
else
    fprintf('%-6s  %-10s  %-10s  %-10s  %-10s\n', 'pulse', 'name', 'xdata', 'length', 'format')
    fprintf('--------------------------------------------------\n');
    for i = 1:length(rng)
        if(~isempty(pulses(i).format))
          switch pulses(i).format
              case 'wf'
                  len = size(pulses(i).data.wf, 2)/(plsdata.tbase);
              case 'tab'
                  len = max(pulses(i).data.pulsetab(1,:));
              otherwise
                  len = nan; % not implemented
          end
          if isempty(pulses(i).xval)
              pulses(i).xval = nan;
          end
          fprintf('%6d  %-10s  %10g  %10g  %-10s\n', rng(i), pulses(i).name, pulses(i).xval(1), len, pulses(i).format);
          
        end
    end
end