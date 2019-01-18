function sleep(opts)
% Set experiment in state where it can be left alone by turning off RF, saving structs, 
% backing up structs. 
% function sleep(opts)
% opts:
%   'fast': doesn't save smadata, so it takes a lot less time.

global smdata; global tuneData; global fbdata; global qdata; global scandata;
if ~exist('opts','var'), opts=''; end

if now-smdata.backup > 24 * 60 * 60
    backup = true;
    smdata.backup = now;
else
    backup = 0;
end % once a day save the back up.
try % Turn off pulses
    smset('PulseLine',awgseqind('all_off_LR'));
catch
    fprintf('Error setting pulseline\n');
end
if exist('scandata','var') && ~isempty(scandata)
    save(scandata.file,'scandata');
    if backup
        save([scandata.file '_backup'],'scandata');
    end
end
if exist('smdata','var') && ~isempty(smdata) && isempty(strfind(opts,'fast')) %#ok<STREMP>
    save(smdata.files.smdata,'smdata')
    if backup
        save([smdata.files.smdata '_backup'],'smdata')
    end
elseif ~isempty(strfind(opts,'fast')) %#ok<STREMP>
    fprintf('snoozing...');
else
    warning('IDIOT')
end
if exist('tuneData','var') && ~isempty(tuneData)
    save(tuneData.file,'tuneData');
    if backup
        save([tuneData.file '_backup'],'tuneData');
    end
else
    warning('Finely tuned IDIOT');
end
if exist('fbdata','var') && ~isempty(fbdata)
    save(fbdata.file,'fbdata');
    if backup
        save([fbdata.file '_backup'],'fbdata');
    end
else
    warning('feedbacky IDIOT');
end
try
    brick = inl('LabBrick');
    for i = 1:length(brick)
        smaLabBrick(brick(i),'save');
    end % Save lab brick states
catch
    warning('LabBricks not responding')
end
fprintf('zzzzzzz\n');
end