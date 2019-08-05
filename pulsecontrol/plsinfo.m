function val = plsinfo(ctrl, group)
% Return information on pulse. 
% val = plsinfo(ctrl, group, ind, time)
% ctrl: xval, zl, gd, log(=pl), sl, ls, rev, stale, params, ro
% ls: List saved pulsegroups.  group can be an optional mask like 'dBz*'
% stale: return 1 if group not uploaded, or has been updated since uploading. 

global plsdata; global awgdata;
if exist('group','var') && ~ischar(group)
    group = awgdata(1).pulsegroups(group).name;
end

switch ctrl
    case 'stale' % This doesn't work correctly for sequence combined groups
        ind=awggrpind(group);
        if isnan(ind)
            val = true;
        elseif awgdata(1).pulsegroups(ind).changed == true
            val = true;
        else
            val = false;
        end
    case 'ls'
        currFolder = pwd;
        cd(plsdata.grpdir);
        if ~exist('group','var')
            group = '*';
        end
        val=dir(['pg_', group, '.mat']);
        val = {val.name};
        val=regexprep(val,'^pg_','');
        val=regexprep(val,'.mat$','');
        if nargin > 1
            for i=1:length(val)
                fprintf('%s\n',val{i});
            end
            clear val;
        end
        cd(currFolder);
end
end