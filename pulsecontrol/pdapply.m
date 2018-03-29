function [pulse, changedout]=pdapply(pd,pulse,ind)
% Apply a pulse dictionary to a pulse. Return the new pulse. Changed is true if the application was non-trivial
% function [pulse,changedout]=pdapply(pd, pulse, ind)
% new feature; an entry like '@1;foo,2;bar,...' will expand to foo for the first pulse, 2 for the second, etc.
% The first dictionary in the list has priority. 
if ~exist('ind','var'), ind=''; end
changedout = 0; 
if ~strcmp(pulse.format,'elem')
elseif iscell(pd) % If pd is a cell array, apply each dictionary in sequence.  This allows for some *neat* effects. :p
    while 1
        changed = 0;
        for i=1:length(pd)
            [pulse,changed2] = pdapply(pd{i},pulse,ind);
            changed = changed || changed2;
        end
        changedout = changed | changedout;
        if ~changed
            break;
        end
    end
else
    if ischar(pd),   pd=pdload(pd); end
    changed = 1; startInd=1;
    while changed
        changed=0;
        for i=startInd:length(pulse.data)
            if ~changed && (pulse.data(i).type(1) == '#')
                entries=regexp(pulse.data(i).type(2:end),'(\d*;)?([^,])*','tokens');
                for e=1:length(entries)
                    if isempty(ind) || isempty(entries{e}{1}) || (ind == str2double(entries{e}{1}))
                        changed = 1;
                        pulse.data(i).type = entries{e}{2};
                    end
                end
            end
            if ~changed && pulse.data(i).type(1) == '@'
                elemName = pulse.data(i).type(2:end);
                if isfield(pd,elemName) % check that element is in dictionary 
                    elemInfo=pd.(elemName);
                    if ischar(elemInfo)
                        elemInfo={elemInfo};
                    end
                    if iscell(elemInfo)
                        elemInfo=struct('type',elemInfo,'time',[],'val',[]);
                    end
                    plsInfo=pulse.data(i);
                    pulse.data(i)=[];
                    timeReal = ~isnan(plsInfo.time);
                    valReal = ~isnan(plsInfo.val);
                    for j=1:length(elemInfo) % replace default element with anything specified in pulse. 
                        if ischar(elemInfo(j))
                            elemInfo(j) = struct('type',elemInfo(j),'time',[],'val',[]);
                        end
                        elemInfo(j).time(timeReal)=plsInfo.time(timeReal);
                        elemInfo(j).val(valReal)=plsInfo.val(valReal);
                    end
                    pulse.data = [pulse.data(1:i-1) orderfields(elemInfo,pulse.data(1)) pulse.data(i:end)]; % replace element of data with new version 
                    changed=1; changedout=1;
                    if isfield(pulse,'pardef') && ~isempty(pulse.pardef) % this also expands out dictionary elements with length > 1
                        pulse.pardef = bump_pardef(pulse.pardef,i,length(elemInfo)-1);
                    end
                    break;
                end
            end
        end
        startInd=i;
    end
end
end

function pardef = bump_pardef(pardef, from, by)
tobump=find(pardef(:,1) > from);
pardef(tobump,1)=pardef(tobump,1)+by;
end