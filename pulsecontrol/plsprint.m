function plsprint(pls)
% print out elems of pls. Takes number or char.
global plsdata
if ischar(pls)
    plsNames = {plsdata.pulses.name};
    pls = find(strcmpi(plsNames,pls));
elseif isnumeric(pls)
else
    error('Format not supported')
end

if length(pls) > 1
    for i = 1:length(pls)
        plsprint(pls(i))
    end
    return
end

plsInfo = plsdata.pulses(pls);
fprintf('%d: %s \n', pls, plsInfo.name);
if ~isempty(plsInfo.trafofn) && ~isempty(plsInfo.pardef)
    trafo = func2str(plsInfo.trafofn);
    stBrc = strfind(trafo,'[');
    endBrc = strfind(trafo,']');
    if ~isempty(stBrc)
        trafo = trafo(stBrc(1)+1:endBrc(end)-1);
        inStBrc = strfind(trafo,'[');
        inEndBrc = strfind(trafo,']');
    end
    commas = strfind(trafo,',');    
    par = 1;
    if length(commas) ~= size(plsInfo.pardef,1)-1 && ~isempty(commas)
        fprintf('Nonstandard trafofn, beware \n %s \n',trafo)
        par = 0;
    end
else
    par = 0;
end

if strcmp(plsInfo.format,'elem')
    for i = 1:length(plsInfo.data)
        fprintf('%s \t',plsInfo.data(i).type);
        if ~isempty(plsInfo.data(i).time)
            fprintf('%2.2f ns ', plsInfo.data(i).time*1e3)
        end
        if par
            timePardef = plsInfo.pardef(:,1)==i & plsInfo.pardef(:,2)<0;
            if any(timePardef)
                timeInds = find(timePardef);
                for j = 1:length(timeInds)
                    n = timeInds(j);
                    if isempty(commas) || isempty(stBrc) || isempty(endBrc)
                        fprintf('T%d: %s. ',-plsInfo.pardef(n,2),trafo);
                    else
                        if n == 1
                            fprintf('T%d: %s. ',-plsInfo.pardef(n,2),trafo(1:commas(1)-1));
                        elseif n == size(plsInfo.pardef,1)
                            fprintf('T%d: %s. ',-plsInfo.pardef(n,2),trafo(commas(end)+1:end));
                        else
                            fprintf('T%d: %s. ',-plsInfo.pardef(n,2),trafo(commas(n-1)+1:commas(n)-1));
                        end
                    end
                end
            end
        end
        if ~isempty(plsInfo.data(i).val)
            fprintf('at %2.2 mV',plsInfo.data(i).val)
        end
        if par
            valPardef = plsInfo.pardef(:,1)==i & plsInfo.pardef(:,2)>0;
            if any(valPardef)
                valInds = find(valPardef);
                done=false;
                for j = 1:length(valInds)
                    if done
                        done = false;
                        continue
                    end
                    n = valInds(j);
                    if isempty(commas) || isempty(stBrc) || isempty(endBrc)
                        fprintf('T%d: %s. ',-plsInfo.pardef(n,2),trafo);
                    else
                        if n == 1
                            fprintf('V%d: %s. ',plsInfo.pardef(n,2),trafo(1:commas(1)-1));
                        elseif n == size(plsInfo.pardef,1)
                            fprintf('V%d: %s. ',plsInfo.pardef(n,2),trafo(commas(end)+1:end));
                        else
                            if ~isempty(inStBrc) && any(commas(n) > inStBrc & commas(n) < inEndBrc)
                                done = true;
                                brc = commas(n) > inStBrc & commas(n) < inEndBrc;
                                fprintf('V%d: %s. ',plsInfo.pardef(n,2),trafo(commas(n-1)+1:inEndBrc(brc)));
                            else
                                fprintf('V%d: %s. ',plsInfo.pardef(n,2),trafo(commas(n-1)+1:commas(n)-1));
                            end
                        end
                    end
                end
            end
        end
        fprintf('\n');
    end
end
fprintf('\n')
end