function grpdef2 = plsmakegrp(name, ctrl, pulseInds)
% Convert pulses in pulsegroup to wf format.
% grpdef = plsmakegrp(name, ctrl, ind)
% name: group name.
% ctrl: 'plot', 'check', 'upload'
%       for maintenance/debugging: 'clrzero', 'local'.
%       These may mess with the upload logging, so use with care.
% ind: optional pulse index. The default is pulseind or all pulses.

global plsdata; global awgdata;
if ~exist('ctrl','var'), ctrl=''; end
if ~exist('pulseInds','var'), pulseInds=[]; end
if ~iscell(name), name = {name}; end         
for i = 1:length(name)
    if ~isstruct(name{i}) % if group name is a string, then load file with grpdef.
        load([plsdata.grpdir, 'pg_', name{i}]);
    else
        grpdef=name{i};
        if isstruct(name{i}.pulses) && ~isfield(name{i}.pulses(1).data,'marker')
            d=load([plsdata.grpdir, 'pg_', grpdef.name]);            
            grpdef.pulses = d.grpdef.pulses;
        end
    end
    if isopt(grpdef.ctrl, 'seq') % Seems like seq is a cell of groups in a group that gets loaded together. 
        for m = 1:length(grpdef.pulses.groups)
           plsmakegrp(grpdef.pulses.groups{m},ctrl,pulseInds);
        end    
        return;
    end           
    if ~isfield(grpdef, 'varpar'), grpdef.varpar = []; end
    if ~isfield(grpdef, 'params'), grpdef.params = []; end      
    makeType = strtok(grpdef.ctrl); % just the first word of ctrl. pls or grp.
    switch makeType
        case 'pls' % Used for single group
            if isfield(grpdef,'varpar') && size(grpdef.varpar,1)==1 && size(grpdef.varpar,2)>1 % Needs to be a column, not row.
                error('varpar is not a column vector. Did you want to transpose it?!?');
            end
            npars = max(1, size(grpdef.varpar, 1)); % number of varpar values
            origPulses = grpdef.pulses; 
            grpdef.pulses = plsdefault(grpdef.pulses); % Returns pulse info.
            plsDef = grpdef.pulses;
            if ~exist('pulseInds','var') || isempty(pulseInds) % list of pulses?
                if isfield(grpdef, 'pulseind')
                    pulseInds = unique(grpdef.pulseind);
                else
                    pulseInds = 1:length(plsDef)*npars;
                end
            end
            if isfield(grpdef,'dict') && ~isempty(grpdef.dict) % Load dictionary and put in grpdef.
                origDict = grpdef.dict; 
                grpdef.dict=pdpreload(grpdef.dict);
            end
            for j = 1:length(pulseInds) % most of the time, this will give the same pulse for all...
                plsOrd = floor((pulseInds(j)-1)/npars)+1;
                plsDef(j)=grpdef.pulses(plsOrd);
            end
            grpdef = rmfield(grpdef,'pulses');
            for j = 1:length(pulseInds)
                pulseNum = mod(pulseInds(j)-1, npars)+1; % most of the time will be j
                params = grpdef.params; % transfer valid pulse dependent parameters. Interpretation of nan may not be so useful here, but should not hurt.
                if ~isempty(grpdef.varpar) % assign varpar for this group
                    paramMask = ~isnan(grpdef.varpar(pulseNum, :));
                    nVars = size(grpdef.varpar,2);
                    params(end-nVars + find(paramMask)) = grpdef.varpar(pulseNum, paramMask);
                end
                if ~isempty(plsDef(j).trafofn) % assign trafofn for these params
                    params = plsDef(j).trafofn(params);
                end
                if isfield(grpdef,'dict') && ~isempty(grpdef.dict) && strcmp(plsDef(j).format,'elem') % Apply dictionary to pulsedata before params.
                    plsDef(j)=pdapply(grpdef.dict, plsDef(j),pulseInds(j));
                end
                plsOrd = floor((pulseInds(j)-1)/npars)+1;
                if ~isempty(plsDef(j).pardef) % update parameters - could move to plstowf
                    pardef = plsDef(j).pardef;
                    switch plsDef(j).format % applies params to the pardef in pls.
                        case 'elem'
                            for k = 1:size(pardef, 1)
                                if isnan(params(k))
                                    continue;
                                end
                                if pardef(k, 2) < 0
                                    plsDef(j).data(pardef(k, 1)).time(-pardef(k, 2)) = params(k);
                                else
                                    plsDef(j).data(pardef(k, 1)).val(pardef(k, 2)) = params(k);
                                end
                            end
                        case 'tab'
                            for k = 1:size(pardef, 1)
                                if isnan(params(k))
                                    continue;
                                end
                                if pardef(k, 1) < 0
                                    plsDef(j).data.marktab(pardef(k, 2), -pardef(k, 1)) = params(k);
                                else
                                    plsDef(j).data.pulsetab(pardef(k, 2), pardef(k, 1)) = params(k);
                                end
                            end
                        otherwise
                            error('Parametrization of pulse format ''%s'' not implemented yet.', plsDef(plsOrd).format);
                    end
                end
                grpdef.pulses(j) = plstowf(plsDef(j));
                if iscell(grpdef.varpar)
                    grpdef.pulses(j).xval = [grpdef.varpar{plsOrd}(pulseNum, end:-1:1), grpdef.pulses(j).xval];
                elseif ~isempty(grpdef.varpar)
                    grpdef.pulses(j).xval = [grpdef.varpar(pulseNum, end:-1:1), grpdef.pulses(j).xval];
                end
            end
        case 'grp' % Used for pulsing on both sides. 
            groupdef = grpdef.pulses;
            grpdef.pulses = struct([]);
            nchan = size(grpdef.matrix, 2);
            % To allow for virtual awg channels, make a list of all channels in groupdef.groups first. "chan" does double duty--it isboth a list of channels and a vector of indices-JMN
            chanList=[];
            for j = 1:length(groupdef.groups)
                pgList{j}=plsmakegrp(groupdef.groups{j},'upload local');
                chanList=[chanList pgList{j}.chan];
            end
            grpdef.pgList = pgList; 
            for j = 1:length(groupdef.groups) % Across groups
                pg=pgList{j};
                if ~isfield(pg, 'pulseind') % Some flag set to apply pulseind after adding, same for all groups
                    pg.pulseind = 1:length(pg.pulses);
                end
                % target channels for j-th group
                if isfield(groupdef, 'chan') % chan assumed to be an index
                    chan=groupdef.chan(1,j);
                else % chan assumed to be a channel
                    chan = pg.chan;
                    % Find the indices of the virtual channels for the current group from the list. chan now becomes a list of indices-JMN
                    [~, chan]=ismember(chan,chanList); 
                end
                mask = chan > 0;
                chan(~mask) = [];
                
                if ~isfield(groupdef, 'markchan') % target channels for markers
                    markchan = chan;
                else
                    markchan = groupdef.markchan(j, :);
                end
                markmask = markchan > 0;
                markchan(~markmask) = [];
                %ind not given to recursive call above, so plsmagegrp makes all pulses, specified by pulseind or default
                % of source group. Need to reconstruct indices as used for file names by inverting unique
                [~, ~, pind] = unique(pg.pulseind(min(j,end),:));
                for k = 1:length(pg.pulses(1).data)
                    for m = 1:length(pg.pulseind) % Across pulses
                        if j == 1 % first pf determines size
                            grpdef.pulses(m).data(k).wf = zeros(nchan, size(pg.pulses(pind(m)).data(k).wf, 2));
                            grpdef.pulses(m).data(k).marker = zeros(nchan, size(pg.pulses(pind(m)).data(k).wf, 2), 'uint8');
                            grpdef.pulses(m).data(k).readout = pg.pulses(pind(m)).data(k).readout; % a bit of a hack.
                            grpdef.pulses(m).data(k).clk = pg.pulses(pind(m)).data(k).clk;
                            grpdef.pulses(m).xval = [];
                        else
                            for n=1:size(pg.pulses(pind(m)).data(k).readout,1)
                                if isempty(grpdef.pulses(pind(m)).data(k).readout)
                                    grpdef.pulses(m).data(k).readout = pg.pulses(pind(m)).data(k).readout; % a bit of a hack.
                                else
                                    roi=find(grpdef.pulses(pind(m)).data(k).readout(:,1) == pg.pulses(pind(m)).data(k).readout(n,1));
                                    if ~isempty(roi)
                                        fprintf('Overwriting readout window\n');
                                        grpdef.pulses(pind(m)).data(k).readout(roi(1),2:3) = pg.pulses(pind(m)).data(k).readout(n,2:3);
                                    else
                                        grpdef.pulses(pind(m)).data(k).readout(end+1,:) = pg.pulses(pind(m)).data(k).readout(n,1:3);
                                    end
                                end
                            end
                        end
                        grpdef.pulses(m).data(k).wf(chan, :) = grpdef.pulses(m).data(k).wf(chan, :) + pg.pulses(pind(m)).data(k).wf(mask, :);
                        grpdef.pulses(m).data(k).marker(markchan, :) = bitor(grpdef.pulses(m).data(k).marker(markchan, :), ...
                            pg.pulses(pind(m)).data(k).marker(markmask, :));
                        grpdef.pulses(m).xval = [grpdef.pulses(m).xval, pg.pulses(pind(m)).xval];
                    end
                end
            end
            if nargin < 3 || isempty(pulseInds)
                pulseInds = 1:length(grpdef.pulses);
            else
                grpdef.pulses = grpdef.pulses(pulseInds);
            end
            [grpdef.pulses.format] = deal('wf');
        otherwise
            error('Group control %s not understood.\n',grpdef.ctrl);
    end
    if isfield(grpdef, 'xval') && ~isempty(grpdef.xval) % add groups xval to pulse xval.
        for j=1:length(grpdef.pulses)
            grpdef.pulses(j).xval = [grpdef.xval grpdef.pulses(j).xval];
        end
    end
    for j = 1:length(pulseInds) % Inc. matrix, offset, trafofn
        for k=1:length(grpdef.pulses(j).data) % multiply by matrix and add offset
            npulses = size(grpdef.pulses(j).data(k).wf, 2);
            grpdef.pulses(j).data(k).wf = grpdef.matrix * (grpdef.pulses(j).data(k).wf + repmat(grpdef.offset, 1, npulses));
            if isfield(grpdef, 'trafofn') && ~isempty(grpdef.trafofn) % transform data
                wf=grpdef.pulses(j).data(k).wf;
                for m=1:length(grpdef.trafofn)
                    fn=grpdef.trafofn(m).func;
                    args=grpdef.trafofn(m).args;
                    if ~iscell(args)
                        args={args};
                    end
                    for q=1:size(wf,1)
                        wf(q,:) = fn(wf,q,args{:}); %JMN, 2016/05/11 this will break a lot of stuff.
                    end
                end
                grpdef.pulses(j).data(k).wf=wf;
            end
            if isfield(grpdef, 'markmap')
                markerData = grpdef.pulses(j).data(k).marker;
                grpdef.pulses(j).data(k).marker = zeros(size(grpdef.matrix, 1), size(markerData, 2), 'uint8');
                grpdef.pulses(j).data(k).marker(grpdef.markmap(2, :), :) = markerData(grpdef.markmap(1, :), :);
            end
        end
    end
    grpdef.ctrl = ['pls', grpdef.ctrl(find(grpdef.ctrl == ' ', 1):end)]; % now make pls the first control..
    firstCtrl = strtok(ctrl);
    grpdef2 = grpdef;
    switch firstCtrl
        case 'plot'
            if isfield(grpdef,'dict') && ~isempty(grpdef.dict)
                plsplot(grpdef.pulses,grpdef.dict,ctrl);
            else
                plsplot(grpdef.pulses,[],ctrl);
            end
        case 'check'
            for j = 1:length(pulseInds)
                over=0;
                for a=1:length(awgdata)
                    for k=1:length(grpdef.pulses(j).data)
                        s=min(awgdata(a).scale);
                        over = over || any(abs(grpdef.pulses(j).data(k).wf(l,:)) > s);
                    end
                end
                if over
                    fprintf('Pulse %d is too large\n',j);
                end
            end
        case 'upload'
            if isopt(ctrl, 'force') || plsinfo('stale',grpdef.name) % Modified since last upload (or upload forced)
                if isopt(grpdef.ctrl,'pack')  % A little naughty; secretly pack all the pulse waveforms together for load...
                    if any(~strcmp('wf',{grpdef.pulses.format}))
                        error('Pack can only deal with waveforms.');
                    end
                    packdef = grpdef;
                    packdef.pulses=[];
                    packdef.pulses(1).format='wf';
                    for j=1:length(grpdef.pulses(1).data)
                        data=vertcat(grpdef.pulses.data);
                        data=data(:,j);
                        packdef.pulses(1).data(j).marker = [data.marker];
                        packdef.pulses(1).data(j).wf = [data.wf];
                        packdef.pulses(1).data(j).clk = data(1).clk;
                    end                    
                else
                    packdef = grpdef;
                end
                if ~isopt(ctrl, 'local') % Actually handle the upload... % awgload/zero doesn't use anything else.
                    zerolen = awgload(packdef, pulseInds);
                else
                    zerolen = awgzero(packdef, pulseInds);
                end
                if isopt(grpdef.ctrl,'pack')
                    for j=1:length(zerolen)
                        zerolen{j}=repmat(zerolen{j}(1,:)/length(grpdef.pulses),length(grpdef.pulses),1);
                    end
                end
                if isfield(grpdef.pulses(1).data,'readout')
                    readout=[];
                    for j=1:length(grpdef.pulses)                        
                        if ~isempty(grpdef.pulses(j).data(1).readout) % each row is different readout event (most pulses will only have one).
                            readout(:,:,j) = grpdef.pulses(j).data(1).readout; %#ok<*AGROW>
                            plens(j) = size(grpdef.pulses(j).data(1).marker,2);
                        end
                    end
                    chg=[];
                    for j=1:size(readout) % first dimension?
                        tmp=find(squeeze(any(abs(diff(readout(j,:,:),[],3))>1e-4)));
                        chg=union(chg,tmp);
                    end
                    chg=sort(chg);
                    if isempty(chg)
                        rd.readout = readout(:,:,1);
                        rd.reps = size(readout,3); % number of pulses?
                        rd.plens = size(grpdef.pulses(1).data(1).marker,2); %for whatever eras
                    else
                        rd.readout=readout;
                        rd.reps=ones(1,size(readout,3));
                        rd.plens=plens;
                    end
                    if any(squeeze(any(abs(diff(readout,[],3)) > 1e-10)))
                        warning('plsmakegrp:readoutTime','Readout changes between pulses in %s',name{i});
                    end
                    grpdef.readout=rd;
                end
                for j=1:length(grpdef.pulses)
                    l=length(grpdef.pulses(j).xval(:));
                end
                grpdef.ind = pulseInds;
                grpdef.zerolen = zerolen;                                
                grpdef2 = grpdef;  
                grpdef.pulses=[];
                if exist('origPulses','var')
                    grpdef.pulses = origPulses; 
                else
                    %grpdef.pulses=pulseNum;%origPulses;
                end                                                
                if ~isempty(grpdef.dict)
                    grpdef.dict = origDict;
                end
                save([plsdata.grpdir, 'pg_', name{i}], 'grpdef'); %2012/12/26
                logentry('Uploaded group %s', grpdef.name);
            else
                % Save dictionary 
                grpdef.dict = origDict; %grpdef.pulses = origPulses; 
                %grpdef.pulses = pulseNum;
                grpdef2 = grpdef;
            end
    end
end
end