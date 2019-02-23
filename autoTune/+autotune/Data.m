classdef Data < dynamicprops
    %autotune.Data < dynamic props. represents tuneadta
    %   see help for individual properties and methods.
    %   e.g. >> help autotune.Data.activeSetName
    %   e.g. >> help autotune.Data.newRun
    
    properties
        dir=''; %directory to save information
        activeSetName = 'left'; % name of active set
        alternates; %array of alternates to swap in
        basis = diag([-ones(1,4),ones(1,12)]); % gradient matrix
        baseNames={'XL','YL','XR','YR','Lead1','Lead2','Lead3','Lead4','N12','N34','T12','T34','pXL','pYL','pXR','pYR'}; 
        gateChans = {'2a','1a','3a','4a','1b','2b','3b','4b','N12','N34','T12','T34','PlsRamp2','PlsRamp1','PlsRamp3','PlsRamp4'};        
        measPt = [0,0]; %current meas point
        sepDir = [1 -1]; %direction of sep. should be dict.sep.val(1:2_
        displayFig = 3; % where to plot things
        axes; 
        numAxes = [3,3]; % will call subplot(numAxes(1),numAxes(2),xx)
        dataChan = 'DAQ1'; %default getchan for all scans
        xyChan = {'PlsRamp2','PlsRamp1'}; % X and Y axes in charge scan
        xyBasis = {'XL','YL'}; %basis vectors for X and Y directions (for centering)
        figHandle = 3; %where to plot stuff   
        basisStore; 
        file; 
    end
    properties (Dependent = true)
        runNumber; %current run number. dependent property get calculated for you
    end    
    properties (Constant, Hidden)
        doNotSwap = {'dir','alternates','basis','baseNames','runNumber','gateChans','file'}; %properties not to swap
    end    
    methods
        function this = Data(side)
            if exist('side','var') && strcmp(side,'right')
                 this.activeSetName = 'right';
                 this.displayFig = 4;
                 this.dataChan = 'DAQ2';
                 this.xyChan = {'PlsRamp3','PlsRamp4'};
                 this.xyBasis = {'XR','YR'};
            end
        end
       
        function num = get.runNumber(this) 
            %function num = get.runNumber(this)
            %get the run number. used internally for dependent property
            if isprop(this,'chrg') && ~isempty(this.chrg)
               num = size(this.chrg.trTriple,1);
            else
                fprintf('Confused about run number. Is tuneData empty?');
                num = 0;
            end
        end
        
        function equalizeRuns(this)
            % If the number of runs becomes different for different runs
            % (e.g. if something crashed), run this to fix. 
            runNumber = this.runNumber; 
            if size(this.lead.timeX,1)>runNumber
                this.lead.timeX(runNumber+1:end,:) = []; %#ok<*MCNPR>
                this.lead.posX(runNumber+1,:) = [];
                this.lead.timeY(runNumber+1,:)=[];
                this.lead.posY(runNumber+1,:)=[];
            elseif size(this.lead.timeX,1)<runNumber
                for i = size(this.lead.posY,1):runNumber
                    this.lead.timeX(i,:) = nan;
                    this.lead.posX(i,:) = nan;
                    this.lead.timeY(i,:)=nan;
                    this.lead.posY(i,:)=nan;
                end
            end            
                        
            if size(this.zoom.measPt,1)>runNumber
                this.zoom.measPt(runNumber+1:end,:) = [];
            elseif size(this.zoom.measPt,1)<runNumber
                for i = size(this.zoom.measPt,1):runNumber
                    this.zoom.measPt(i,:) = nan;
                end
            end       
            
            if length(this.tl.location)>runNumber
                this.tl.location(runNumber+1:end) = [];
                this.tl.width(runNumber+1:end) = [];
                this.tl.widtherr(runNumber+1:end) = [];
            elseif length(this.tl.location)<runNumber
                for i = length(this.tl.location):runNumber
                    this.tl.location(i) = nan;
                    this.tl.width(i) = nan;
                    this.tl.widtherr(i) = nan;
                end
            end
            
            if length(this.stp.location)>runNumber
                this.stp.location(runNumber+1:end) = [];
                this.stp.width(runNumber+1:end) = [];
                this.stp.widtherr(runNumber+1:end) = [];
            elseif length(this.stp.location)<runNumber
                for i = length(this.stp.location):runNumber
                    this.stp.location(i) = nan;
                    this.stp.width(i) = nan;
                    this.stp.widtherr(i) = nan;
                end
            end
            
            if length(this.t1.t1)>runNumber
                this.t1.t1(runNumber+1:end) = [];
            elseif length(this.t1.t1)<runNumber
                for i = length(this.t1.t1):runNumber
                    this.t1.t1(i) = nan;
                end
            end
            
            if length(this.loadTime.time)>runNumber
                this.loadTime.time(runNumber+1:end) = [];
                this.loadTime.amp(runNumber+1:end) = [];
            elseif length(this.loadTime.time)<runNumber
                for i = length(this.loadTime.time):runNumber
                    this.loadTime.time(i) = nan;
                    this.loadTime.amp(i) = nan;
                end
            end
            
            if length(this.loadPos.width)>runNumber
                this.loadPos.width(runNumber+1:end) = [];
            elseif length(this.loadPos.width)<runNumber
                for i = length(this.loadPos.width):runNumber
                    this.loadPos.width(i) = nan;
                end
            end
            
            if length(this.line.width)>runNumber
                this.line.width(runNumber+1:end) = [];
            elseif length(this.line.width)<runNumber
                for i = length(this.line.width):runNumber
                    this.line.width(i) = nan;
                end
            end
        end
        
        function newRun(this)
            % function newRun(this)
            % make a new run
            % go through all of the properites that are 'autotune.Op' and
            % ask each to make a new run
            % set measPt to 0,0; clear figure
            propList =  properties(this);
            newRunNumber = this.runNumber+1;
            for j = 1:length(propList)
                if isa(this.(propList{j}),'autotune.Op')
                    this.(propList{j}).makeNewRun(newRunNumber);
                end
            end
            this.measPt = [0,0];
            figure(this.figHandle); clf;
            this.rePlot; 
            fprintf('Making new tune run %i %s: \n',this.runNumber,this.activeSetName);
        end
        
        function items = rePlot(this,num,opts)            
            %function items = rePlot(this,num,opts)
            % Create a new axis for plots. Replot all the most recent data.
            % If given the pulse option, will also plot the pulses of all
            % the tuneData pulses. 
            if ~exist('opts','var'), opts = ''; end
            plotSpace = {this.numAxes(1),this.numAxes(2), [0.063, 0.073], [0.06 0.045], [0.06, 0.1]};            
            figure(2); clf; 
            if ~exist('num','var') || isempty(num)
                this.chrg.ana('last auto');
                outZoom=this.zoom.ana('last noset');
                figure(this.displayFig); clf;
                this.axes = tight_subplot(plotSpace{:});
                outTL=this.tl.ana('last');
                outSTP=this.stp.ana('last');
                this.line.ana('last');
                this.loadTime.ana('last');
                outload=this.loadPos.ana('last');
                this.lead.ana('last');
                num = this.runNumber; 
            else
                this.chrg.ana('auto',num);
                outZoom=this.zoom.ana('noset',num);
                figure(this.displayFig); clf;
                this.axes = tight_subplot(plotSpace{:});
                outTL=this.tl.ana('',num);
                outSTP=this.stp.ana('',num);
                this.line.ana('',num);
                outload=this.loadPos.ana('',num);
                this.loadTime.ana('',num);
                this.lead.ana('',num);
            end
            if isopt(opts,'pulse')
                figure(10); clf;
                plotSpace = {2,2, [0.063, 0.05], [0.06 0.04], [0.06, 0.06]};
                plsAxes = tight_subplot(plotSpace{:});
                plotSpace2 = {4,2};%, [0.063, 0.05], [0.06 0.04], [0.06, 0.06]};
                figure(11); plsAxes2 = tight_subplot(plotSpace2{:});
                items =false; 
                if ~isempty(fieldnames(outZoom))
                    atplschk(this.zoom.plsGrp,this.activeSetName,struct('pulses',[2,1],'axis',[plsAxes(1); plsAxes2(1:2)],'run','num','title','Zoom'));
                    items = true; 
                end
                if ~isempty(fieldnames(outload))
                    atplschk(this.loadPos.plsGrp,this.activeSetName,struct('pulses',[51,26,1],'axis',[plsAxes(2); plsAxes2(3:4)],'offset',-outload.measPt,'run',num,'title','LoadPos','time',outload.time));                    
                    items =true;
                end                
                if ~isempty(fieldnames(outSTP))
                    atplschk(this.stp.plsGrp,this.activeSetName,struct('pulses',[100,50,1],'axis',[plsAxes(3); plsAxes2(5:6)],'offset',-outSTP.measPt,'run',num,'title','STP')); 
                    items = true; 
                end
                if ~isempty(fieldnames(outTL))
                    atplschk(this.tl.plsgrp,this.activeSetName,struct('pulses',[100,50,1],'axis',[plsAxes(4); plsAxes2(7:8)],'offset',-outTL.measPt,'title','TL')); 
                    items = true; 
                end
            end
        end
        
        function start(this,center) 
            % function start(this,center) 
            % Begin using autotune. Give the center of junction. Turn off
            % autoramping. Set bottom sensor gate to correct value. Run a
            % sensor and charge scan. 
            global scandata; 
            smset([scandata.sens.loops(1).setchan(1),scandata.sens.loops(2).setchan(1)],center)
            smset(scandata.sensor.loops(2).setchan{1},autoscan('sens',struct('sensorVal',center)));       
            scandata.autoramp = 0; 
            this.sensor.run; 
            this.chrg.run;                      
        end
        
        function center(this,opts) 
            % function center(this,opts) (replaces atcenter)
            % center things so measPt = [0,0]
            % opts: 'quiet' to supress output
            %       'noconfirm'
            %       'offset' will center things by change awgdata.offset to be the right value. 
            if ~exist('opts','var'), opts = ''; end            
            if isopt(opts,'offset')
                global awgdata; %#ok<TLEV>
                c1=str2double(this.xyChan{1}(end)); %figure out which plsramp
                c2=str2double(this.xyChan{2}(end));
                newoffset=awgdata(1).offset;
                fprintf('Old awg offsets are [%s]\n',sprintf('%g ',newoffset));
                newoffset([c1 c2]) = newoffset([c1 c2]) + 1e3*this.measPt;
                fprintf('New awg offsets will be [%s]\n',sprintf('%g ',newoffset));
                if strcmp(input('Accept (yes/[no])? ','s'), 'yes') == 0
                    fprintf('not centering \n');
                    return;
                end
                
                for j = 1:length(awgdata)
                    awgdata(j).offset(awgdata(j).chans)=newoffset(awgdata(j).chans);
                end
                if ~isopt(opts,'quiet')
                   fprintf('Consider repacking waveforms') 
                end
            else
                cx = -this.measPt(1); 
                cy = -this.measPt(2); 
                if all([cx,cy] == 0)
                    fprintf('No change\n');
                    return;
                end
                if ~isopt(opts,'quiet')
                    fprintf('Step of %g mV on %s, %g mV on %s\n', cx *1e3, this.xyBasis{1},cy * 1e3,this.xyBasis{2});
                end
                if ~isopt(opts,'noconfirm') && (strcmp(input('Accept (yes/[no])? ','s'), 'yes') == 0)
                    fprintf('Not centering \n');
                    return;
                end
                if ~isopt(opts,'quiet')
                    fprintf('atchg(''%s'',%g)\n',this.xyBasis{1},cx);
                    fprintf('atchg(''%s'',%g)\n',this.xyBasis{2},cy);
                end
                this.tmp.lastCenter(1) = cx; 
                this.tmp.lastCenter(2) = cy; 
                this.change(this.xyBasis{1},cx);
                this.change(this.xyBasis{2},cy);
                this.measPt= 0*this.measPt; % keep the size the same;
            end
        end

        function undoCenter(this)
            %function undoCenter(this)
            % If yout don't like what the last center did, undo it. 
            this.change(this.xyBasis{1},-this.tmp.lastCenter(1));
            this.change(this.xyBasis{2},-this.tmp.lastCenter(2));
        end
        
        function runInfo = getRunInfo(this,ind)
           %function runInfo = getRunInfo(this,ind)
           %return a struct will information from run ind
            runInfo = struct(); 
           propList =  properties(this);
           for j = 1:length(propList)
               if isa(this.(propList{j}),'autotune.Op')
                   runInfo.(propList{j}) = this.(propList{j}).getData(ind);
               end
           end
        end
        
        function runAll(this,opts) 
         if ~exist('opts','var'), opts = ''; end   
            this.chrg.run; 
            this.line.run; 
            this.lead.run; 
            this.zoom.run 
            if ~isopt(opts,'nomeas')
                this.loadPos.run; 
                this.loadTime.run; 
                this.stp.run;
                this.tl.run;
                this.t1.run;                 
            end
        end
        
        function updateAll(this,opts)
            % Update zoom, loadPos, loadTime, stp, tl, t1 scans. Add feedback groups. Turn on AWG. 
            % If you don't give the nodict option, will update dictionaries. 
            if ~exist('opts','var'), opts = ''; end
            if isopt(opts,'all') 
                global awgdata
                awgdata.zeropls = [];
                awgrm('all'); awgclear('unusued'); 
                awgadd('all_off_LR'); 
                awgadd('chrg_1_L'); 
                awgadd('sqrX_L'); 
                awgadd('sqrY_L');                 
            end
            try
                awgrm(4,'after');
                awgclear('unused');
            catch
            end
            if this.sepDir(1) == 1
                side = 'TL';
            else
                side = 'BR'; 
            end
            if ~isopt(opts,'nodict')
                updateExch(side,struct('opts','all'));
                this.loadPos.updateGroup('target');
                this.stp.updateGroup('target');
                this.tl.updateGroup('target');
            else
               this.loadPos.updateGroup; 
               this.stp.updateGroup;
                this.tl.updateGroup;
            end
            this.zoom.updateGroup                        
            this.loadTime.updateGroup;            
            this.t1.updateGroup;    
            feedbackTest('all'); 
            awgcntrl('on start amp err');
        end
        
        function addGroups(this,opts)
            %function addGroups(this,opts) 
            % opts: start 
            % Adds all_off, chrg, sqr pulses. 
            if ~exist('opts','var'), opts = ''; end
            pg.ctrl = 'loop';
            pg.nrep = 0;
            if isopt(opts,'start')
                namepat = 'all_off_';
                pg.pulses = 1;
                pg.params=[];
                pg.varpar = [];
                pg.chan = [2 1];
                pg.dict = 'left';                
                pg.name=[namepat 'LR'];
                plsdefgrp(pg);
                awgadd(pg.name);
            end
            side = upper(this.activeSetName(1));
            chan=[str2double(char(regexp(this.xyChan{1},'\d+','match'))),str2double(char(regexp(this.xyChan{2},'\d+','match')))];
            % chrg scan             
            
            pg.name = ['chrg_1_' side];
            pg.chan = chan; 
            pg.pulses =2;                        
            pg.dict={this.activeSetName};
            plsdefgrp(pg);
            awgadd(pg.name);                        
            
            clear pg            
            pg.name = ['sqrX_' side];
            pg.pulses = 3;
            pg.chan = chan; 
            plsdefgrp(pg);            
            awgadd(pg.name); 
            
            pg.name = ['sqrY_' side];
            pg.pulses = 4;
            pg.chan = chan; 
            plsdefgrp(pg);            
            awgadd(pg.name); 
        end
        
        function new = copy(this)
            % function new = copy(this)
            % Instantiate new object of the same class.
            % this is used because autotune.Data is a handle class
            % a = copy(tuneData); will make a completely new copy
            % a = tuneData; if you now modify a, tuneData is modified;
            new(length(this)) = feval(class(this));          
            for ind = 1:length(this)
                % Copy all non-hidden properties.
                %p = setdiff(properties(this(ind)),properties(class(this)));
                p = setdiff(properties(this(ind)),{'runNumber','doNotSwap'});
                for i = 1:length(p)
                    if ~isprop(new(ind),p{i})
                        addprop(new(ind),p{i});
                    end
                    new(ind).(p{i}) = this(ind).(p{i});
                end
            end
        end
        
        function new = swapCopy(this)
            %function new = swapCopy(this)
            %used internally. makes a copy with only some of the properties
            % see autotune.Data.Copy
            new(length(this)) = feval(class(this));          
            for ind = 1:length(this)
                p = setdiff(properties(this(ind)),autotune.Data.doNotSwap);
                %p = setdiff(properties(this(ind)),properties(class(this)));
                %p = [p',{'activeSetName','measPt','displayFig'}];
                for i = 1:length(p)
                    if ~isprop(new(ind),p{i})
                        addprop(new(ind),p{i});
                    end
                    new(ind).(p{i}) = this(ind).(p{i});
                end
            end
        end
        
        function swap(~,newActiveSet)
            error('broken. Use autotune.swap(''%s'')',newActiveSet);            
        end
        
        function copyRun(this,runNum,target)
            % copy relevant information from run number runNum to target
            % if target not given assumed = tuneData.runNumber+1
            % if runNum also not given assumed = tuneData.runNumber
            switch nargin
                case 3 % runNum and target given

                case 2 % only runNum given;
                    target = this.runNumber+1;
                    this.newRun();
                case 1 %neither runNum not target given
                    target = this.runNumber+1;
                    runNum = this.runNumber;
                    this.newRun;
                otherwise
                    error('improper usage')
            end
            if numel(target)>1
               error('multi-dimensional copying not implemented') 
            end
            if target == this.runNumber+1
                this.newRun();
            end
            data = this.getRunInfo(runNum);
            ops = fieldnames(data); %different operations
            for o = 1:length(ops)
                fn = fieldnames(data.(ops{o}));
                for f = 1:length(fn)
                    if all(size(data.(ops{o}).(fn{f}))~=1)
                       error('multi-dimensional indexing not implemented') 
                    end
                    if numel(data.(ops{o}).(fn{f}))==1
                        this.(ops{o}).(fn{f})(target) = data.(ops{o}).(fn{f});
                    else %fixme, this will break if the field is something x Nrun istead of Nrun x something
                        this.(ops{o}).(fn{f})(target,:) = data.(ops{o}).(fn{f});
                    end
                end
            end
        end
        
        function change(this,gate,dv)
            %function change(this,gate,dv)            
            % changes the gates along a basis direction indicated by gate by an amount dv. gate can be number (index into
            % tuneData.basis) or name (must be part of tuneData.baseNames
            % example: tuneData.change('Lead2',5e-3);
            
            if abs(dv) > 25e-3 || isnan(dv)
                fprintf('Crazy big change.  I don''t believe you; ignoring.\n');
                return;
            end            
            if ischar(gate)
                gatename = gate;
                gate = strcmp(gate, this.baseNames);
                if isempty(find(gate,1))
                    error('Cannot find basis direction named %s',gatename);
                end
            end            
            v = smget(this.gateChans);
            smset(this.gateChans, [v{:}] + this.basis(:, gate)' * dv);
        end
                
        function print(this,opts)
            % Print the basis. If an option of a side given, only prints
            % one side. 
            if ~exist('opts','var'), opts = ''; end
            if ~isempty(opts)
                if isopt(opts,'right')
                    gates = [3,4,7,8,10,12,15,16];
                elseif isopt(opts,'left')
                    gates = [1,2,5,6,9,11,13,14]; 
                end
                nGates = length(gates);
                    fprintf(['%-12s',repmat('%-9s', 1, nGates+1), '\n'], '',this.gateChans{gates});
                    fprintf('\n');
                    fprintf('------------------------------------------------------------------------------------------------------------------------\n')
                    for i = 1: nGates
                        n = gates(i);
                        fprintf(['%-9s:', repmat('%8.3g', 1, nGates), '\n'], this.baseNames{n}, this.basis(gates, n));
                    end                
            else
                nGates = length(this.gateChans);
                fprintf(['%-10s',repmat('%-8s', 1, nGates+1), '\n'], '',this.gateChans{:});
                fprintf('\n');
                fprintf('------------------------------------------------------------------------------------------------------------------------\n')
                for i = 1:size(this.basis,2)
                    fprintf(['%-9s:', repmat('%8.3g', 1, nGates), '\n'], this.baseNames{i}, this.basis(:, i));
                end
            end                
        end
        
        function dotCenter(this)
            % Run center repeatedly until centered. 
            dv=60e-6;
            count = 0;
            oldMeasPt = this.measPt;           
            while (norm(oldMeasPt) > dv && count < 25) || count == 0 
               this.stp.run();
               this.tl.run();
               this.tmp.run();                                              
                 if this.stp.foundSTP == 0 || this.tl.foundTL==0
%                     this.chrg.run([],'auto'); % add error checking for these. 
                      fprintf('STP or TL did not fit \n'); 
                      break 
%                     if ~this.chrg.fitTrip 
%                         %possibly write in if there is an issue with the sensor -- what sensitivity 
%                         fprintf('Charge scan did not fit \n')                        
%                         sleep 
%                         return
                  end
%                     this.center;                     
%                     this.dotCenter; 
%                 end
                oldMeasPt = this.measPt; 
                if norm(this.measPt) > 3*norm(oldMeasPt) && count > 0 
                    fprintf('This is not centering. Consider remaking basis. \n');
                    sleep
                    return
                end                
                if norm(this.measPt) < dv
                    fprintf('Well centered after %01d runs \n',count);                                        
                end               
                this.center('quiet noconfirm');
                count = count+1;                
            end
            sleep
        end        

        function restore(this,num)
            % Restore gate vals to a given charge run.
            global tuneData; 
            if exist('num','var') && ~isempty(num)
                sprintf('%s/sm_chrg%s_%03i', this.dir, upper(tuneData.activeSetName(1)),num);
            else
                file = uigetfile([this.dir '\*']);
            end
            load([this.dir '\' file], 'configvals', 'configch');
            if ~empty(tuneData.basisStore{num})
                tuneData.basis = tuneData.basisStore{num}; 
            end
%             configch = smchanlookup(configch); %#ok<*NODEF>
%             channels = smchanlookup(this.gateChans);
%             if ~all(ismember(channels, configch))
%                 warning('WARNING: some channel values not found.\n');
%             end
%             mask = ismember(configch, channels);
%             configch = configch(mask);
%             configvals = configvals(mask);
%             
            smset(configch, configvals);
        end
    end
end
