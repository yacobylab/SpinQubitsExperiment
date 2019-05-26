classdef Tl < autotune.Op
    %TL < autotune.Op % represents the search for the T+ load
    %   see help for individual properties and methods.
    %   e.g. >> help autotune.TL.plsgrp
    %   e.g. >> help autotune.TL.Run
    
    properties
        plsgrp = 'topLead_1_L'; %plsgrp for scan 
        nLoop = 500; %nloop for scan
        nRep = 10; %nrep for scan
        target = 50;
        search = struct('range',2e3,'points',100,'time',200); %struct to describe parameters of TL search
        sweep = struct('range',300,'time',0.3,'offset',50); %struct to describe the params of the sweep (for FB)
        slope = -.3;%1e-6*[.2, 1]; %copy from tuneData.chrg.(x,y)LeadSlope
        dist = 1; %distance in mV away from triple pt to sweep
        subPlot = 7; %for plotting in tuneData.figHandle
        fitFunc = '@(p,x) p(1)+p(7)*x+p(2)*(tanh((x-p(3))/p(4))+1)/2 - p(5)*(tanh((x-p(6))/p(4))+1)/2';
        fineIndex = 1; 
    end
    
    properties (SetAccess= {?autotune.Data, ?autotune.Op})
        location = 1000; %location in uV  of tl peak, (Nrun x 1) double 
        width; %width of peak in uV, (Nrun x 1) double
        widtherr;
        foundTL;
    end
    
    methods
        function this = Tl
            global tuneData
            if strcmp(tuneData.activeSetName,'right')
                this.plsgrp = 'topLead_2_R';                
            end
        end     
        
        function out = getData(this,runNumber)
            out = struct();
            out.location = this.location(runNumber);
            out.width= this.width(runNumber);
        end
        
        function makeNewRun(this,runNumber)
            if runNumber > 2
                if runNumber ~= length(this.location)+1 || runNumber ~= length(this.width)+1
                    warning('runNumber not consistent with know chrg runs');
                end
                this.location(end+1) = nan;
                this.width(end+1) = nan;
                this.widtherr(end+1) = nan;
                this.fineIndex = 1;           
            end
        end
        
        function run(this)
            global tuneData;
            if ~awgcntrl('ison') 
                awgcntrl('on start wait err');
            end
            scan = fConfSeq(this.plsgrp,{'nloop',this.nLoop,'nrep',this.nRep, 'datachan',tuneData.dataChan,'opts','ampok'});
            scan.consts(end+1).setchan=tuneData.xyChan{1};
            scan.consts(end).val=tuneData.measPt(1);
            scan.consts(end+1).setchan=tuneData.xyChan{2};
            scan.consts(end).val=tuneData.measPt(2);
            scan.loops(1).stream = 1; 
            file = sprintf('%s/sm_tl%s_%04i_%03i',tuneData.dir, upper(tuneData.activeSetName(1)),tuneData.runNumber,this.fineIndex);
            this.fineIndex = this.fineIndex+1; 
            data = smrun(scan, file);
            data = data{1};            
            this.ana('',data,scan);
        end
        
        function out=ana(this,opts,data,scan)
            global tuneData; global fbdata; 
            if ~exist('opts','var'), opts = ''; end
            runNumber = tuneData.runNumber;
            if ~exist('data','var') || isempty(data) || ischar(data) || numel(data)==1 % Check if loading old scan or new data
                if (~exist('data','var') || isempty(data)) && ~isopt(opts,'last')
                    [data,scan,~,time]=loadAna('sm_tl*');
                    out.time = time;
                elseif exist('data','var') && ~isempty(data) && ischar(data)
                    [data,scan,~,time]=loadAna(data);
                    out.time = time;
                else
                    side = upper(tuneData.activeSetName(1));
                    if isopt(opts,'last')
                        data = tuneData.runNumber;
                    end
                    fileName = sprintf('sm_tl%s_%04.f_001.mat',side,data);
                    [data,scan,~,time]=loadAna(fileName);
                    if isempty(data)
                        out=struct;
                        return
                    else
                        out.time = time;
                    end
                end
                anaData=1;
                out.scan = scan; 
            else
                anaData=0;
            end
            eps = scan.data.pulsegroups.varpar'; % tl scan sweeps epsilon value (along TL curve)
            data=mean(data,1); % Average multiple lines of data. 
            fbdata.refval(1) = mean(data); % Set up refval, used for histogramming. 
            [~,ci]=max(data);
            epsMax=eps(ci); % TL is centered around where data largest. 
            axes(tuneData.axes(this.subPlot)); cla; 
            % 1: offset, 2: upward slope, 3: upward center, 4: width both, 5: downward slope, 6: downward center
            beta0 = [min(data), range(data), epsMax-range(eps)/6, range(eps)/3, range(data)/2, epsMax+range(eps)/6,-2e-6];
            try
                params=fitwrap('plinit plfit samefig woff',eps,data,beta0,this.fitFunc);
                [params,~,~,~,~,err]=fitwrap('plinit plfit samefig woff',eps,data,params,this.fitFunc);
                tlPt = (params(3)+params(6))/2; % TL point is center point between 2 tanh fcns.
                tlWid = params(6)-params(3);
                err = err(1,4,1);                
                if ~anaData % If new data, add fit info to tuneData.
                    this.width(runNumber) = tlWid;
                    this.location(runNumber) = tlPt;
                    this.widtherr(runNumber) = err;
                    this.foundTL = err < 100;
                end
                title(sprintf('TL %3.1f, wdth %3.1f',tlPt,err));
            catch
            end
            a = gca; a.YTickLabelRotation=-30;
            a.YLabel.Position(1) = a.XLim(1) - range(a.XLim)/14;
            a.XLabel.Position(2) = a.YLim(1) - range(a.YLim)/7;
            a.XLim = [min(eps),max(eps)];            
            %            if params(6) > params(3)+ 800
            %               tlpt = eps(end)-100;
            %              fprintf('Broad peak. Using left edge + 100.\n')
            %         else            
            %         end            
            figure(tuneData.chrg.figHandle); hold on;
            tlParams=scan.data.pulsegroups.params; % Has format tlcenter, tlDir.
            tlRng=[tlParams(1:2)+tlParams(3:4)*max(eps)*1e-3; tlParams(1:2)+tlParams(3:4)*min(eps)*1e-3];
            plot(tuneData.measPt(1)+tlRng(:,1)*1e-3,tuneData.measPt(2)+tlRng(:,2)*1e-3,'g-');
            
            tlptEps=tlParams(1:2)+tlParams(3:4)*tlPt*1e-3;
            plot(tuneData.measPt(1)+tlptEps(1)*1e-3,tuneData.measPt(2)+tlptEps(2)*1e-3,'kx');            
            for i=1:2
                for j = 1:length(scan.consts)
                    if strcmp(tuneData.xyChan{i},scan.consts(j).setchan)
                        out.measPt(i) = scan.consts(j).val;
                    end
                end
            end
        end
                
        function updateGroup(this,opts,config)
            %function updateGroup(this,config)
            % if config is a pulsegroup struct, uses that.
            % otherwise: make tlsweep entry in dictionary based on this.sweep,
            % update group based on this.search. add to awg memory. 
            % this will NOT update the feedback group.
            % Updating group correctly requires that triple points and default Offset are set correctly. 
            global tuneData;
            if ~exist('config','var'), config = []; end
            if ~exist('opts','var'), opts = ''; end
            if isstruct(config) %regular pulse group struct
                if isfield(config,'name') && ~strcmp(config.name,this.plsgrp)
                    error('pg.name = %s, tl plsgrp.name = %s\n',config.name,this.plsgrp);
                end
                plsdefgrp(config);
                awgadd(config.name);
                awgcntrl('on start wait err');
            elseif isempty(config) %make default group;
                if isopt(opts,'target'), this.target = this.location(end); end
                if tuneData.sepDir(1)>0  % TL meas point. Set up direction of lead. 
                    leadDir = [1 this.slope];
                    tlDir = [1 1./leadDir(2)];  tlDir=-tlDir/(norm(tlDir));
                else  % BR meas point
                    leadDir = -[1 this.slope];
                    tlDir = [1 1./leadDir(2)];  tlDir=tlDir/(norm(tlDir));
                end                
                leadDir = leadDir / norm(leadDir); % lead Direction, towards the right
                tlCenter=1e3*(1/2*(tuneData.chrg.trTriple(end,:)-tuneData.chrg.blTriple(end,:))-tuneData.chrg.defaultOffset) + leadDir*this.dist;
                tlSweep(1).type='wait';
                tlSweep(1).val = tlCenter + (this.target+this.sweep.range/2)*1e-3*tlDir;
                tlSweep(1).time=2e-3;
                tlSweep(2).type='ramp';
                tlSweep(2).val = tlCenter + (this.target-this.sweep.range/2)*1e-3*tlDir;
                tlSweep(2).time=this.sweep.time;
                
                dict = pdload(tuneData.activeSetName);
                dict.tlsweep = tlSweep;
                pdsave(tuneData.activeSetName, dict);
                
                %now make the group
                pg.name = this.plsgrp;
                pg.ctrl = 'pack loop';
                pg.dict={tuneData.activeSetName};
                pg.pulses = 7;
                pg.params= [tlCenter, tlDir, 0];
                pg.varpar = linspace(-this.search.range/2, this.search.range/2, this.search.points)' + this.target;
                pg.chan=[str2double(char(regexp(tuneData.xyChan{1},'\d+','match'))),str2double(char(regexp(tuneData.xyChan{2},'\d+','match')))];
                plsdefgrp(pg);
                awgadd(pg.name);
            else
                error('Must pass pulsegroup struct or empty to updateGroup')
            end
        end
    end
end