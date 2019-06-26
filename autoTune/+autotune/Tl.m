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
        filePat = 'tl'; 
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
            scan = measAmp(scan);
            side = upper(tuneData.activeSetName(1)); 
            file = sprintf('%s/sm_%s%s_%04i_%03i',tuneData.dir,this.filePat, side,tuneData.runNumber,this.fineIndex);
            this.fineIndex = this.fineIndex+1; 
            data = smrun(scan, file);
            data = data{1};            
            this.ana('',data,scan);
        end
        
        function out=ana(this,opts,data,scan)
            global tuneData; global fbdata; 
            if ~exist('opts','var'), opts = ''; end            
            if ~exist('data','var'), data = []; end
            [data,out] = loadTunes(data,opts,this.filePat);                 
            if isempty(data), return; end
            if ~isfield(out,'scan'), out.scan = scan; end
            
            eps = out.scan.data.pulsegroups.varpar'; % tl scan sweeps epsilon value (along TL curve)
            data=1e3*mean(data,1); % Average multiple lines of data. 
            % Set up refval, used for histogramming.
            fbdata.refval(str2double(tuneData.dataChan(end))) = mean(data)/1000;  
            [~,maxTlInd]=max(data);
            epsMax=eps(maxTlInd); % TL is centered around where data largest. 
            %       1: offset, 2: upward slope, 3: upward center, 4: width both, 5: downward slope, 6: downward center
            beta0 = [min(data), range(data), epsMax-range(eps)/6, 150, range(data)/2, epsMax+range(eps)/6,-2e-6];
            axes(tuneData.axes(this.subPlot)); cla; 
            try
                params=fitwrap('plinit samefig woff',eps,data,beta0,this.fitFunc);
                [params,~,~,~,~,err]=fitwrap('plinit plfit samefig woff',eps,data,params,this.fitFunc);                
                tlPt = (params(3)+params(6))/2; % TL point is center point between 2 tanh fcns.
                tlWid = params(6)-params(3);
                err = err(1,4,1);                
                if ~out.anaData % If new data, add fit info to tuneData.
                    this.width(tuneData.runNumber) = tlWid;
                    this.location(tuneData.runNumber) = tlPt;
                    this.widtherr(tuneData.runNumber) = err;
                    this.foundTL = err < 100;
                end
                title(sprintf('TL %3.1f, wdth %3.1f',tlPt,err));
            catch
                warning('TL fit failed'); 
            end
            fitFunc = str2func(this.fitFunc); 
            plot(tlPt,fitFunc(params,tlPt),'x'); 
            plot(params(3),fitFunc(params,params(3)),'x'); 
            plot(params(6),fitFunc(params,params(6)),'x'); 
            a = gca; a.YTickLabelRotation=-30;
            a.XLim = [min(eps),max(eps)];            
            %            if params(6) > params(3)+ 800
            %               tlpt = eps(end)-100;
            %              fprintf('Broad peak. Using left edge + 100.\n')
            %         else            
            %         end            
            tlParams=out.scan.data.pulsegroups.params; % Has format tlcenter, tlDir.
            tlRng=[tlParams(1:2)+tlParams(3:4)*max(eps)*1e-3; tlParams(1:2)+tlParams(3:4)*min(eps)*1e-3];
            tlptEps=tlParams(1:2)+tlParams(3:4)*tlPt*1e-3;
            
            % Plot range of scan in green, start point in white, end black. STP pt x            
            figure(tuneData.chrg.figHandle); hold on;
            plot(tuneData.measPt(1)+tlRng(:,1)*1e-3,tuneData.measPt(2)+tlRng(:,2)*1e-3,'g-');
            plot(tuneData.measPt(1)+tlRng(1,1)*1e-3,tuneData.measPt(2)+tlRng(1,2)*1e-3,'w.');
            plot(tuneData.measPt(1)+tlRng(end,1)*1e-3,tuneData.measPt(2)+tlRng(end,2)*1e-3,'k.');                        
            plot(tuneData.measPt(1)+tlptEps(1)*1e-3,tuneData.measPt(2)+tlptEps(2)*1e-3,'kx');            
            for i=1:2
                for j = 1:length(out.scan.consts)
                    if strcmp(tuneData.xyChan{i},out.scan.consts(j).setchan)
                        out.measPt(i) = out.scan.consts(j).val;
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
                    tlDir = [-1 1./leadDir(2)];  tlDir=-tlDir/(norm(tlDir));
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
                
                % TL search group 
                % Length hard coded at 4 us, meas at 1 us.                 
                pg.name = this.plsgrp;
                pg.ctrl = 'pack loop';
                pg.dict={tuneData.activeSetName};
                pg.pulses = 7;
                pg.params= [tlCenter, tlDir, 0];
                pg.varpar = linspace(-this.search.range/2, this.search.range/2, this.search.points)' + this.target;
                pg.chan=[getNum(tuneData.xyChan{1}),getNum(tuneData.xyChan{2})];
                plsdefgrp(pg);
                awgadd(pg.name);
            else
                error('Must pass pulsegroup struct or empty to updateGroup')
            end
        end
    end
end