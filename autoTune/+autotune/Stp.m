classdef Stp < autotune.Op
    %Stp < autotune.Op %represents the search for the S-T+ transition
    %   see help for individual properties and methods.
    %   e.g. >> help autotune.Stp.plsGrp
    %   e.g. >> help autotune.Stp.Run
    
    properties
        plsGrp = 'STP_1_L'; % plsGrp for scan
        nLoop = 1000; % nloop for scan
        nRep = 10; % nrep for scan
        target = 1000; % target stp
        search = struct('range',1e3,'points',100,'time',200); %Parameters of stp search
        sweep = struct('range',300,'time',.3); % Params of sweep (for FB)
        fitFn = '@(p,x) p(1)+p(2)*exp(-(x-p(3)).^2/(2*p(4)^2)) + p(5)*x';
        subPlot = 6; % For plotting in tuneData.figHandle
        fineIndex = 1;
        filePat = 'stp';        
    end
    
    properties (SetAccess= {?autotune.Data, ?autotune.Op})
        location = 1000; % center where stp is found in uV of epsilon (Nrun x 1 double)
        width; % width of peak in uV (Nrun x 1 double)
        widtherr;
        foundSTP;
        slope;
    end
    
    methods
        function this = Stp
            global tuneData
            if strcmp(tuneData.activeSetName,'right')
                this.plsGrp ='STP_1_R';
            end
        end
        
        function out = getData(this,runNumber)
            out = struct();
            out.location = this.location(runNumber);
            out.width= this.width(runNumber);
        end
        
        function makeNewRun(this,runNumber) % what if you are redoing??
            if runNumber ~= length(this.location)+1 || runNumber ~= length(this.width)+1
                warning('runNumber not consistent with know chrg runs');
            end
            this.location(runNumber) = nan;
            this.width(runNumber) = nan;
            this.widtherr(runNumber) = nan;
            this.fineIndex = 1;
        end
        
        function run(this)
            global tuneData;
            if ~awgcntrl('ison'), awgcntrl('on start wait err'); end
            scan = fConfSeq(this.plsGrp,{'nloop',this.nLoop,'nrep',this.nRep, 'datachan',tuneData.dataChan,'opts','ampok'});
            scan = measAmp(scan);
            side = upper(tuneData.activeSetName(1)); 
            file = sprintf('%s/sm_%s%s_%04i_%03i',tuneData.dir, this.filePat,side,tuneData.runNumber,this.fineIndex);
            this.fineIndex = this.fineIndex+1;
            data = smrun(scan, file);
            this.ana('',data{1},scan);
        end
        
        function out=ana(this,opts,data,scan)
            global tuneData
            if ~exist('opts','var'), opts = ''; end
            runNumber = tuneData.runNumber;
            if ~exist('data','var'), data = []; end
            [data,out] = loadTunes(data,opts,this.filePat);                 
            if isempty(data), return; end
            if ~isfield(out,'scan'), out.scan = scan; end
            data = 1e3*nanmean(data,1);
            
            eps = out.scan.data.pulsegroups.varpar'*1e3;
            pf=polyfit(eps,data,1); % Remove offset, linear slope
            dataLin=smooth(data-pf(1)*eps - pf(2));
            
            ign=5; % Points at start and end to ignore.
            [maxSTP,maxSTPind]=max(dataLin(ign:end-ign));
            %       offset, STP height, STP loc,     STP width,   linslope
            beta0 = [pf(2),maxSTP,eps(maxSTPind+ign),range(eps)/8,pf(1)];
            mask = [0 1 1 1 0];
            
            axes(tuneData.axes(this.subPlot)); cla;
            params=fitwrap('plfit plinit samefig woff noplot',eps,data,beta0,this.fitFn,mask);
            [params,~,~,~,~,err]=fitwrap('plfit plinit samefig woff',eps,data,params,this.fitFn);
            stpPt=params(3)*this.slope*1e-6;
            a = gca; a.YTickLabelRotation=-30;
            a.XLim = [min(eps),max(eps)];
            title(sprintf('ST+: %3.1f; wdth %3.1f',params(3),abs(params(4))));
            
            figure(tuneData.chrg.figHandle); hold on;            
            stpRng= this.slope .* eps'*1e-6;
            % Plot range of scan in green, start point in white, end black. STP pt x            
            plot(tuneData.measPt(1)+stpRng(:,1),tuneData.measPt(2)+stpRng(:,2),'g-');            
            plot(tuneData.measPt(1)+stpRng(1,1),tuneData.measPt(2)+stpRng(1,2),'w.');
            plot(tuneData.measPt(1)+stpRng(end,1),tuneData.measPt(2)+stpRng(end,2),'k.');            
            plot(tuneData.measPt(1)+stpPt(1),tuneData.measPt(2)+stpPt(2),'kx'); 
            
            if ~out.anaData
                this.location(runNumber)=params(3);
                this.width(runNumber)=params(4);
                this.widtherr(runNumber) = err(1,4,1); %Check me!
                this.foundSTP = this.widtherr(end)<6;
            end
            for i=1:2
                for j = 1:length(out.scan.consts)
                    if strcmp(tuneData.xyChan{i},out.scan.consts(j).setchan)
                        out.measPt(i) = out.scan.consts(j).val;
                    end
                end
            end
        end
        
        function updateGroup(this,opts,config)
            % function updateGroup(this,config)
            % if config is a pulsegroup struct uses that.
            % otherwise: make stpsweep entry in dictionary based on this.sweep,
            % update group based on this.search.
            % add to awg memory. this will NOT update the feedback group.
            global tuneData
            
            if ~exist('config','var'), config = []; end
            if ~exist('opts','var'), opts = ''; end
            if isstruct(config) % Regular pulsegroup struct
                if isfield(config,'name') && ~strcmp(config.name,this.plsGrp)
                    error('pg.name = %s, stp plsGrp.name = %s\n',config.name,this.plsGrp);
                end
                plsdefgrp(config);
                awgadd(config.name);
                awgcntrl('on start wait err');
            elseif isempty(config) %make default group;
                if isopt(opts,'target'), this.target = this.location(end); end                
                
                % Make dictionary element for feedback. 
                dict = pdload(tuneData.activeSetName);
                exch = dict.exch.val(1:2);
                % Wait at 0,0 for 5 ns 
                stpSweep(1).type='@wait';
                stpSweep(1).time=0.005; % overshoot
                % Wait at point before STP for 1 ns. 
                stpSweep(end+1).type='wait';
                stpSweep(end).val = [exch (this.target-this.sweep.range/2)*1e-3];
                stpSweep(end).time = 1e-3;
                % Ramp to point past STP during search time. 
                stpSweep(end+1).type = 'ramp';
                stpSweep(end).val = [exch (this.target+this.sweep.range/2)*1e-3];
                stpSweep(end).time = this.sweep.time;
                this.slope = exch;
                dict.stpsweep = stpSweep;
                pdsave(tuneData.activeSetName, dict);
                
                % Make the group for stp scan 
                % Meas is hard coded at 1 us, may want to change. 
                pg.name = this.plsGrp;
                pg.ctrl='pack loop';
                pg.dict={tuneData.activeSetName};
                pg.pulses = 8;
                % Sit at exch magnitude given in varpar for time
                % this.search.time.
                pg.params = [this.search.time, 0];
                neps = this.search.points;
                % min step size.  Fixme to take awgdata.scale, sep dir into account
                minStep =6e5/142*2e-13;
                deps = round(this.search.range/(neps * minStep)) * minStep;
                pg.varpar= 1e-3*((this.target)+(-neps/2:neps/2-1)'*deps);
                pg.chan=[getNum(tuneData.xyChan{1}),getNum(tuneData.xyChan{2})];
                
                plsdefgrp(pg);
                awgadd(pg.name);
            else
                error('must past pulsegroup struct or empty to updateGroup')
            end
        end
    end
end