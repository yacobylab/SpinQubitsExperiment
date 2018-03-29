classdef Line < autotune.Op
   % Run scan to across junction to find tunnel coupling / temperature. 
    
    properties
        scan; % scan to take data
        fitFn = '@(p, x)p(1)+p(2)*(x)-p(4)*(1+tanh((x-p(3))./p(5)))/2'; %fit function for fitting
        perpDist = 0.55; %distance parallel to jcn to scan. 0 = blTriple, 1 = trTriple. 0.5 = middle
        subPlot = 3; % for plotting in tuneData.figHandle
        nrep = 5; 
    end
    
    properties (SetAccess= {?autotune.Data, ?autotune.Op})
       width; % double, size Nrun 
    end
    
    methods         
        function this = Line            
            this.makeNewScan;             
        end
        
        function out = getData(this,runNumber)
            out = struct();
            out.width = this.width(runNumber);
        end
        
        function makeNewRun(this,runNumber)
            if runNumber ~= length(this.width)+1 
                warning('runNumber not consistent with know chrg runs');
            end
            this.width(end+1) = nan;
        end
        
        function run(this)
            global tuneData;
            file = sprintf('%s/sm_line%s_%04i', tuneData.dir, upper(tuneData.activeSetName(1)),tuneData.runNumber);
            this.scan.consts = tuneData.chrg.scan.consts;
            smset({this.scan.consts.setchan}, [this.scan.consts.val]);
            % Permute this away from the center a bit to avoid the mousebite.
            chrg = [tuneData.chrg.blTriple(tuneData.runNumber,:),tuneData.chrg.trTriple(tuneData.runNumber,:)];
            cntr = [chrg(1)*this.perpDist+chrg(3)*(1-this.perpDist), chrg(2)*this.perpDist+chrg(4)*(1-this.perpDist)];
            this.scan.loops(1).trafofn(1).args = {cntr(1)};
            this.scan.loops(1).trafofn(2).args = {cntr(2)};
            data = smrun(this.scan, file); data = data{1}; 
            this.ana('',data);
        end
        
        function ana(this,opts,data)  
            global tuneData;   
            if ~exist('opts','var') 
                 opts = '';
            end 
            runNumber = tuneData.runNumber; 
            if ~exist('data','var') || isempty(data) || ischar(data) || numel(data)==1
                if (~exist('data','var') || isempty(data)) && ~isopt(opts,'last')
                    [data,scan]=loadAna('sm_line*');
                elseif exist('data','var') && ~isempty(data) && ischar(data)
                    [data,scan]=loadAna(data);
                else
                    side = upper(tuneData.activeSetName(1));
                    if isopt(opts,'last')
                        data = tuneData.runNumber;
                    end
                    fileName = sprintf('sm_line%s_%04.f.mat',side,data);
                    [data,scan]=loadAna(fileName);
                    if isempty(data)
                        return
                    end
                end
            else
                scan = this.scan;
            end
            eps = linspace(scan.loops(1).rng(1), scan.loops(1).rng(2), scan.loops(1).npoints);  
            data = mean(data);%FIXME
         
            slp=median(smooth(diff(data)))/(eps(2)-eps(1)); % do a linear fit for initial guess
            mn=mean(data);     
            ign = 5; 
            dataFlt=data-slp(1)*eps; % Look at the skew of the derivative to guess a sign for the atan.
            datadiff = abs(smooth(diff(dataFlt)));            
            [~, mi] = max(datadiff(ign:end-ign)); % Find the step
            stepCent=eps(mi+ign);
            
            if(mean(dataFlt(1:mi+ign)) < mean(dataFlt(mi+ign:end)))
                stepSlp=-range(dataFlt);
            else
                stepSlp=range(dataFlt);
            end
               
            beta0 = [mn,slp, stepCent, stepSlp, range(eps)/16.0];
            axes(tuneData.axes(this.subPlot)); 
            fitPars = fitwrap('plfit woff fine samefig', eps, data, beta0,this.fitFn,[1 1 0 0 0]);
            fitPars = fitwrap('plinit plfit woff fine samefig', eps, data, fitPars,this.fitFn,[1 1 1 1 1]);
            title(sprintf('Line: %g mV',fitPars(5)*1e3));
            a = gca; a.YTickLabelRotation=-30;
            a.YLabel.Position(1) = a.XLim(1) - range(a.XLim)/14;
            a.XLabel.Position(2) = a.YLim(1) - range(a.YLim)/7;
            fprintf('Line: %g mV\n',fitPars(5)*1e3);
            this.width(runNumber) = fitPars(5);
        end
        
        function makeNewScan(this)
            global tuneData
            scan = defScan;  %#ok<*PROP,*PROPLC>
            scan.loops(2).npoints = this.nrep;
            scan.loops(1).ramptime = -5e-3;
            scan.loops(1).rng = [-1e-3 1e-3];
            scan.loops(1).npoints = 100;
            scan.consts(2).val = awgseqind(tuneData.chrg.pls);
            scan.loops(2).getchan = tuneData.dataChan;
            scan.loops(1).setchan = tuneData.xyChan;
            scan.loops(1).trafofn(1).fn = @(x,y,cntr)x(1)+cntr;
            scan.loops(1).trafofn(2).fn = @(x,y,cntr)-x(1)+cntr;
            scan.loops(1).stream = 1; 
            this.scan = scan;
        end
    end    
end