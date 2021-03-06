classdef Sensor < autotune.Op
    %Zoom < autotune.Op % represents a zoom scan to find Pauli Blockade
    %   see help for individual properties and methods.
    %   e.g. >> help autotune.Zoom.plsGrp
    %   e.g. >> help autotune.Zoom.Run
    
    properties
        scan; % scan for taking data
        subPlot = nan; %zoom scan gets its own plot :)
        sensThresh = 5; 
        peakSep = 0.02;
        sensorInd
        offset = 0.1; 
        failNum =0; % Nrun x 2 double
        fineIndex=1;
    end
    
    properties (SetAccess= {?autotune.Data, ?autotune.Op})        
        newVal;        
        fail;
        sens;
    end
    
    properties (Constant)
        figHandle = 2; %fot plotting zoom scan
    end
    
    methods
        function this = Sensor    
            global tuneData
            scan = defScan('chrg',tuneData.activeSetName); 
            scan.loops(1).npoints = 200; 
            scan.loops(1).ramptime = -10e-3; 
            scan.loops(1).rng = [-.42 -.57];             
            scan.loops(2).npoints = 2;
            if strcmp(tuneData.activeSetName,'right')
                scan.loops(1).setchan = {'SD4top'};
            else
                scan.loops(1).setchan = {'SD1top'};
            end
            this.scan = scan; 
        end
        
        function out = getData(this,runNumber)
        end
        
        function makeNewRun(this,runNumber)
            this.fineIndex = 1;
        end
        
        function run(this,opts)                                                
            % Just 1D scan for now. Selects most sensitive point.
            % If most sensitive point on either end, shifts scan by
            % tuneData.sensor.offset. If most sensitive point < 5, fails.
            % If greatest sensitivity negative, changes DAQ multiplier by -1.
            % If fails to find optimum 5 times, quits.
            % In future, can add auto 2D scan when this happens.            
            
            global tuneData;  
            if ~exist('opts','var'), opts = ''; end
            if this.failNum > 5
                fprintf('Sensor failing to converge. Please retune manually \n')                
                this.failNum = 0; 
            end
            if length(this.sensorInd) == length(tuneData.runNumber) 
                this.sensorInd(tuneData.runNumber)= this.sensorInd(tuneData.runNumber)+1; 
            else 
                this.sensorInd(tuneData.runNumber) = 1; 
            end
            file = sprintf('%s/sm_sensor%s_%04i_%03i', tuneData.dir,...
                upper(tuneData.activeSetName(1)), tuneData.runNumber,this.fineIndex);
            oldVal=cell2mat(smget(tuneData.sensor.scan.loops(1).setchan));
            d=smrun(tuneData.sensor.scan,file);
            
            %Find slope of data
            xvals=scanRng(this.scan,1); 
            data = nanmean(d{1}); 
            diffData = diff(data)./(xvals(2)-xvals(1));
            xDiff = (xvals(1:end-1)+xvals(2:end))/2;                        
            
            % Find maximum slope of new point, old point.
            [maxD,maxInd]=max(abs(diffData));
            newVal = xDiff(maxInd);
            if abs(newVal-oldVal)>this.peakSep
                nearInds = find(abs(xDiff-oldVal)<this.peakSep/2);
                if (maxD-max(abs(diffData(nearInds))))<5
                    [maxD,maxInd]=max(abs(diffData(nearInds)));
                    maxInd = nearInds(maxInd); 
                    newVal = xDiff(maxInd); %#ok<*PROPLC>
                end
            end
            [~,oldInd] = min(abs(oldVal-xDiff));
            oldSens = diffData(oldInd);
            if oldVal < min(xDiff) || oldVal > max(xDiff)
                oldSens = NaN;
            end
            
            % If greatest sensitivity negative, change DAQ multiplier.
            if diffData(maxInd) < 0
                global smdata; %#ok<*TLEV>
                DAQchan = smchanlookup(tuneData.sensor.scan.loops(2).getchan);
                smdata.channels(DAQchan).rangeramp(4) = -1*smdata.channels(DAQchan).rangeramp(4);
                data = -data;
            end
            f=figure(71); clf;
            f.Name = 'Sensor scan'; 
            subplot(2,1,2); plot(xvals,data,'.-');
            xlabel(tuneData.sensor.scan.loops(1).setchan); ylabel('Signal');
            
            subplot(2,1,1); plot(xDiff,diffData,'.-');
            ylabel('Derivative');
            
            % Plot red mark on most sensitive point. 
            hold on; plot(newVal,diffData(maxInd),'.','MarkerSize',10);
            fitAxis;
            subplot(2,1,2); hold on;
            plot(newVal,data(maxInd),'.','MarkerSize',10); fitAxis;
            
            % Set channel to new value
            this.newVal = newVal;
            smset(this.scan.loops(1).setchan,newVal)
            this.sens(tuneData.runNumber) = maxD;
            fprintf('Max sensor diff of %4.4f at %4.4f \n',maxD,newVal);
            fprintf('Shift sensor by %4.4f mV. Sensitivity had been %4.4f  \n',1e3*(newVal-oldVal),oldSens)
            fprintf('Mean sensor value will be %4.4f \n',(data(maxInd)+data(maxInd+1))/2)
            this.fail = 0;
            this.fineIndex = this.fineIndex+1; 
            %Check that points aren't on the end and that sensitivity meets threshold.
            if maxInd == length(diffData)
                fprintf('Most sensitive at edge point, shifting scan \n');
                this.scan.loops(1).rng = this.scan.loops.rng + this.offset;
                this.fail = 1;
                this.failNum =this.failNum+1;
                autotune('sensor')
            elseif maxInd == 1
                fprintf('Most sensitive at edge point, shifting scan \n');
                this.scan.loops(1).rng = this.scan.loops(1).rng - this.offset;
                this.fail = 1;
                this.failNum = this.failNum+1;
                tuneData.sensor.run
            elseif maxD < this.sensThresh
                fprintf('Low sensitivity. Change range or run 2D scan \n')
                this.fail = 1;
                this.failNum = this.failNum+1;
            else
                this.failNum = 0;
            end
         sleep('fast');
         if isopt(opts,'fine')
             scan= this.scan; 
             this.scan.loops(1).rng = [newVal-0.005,newVal + 0.005]; 
             this.scan.loops(1).npoints = floor(this.scan.loops(1).npoints/3); 
             this.scan.loops(1).ramptime = -0.025; 
             try
                this.run; 
             end
             this.scan = scan; 
         end
        end                
    end
end