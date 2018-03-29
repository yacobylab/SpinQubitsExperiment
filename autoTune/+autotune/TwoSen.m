classdef TwoSen < autotune.Op
    %Zoom < autotune.Op % represents a zoom scan to find Pauli Blockade
    %   see help for individual properties and methods.
    %   e.g. >> help autotune.Zoom.plsGrp
    %   e.g. >> help autotune.Zoom.Run
    
    properties
        scan; % scan for taking data
        subPlot = nan; 
        sensThresh = 5; 
        peakSep = 0.02;
        sensorInd
        offset = 0.1; 
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
        end
        
        function out = getData(this,runNumber)
        end
        
        function makeNewRun(this,runNumber)
        end
        
        function run(this)                                                
            % TwoD Scan.             
            
            global tuneData;                      
            ison = awgcntrl('ison')>0; % 0.5 = waiting for trigger
            if any(~ison)
                disp('AWG is off.  <a href="matlab:awgcntrl(''on start wait err'');dbcont">Turn it on?</a> <a href="matlab:disp(''exiting...'');dbquit">Exit?</a> ')
                keyboard
            end            
            
            file = sprintf('%s/sm_SD%s_%04i', tuneData.dir, upper(tuneData.activeSetName(1)), tuneData.runNumber);
            d=smrun(tuneData.twoSen.scan,file);
            diffData = diff(d{1},1,2);  
           	[maxCol,ind2]=max(abs(diffData));
            [~,ind1] = max(abs(maxCol)); 
            ind2 = ind2(ind1);
            loop1Vals = scanRng(this.scan,1); 
            loop2Vals = scanRng(this.scan,2); 
            smset(this.scan.loops(1).setchan{1},loop1Vals(ind1)); 
            smset(this.scan.loops(2).setchan{1},loop2Vals(ind2));                                     
        end                
    end
end