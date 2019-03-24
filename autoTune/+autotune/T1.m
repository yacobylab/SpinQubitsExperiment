classdef T1 < autotune.Op
   % Run scan to across junction to find tunnel coupling / temperature. 
    
    properties
        scan; % scan to take data
        plsGrp = 'dBzT1_2_L';         
        nrep = 50; % Number of repetitions of data
        npulse = 200; % Number of loops
        time = 100; % Time to wait at sep
        samplerate = 4e7; % DAQ sampling rate
        subPlot = NaN; %Uses own figure.
    end
    
    properties (SetAccess= {?autotune.Data, ?autotune.Op})
        fidelity; 
        STdiff; % Separation voltage 
        tMeas; % Opt meas time 
        t1; % t1 time 
    end
    
    methods                 
        function this = T1
            global tuneData
            if strcmp(tuneData.activeSetName,'right')
                this.plsGrp = 'dbZT1_2_R';                
            end
        end
        
        function out = getData(this,runNumber)
            out = struct();
            out.t1 = this.t1(runNumber);
            out.vis = this.vis(runNumber);
            out.tmeas = this.tmeas(runNumber);
        end
        
        function makeNewRun(this,runNumber)
            if runNumber ~= length(this.t1)+1 
                warning('runNumber not consistent with know chrg runs');
            end
            this.t1(end+1) = nan;
        end
        
        function run(this,opts)
            % If nograd is given as an option, will update the group. 
            global tuneData;
            if ~awgcntrl('ison'), awgcntrl('on start wait err');          end
            if ~exist('opts','var'), opts = ''; end
            file = sprintf('%s\\sm_t1%s_%04i', tuneData.dir, upper(tuneData.activeSetName(1)), tuneData.runNumber);            
            getchan = tuneData.dataChan;
            oversamp = this.samplerate*17e-6;                                     
            if isopt(opts,'nograd')
                if tuneData.t1.time < 50
                    tuneData.t1.time = 1000; tuneData.t1.updateGroup; 
                    awgcntrl('on start wait err'); 
                end
                conf=struct('nloop',1,'npulse',this.npulse,'nrep',this.nrep,'oversamp',oversamp,'hwsampler',this.samplerate,'datachan',getchan,'snglshot',0,'opts','nodisp nosave');                
            else
                if tuneData.t1.time > 50
                    tuneData.t1.time = 16; tuneData.t1.updateGroup; 
                    awgcntrl('on start wait err'); 
                end
                conf=struct('nloop',1,'npulse',this.npulse,'nrep',this.nrep,'oversamp',oversamp,'hwsampler',this.samplerate,'datachan',getchan,'snglshot',0,'opts','swfb nodisp nosave');
            end
            this.scan = fConfSeq(this.plsGrp,conf);
            data = smrun(this.scan,file); data = data{1};
            sleep; 
            this.ana(data,opts); 
        end
        
        function ana(this, data,opts)
            global tuneData 
            runNumber = tuneData.runNumber; 
            if ~exist('opts','var'), opts = ''; end
            if ~exist('data','var') || isempty(data) || ischar(data)
                if ~exist('data','var') || isempty(data) || ischar(data)
                    [data,scan]=loadAna('sm_t1*');
                else
                    [data,scan]=loadAna(data);
                end
                anaData=1;
            else
                anaData=0;
                scan = this.scan;
            end
            [fidelity,tMeas,STdiff,t1] = anaT1Meas(opts,scan,data); %#ok<*PROPLC>
            if ~anaData
                this.fidelity(runNumber)=fidelity;
                this.tMeas(runNumber)=tMeas;
                this.STdiff(runNumber)=STdiff;
                this.t1(runNumber)=t1;
            end
        end                
           
        function updateMeas(this)
            global tuneData 
            runNumber = tuneData.runNumber; 
            dict = pdload(tuneData.activeSetName); 
            dict.meas(1).time(1) = this.tMeas(runNumber)*1e6; 
            pdsave(tuneData.activeSetName,dict); 
        end
        
        function updateGroup(this)
            global tuneData; 
            pg.dict={tuneData.activeSetName};
            pg.ctrl='loop';
            pg.pulses = 20;
            % [startVoltage x /y, time)
            pg.params=[0 0 0];
            % Singlet pulse has 0 time at sep, Triplet this.time
            pg.varpar=[0; this.time];
            pg.chan=[str2double(char(regexp(tuneData.xyChan{1},'\d+','match'))),str2double(char(regexp(tuneData.xyChan{2},'\d+','match')))];
            pg.name = this.plsGrp; 
            plsdefgrp(pg);
            awgadd(pg.name); 
        end
    end    
end