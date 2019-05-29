classdef Zoom < autotune.Op
    %Zoom < autotune.Op % represents a zoom scan to find Pauli Blockade
    %   see help for individual properties and methods.
    %   e.g. >> help autotune.Zoom.plsGrp
    %   e.g. >> help autotune.Zoom.Run
    
    properties
        scan; % scan for taking data
        plsGrp = 'zoom_2_L'; %string, pulsegroup fro file
        offSet = [0, 0]; % offset from chg scan
        subPlot = NaN; %zoom scan gets its own plot :)
        res=[]
        rng = [4e-3, 4e-3]; % Range of zoom scan.
        pulseScan;
    end
    
    properties (SetAccess= {?autotune.Data, ?autotune.Op})
        measPt; % Nrun x 2 double
    end
    
    properties (Constant)
        figHandle = 2; % for plotting zoom scan
    end
    
    methods
        function this = Zoom
            global tuneData
            if strcmp(tuneData.activeSetName,'right')
                this.plsGrp ='zoom_1_R';
            end
        end
        
        function out = getData(this,runNumber)
            out = struct('measPt',this.measPt(runNumber));
        end
        
        function makeNewRun(this,runNumber)
            if runNumber ~= size(this.measPt,1)+1
                warning('runNumber not consistent with know chrg runs');
            end
            this.measPt(runNumber,:) = [nan,nan];
        end
        
        function run(this, opts)
            % run(this, opts)
            % opts: narrow: scan run with 4 mV range.
            % Run a scan with rand, run scan with rand + reload different between scans should be readout region
            global tuneData;
            runNumber = tuneData.runNumber;
            if ~exist('opts','var'), opts = ''; end
            awgcntrl('start on amp');
            allData=cell(2,1);
            zoomOffset=tuneData.chrg.defaultOffset+1/2*(tuneData.chrg.blTriple(runNumber,:)+tuneData.chrg.trTriple(runNumber,:)); %center of scan
            if isopt(opts,'narrow')
                this.scan = smscanpar(tuneData.chrg.scan,zoomOffset,[2,2]*1e-3, this.res);
            elseif isopt(opts,'wide')
                this.scan = smscanpar(tuneData.chrg.scan,zoomOffset,[8,8]*1e-3, this.res);
            else
                this.scan = smscanpar(tuneData.chrg.scan,zoomOffset,this.rng, this.res);
            end
            clearMask(tuneData.dataChan);
            for i = 1:2
                file=sprintf('%s/sm_zoom%s_%04i_%d', tuneData.dir, upper(tuneData.activeSetName(1)), runNumber ,i);
                this.scan.consts(2).val = awgseqind(this.plsGrp)+ i;
                data = smrun(this.scan, file);
                allData{i} = data{1};
            end
            this.ana(opts,allData);
        end
        
        function out = ana(this,opts,data)
            % opts: 
            %   noset(don't set the the measurement point) 
            %   setOffset (center zoom along junc)
            %   man (choose zoom point manually)
            %   last: use most recent dataset. 
            % data can be: filename, number or empty (bring up dialog). 
            
            global tuneData
            if ~exist('opts','var'), opts = ''; end            
            % Check if loading old data or analyzing new scan.
            if ~exist('data','var') || isempty(data) || ischar(data) || numel(data)==1 
                if (~exist('data','var') || isempty(data)) && ~isopt(opts, 'last')
                    % Load file through file ui
                    [d1,scan,fname]=loadAna('sm_zoom*1.mat');
                    d2=load([fname(1:end-5) '2.mat'],'data');
                    data={d1 d2.data{1}};
                    runNumber = str2double(fname(end-9:end-6)); 
                elseif exist('data','var') && ~isempty(data) && ischar(data)
                    % File name given
                    [d1,scan,~,time]=loadAna(data,opts);
                    d2 = load([tuneData.dir '\' data(1:end-5),'2.mat'],'data');
                    data ={d1 d2.data(1)};                   
                    out.time = time;
                    runNumber = str2double(data(end-9:end-6)); 
                else
                    % Use number
                    side = upper(tuneData.activeSetName(1));
                    if isopt(opts,'last'), data = tuneData.runNumber; end
                    runNumber = data; 
                    fileName = sprintf('sm_zoom%s_%04.f_1.mat',side,data);
                    [d1,scan,~,time]=loadAna(fileName);
                    d2=loadAna([fileName(1:end-5) '2.mat']);
                    data={d1 d2};
                    if isempty(d1) || isempty(d2)
                        out = struct;
                        return
                    else
                        out.time = time;
                    end
                end
            else
                scan = this.scan; %#ok<*PROPLC> % new data
                runNumber = tuneData.runNumber;
            end
            out.scan = scan; 
            dataDiff = data{2}- data{1};
            grad = gradient(data{1});
            dataDiff(grad > 5 * std(grad(:))) = NaN; %remove outliers.
            
            % Plot
            f=figure(this.figHandle); clf;
            f.Name='Zoom Analysis';
            xvals = scanRng(scan,1); yvals = scanRng(scan,2);
            imagesc(xvals,yvals, dataDiff);
            xlabel(scan.loops(1).setchan{1}); ylabel(scan.loops(2).setchan{1});
            set(gca,'YDir','Normal');
            axis image; hold on;
            % Analyze zoom point
            if ~isopt(opts,'noset') || tuneData.runNumber ~= runNumber
                if isopt(opts,'man')
                    this.onlyPlot;
                else
                    this.anaZoomScan(xvals,yvals,dataDiff);
                end
                if isopt(opts,'setOffset') % Let's you select the zoom point, then assumes that the dist btwn junc and pt is correct, but centers it along junc.
                    juncSlp = tuneData.chrg.trTriple(runNumber,:) - tuneData.chrg.blTriple(runNumber,:); juncSlp = juncSlp/juncSlp(1);
                    epsSlp = [1 -1./juncSlp(2)]; epsSlp = epsSlp/norm(epsSlp);
                    juncCen = (tuneData.chrg.trTriple(runNumber,:) + tuneData.chrg.blTriple(runNumber,:))/2;
                    offSet = tuneData.measPt - juncCen;
                    offSetEps = dot(epsSlp,offSet);
                    tuneData.measPt = juncCen + offSetEps*epsSlp;
                    tuneData.chrg.defaultOffset = offSetEps;
                end
                [~,nearX] = min(abs(xvals-tuneData.measPt(1)));
                [~,nearY] = min(abs(yvals-tuneData.measPt(2)));
                if dataDiff(nearY,nearX)>0
                    flipDAQ(tuneData.dataChan);
                    fprintf('Changing sign of readout card \n');
                end
                tuneData.zoom.measPt(runNumber,:) = tuneData.measPt;
            end
            plot(tuneData.chrg.trTriple(runNumber,1),tuneData.chrg.trTriple(runNumber,2),'r.','MarkerSize',10)
            plot(tuneData.chrg.blTriple(runNumber,1),tuneData.chrg.blTriple(runNumber,2),'r.','MarkerSize',10)
            plot(tuneData.measPt(1), tuneData.measPt(2), 'k.', 'MarkerSize', 20);
            figure(tuneData.chrg.figHandle);
            plot(tuneData.measPt(1), tuneData.measPt(2), 'k.', 'MarkerSize', 20);            
            figure(this.figHandle);
        end
        
        function updateGroup(this,config)
            global tuneData
            % function updateGroup(this,config)
            % if config is a pulsegroup struct uses that.
            % otherwise: make stpsweep entry in dictionary based on this.sweep, update group based on this.search.
            % add to awg memory. this will NOT update the feedback group.
            if ~exist('config','var'), config = []; end
            if isstruct(config) %regular pulsegroup struct
                if isfield(config,'name') && ~strcmp(config.name,this.plsGrp)
                    error('pg.name = %s, stp plsGrp.name = %s\n',config.name,this.plsGrp);
                end
                try
                    plsupdate(config);
                    awgadd(config.name);
                    awgcntrl('on start wait err');
                catch
                    fprintf('Problem making groups %s. quitting...\n',config.name);
                end
            elseif isempty(config) %make default group;
                pg.name = this.plsGrp;
                pg.ctrl='pack loop';
                pg.dict={tuneData.activeSetName};
                pg.pulses = 23;
                pg.params = 0;
                pg.varpar = [0 .5]';
                pg.nrep=0;
                pg.ctrl = '';
                pg.chan=[str2double(char(regexp(tuneData.xyChan{1},'\d+','match'))),str2double(char(regexp(tuneData.xyChan{2},'\d+','match')))];
                
                plsdefgrp(pg);
                awgadd(pg.name);
            else
                error('must past pulsegroup struct or empty to updateGroup')
            end
        end
        
        function measPt = pulsed(this,opts)
            global tuneData;
            if ~exist('opts','var'), opts = ''; end
            rng = 4;
            if ~isopt(opts,'run')
                pg.ctrl='loop pack';
                pg.dict={tuneData.activeSetName};
                pg.pulses = 111;
                pg.varpar = linspace(-rng/2,rng/2,100)';
                pg.chan=[str2double(char(regexp(tuneData.xyChan{1},'\d+','match'))),str2double(char(regexp(tuneData.xyChan{2},'\d+','match')))];
                pg.nrep = 1;
                yvals = linspace(-rng/2,rng/2,30);
                loadTime = [0,0.5];
                for j = 1:2
                    for i = 1:length(yvals)
                        pg.name = sprintf('pulsedZoomL%d_%d',j,i);
                        pg.params = [loadTime(j) yvals(i) 0];
                        plsdefgrp(pg);
                        zoomGrp{i}=pg.name;
                    end
                    awgadd(zoomGrp);
                    this.pulseScan{j} = zoomGrp;
                end
            end
            awgcntrl('on start wait err')
            for j = 1:2
                scan = fConfSeq(this.pulseScan{j},{'nloop',100,'nrep',1, 'datachan',tuneData.dataChan,'opts','ampok'});
                scan.consts(end+1).setchan=tuneData.xyChan{1};
                scan.consts(end).val=tuneData.measPt(1);
                scan.consts(end+1).setchan=tuneData.xyChan{2};
                scan.consts(end).val=tuneData.measPt(2);
                scan.data.measPt = tuneData.measPt;
                d = smrun(scan,smnext(sprintf('pulsedZoom%dL',j)));
                data{j} = d{1};
            end
            figure(tuneData.zoom.figHandle); clf;
            dataDiff = data{2}- data{1};
            imagesc([-rng/2,rng/2],[-rng/2,rng/2], squeeze(dataDiff));
            set(gca,'YDir','Normal')
            measPt = ginput(1);
            rep = input('Update measurement point and load?','s');
            if isopt(rep,'y')
                l = pdload('left');
                l.reload.val = l.reload.val - measPt;
                pdsave('left',l);
                tuneData.loadPos.updateGroup;
                tuneData.measPt = measPt*1e-3+scan.data.measPt;
            end
            if isfield(scan,'data') && isfield(scan.data,'measPt')
                measPt = measPt + scan.data.measPt;
            end
        end
    end
    
    methods (Static)
        function  onlyPlot
            global tuneData
            tuneData.measPt = ginput(1);
        end
        
        function anaZoomScan(xvals,yvals,dataDiff)
            global tuneData
            center = tuneData.zoom.imageRecZoom(dataDiff);
            measPt = [xvals(round(center(1))),yvals(round(center(2)))];
            tuneData.measPt = measPt;
        end
        
        function center = imageRecZoom(data)
            debug = 0;
            dimg=mat2gray(data);
            % Increase contrast. White is 1, Zoom = black.
            % If already has large contrast, dimg has large std.
            % If the median darker than mean, likely very low contrast, since the zoom region is << 1/2 the device size.
            % want to increase the contrast around the zoom region.
            typ = [0.4 0.9];
            cnt = mean(dimg(:));
            cntId=(cnt - 0.7) +typ;
            cntZm = median(dimg(:)) - mean(dimg(:));
            s = sign(cntZm); n=1;
            cntId(1) = cntId(1) - n * s * std(dimg(:));
            cntId(2) = cntId(2) + n * s * std(dimg(:));
            cntId(1) = max(0,cntId(1)); cntId(2) = min(1,cntId(2));
            improcd=imadjust(dimg,cntId,[]);
            bw = im2bw(improcd,0.4); %threshold data
            se = strel('disk', 1);
            bwCl = imclose(bw,se);
            bwOp = imopen(bwCl,se);
            if debug
                figure(11); clf;
                subplot(2,3,1)
                imshow
            end
            zoomData = 1 - dimg;  % for weighting, want zoom region -> 1
            for i = 1:2
                props = regionprops(bwOp,zoomData,'area','centroid','weightedcentroid','PixelList','PixelIdxList');
                [area,n] = max([props.Area]);
                if area < numel(dimg) / 2
                    break
                end
                bwOp = ~bwOp;
            end
            center = props(n).WeightedCentroid;
        end
    end
end