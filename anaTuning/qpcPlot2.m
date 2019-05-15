function [file,out] = qpcPlot2(opts,file,config)
% Plot a bunch of the hysteresis curves.
% function file = qpcPlot2(opts,file,config)
% opts: 
%   offset : shifts the xval by amount config
%   hyst: Find min pinch off point, then plot how far each ramp rate is from there.  
%   mean
%   num

if ~exist('opts','var'),     opts = ''; end
if ~exist('file','var') || isempty(file)
    [file,fpath]=getFiles('sm_*.mat');
else
    fpath = [pwd '\'];
end
if file{1}==0,      return;    end
str='';

nDim = 1;
for i = 1:length(file)
    d{i} = load([fpath, file{i}]);
    nDim = max(nDim,length(d{i}.scan.loops));
    if ~isfield(d{i}.scan,'data') || ~isfield(d{i}.scan,'setchan')
        try
            d{i}.scan.data.setchan = getMultiGateScan(d{i}.scan);
        catch 
            warning('Can''t find channel names'); 
        end
    end
    if ~isfield(d{i},'data')
    	fprintf('%s empty file \n',file{i})
    end
end
if nDim == 3 && ~isopt(opts,'mean')
    zeroVal = nan(length(d),d{1}.scan.loops(2).npoints,d{1}.scan.loops(3).npoints);
else
    zeroVal = nan(length(d),d{1}.scan.loops(2).npoints);
end % preallocate the zerovals. 
condQuant = zeroVal;

fPlot = makeFigure(12); nrow = 2; 
ga= tight_subplot(nrow,3,[0.05 0.03],0.04,0.04);
fHyst = makeFigure(100);
fignums = [12 100]; 
gaa{1} = tight_subplot(2*nrow,3,[0.05 0.03],0.04,0.04);
if isopt(opts,'invis')
    fHyst.Visible = 'off';
    fPlot.Visible = 'off';
else
    fHyst.Visible = 'on';
    fPlot.Visible = 'on';
end
for i = 1:length(d)
    if ~isfield(d{i},'data')
        continue
    end
    if isfield(d{i}.scan,'data') && isfield(d{i}.scan.data,'rangeramp')        
        [~,d{i}.scan.loops(1).rng] = smramprate(d{i}.scan,'save');
    end
    xvals=scanRng(d{i}.scan);
    rampRate(i) = diff(d{i}.scan.loops(1).rng)/d{i}.scan.loops(1).npoints/abs(d{i}.scan.loops(1).ramptime);
    minVal(i) = min(d{i}.scan.loops(1).rng);    
    if isopt(opts, 'offset') && rampRate(i) > 0 % helps to account for bias cooling 
        xvals = xvals + config;
    end
    if ismatrix(d{i}.data{1})
        data = d{i}.data{1};
        if length(d{i}.data)>1
            data2 = d{i}.data{2};
        else, data2 = [];
        end
    else%if isopt(opts,'mean')
        data = squeeze(mean(d{i}.data{1},1));
        if length(d)>1
            data2 = squeeze(mean(d{i}.data{2},1));
        else, data2 = [];
        end
    end
    for j = 1:d{i}.scan.loops(2).npoints
        if isopt(opts,'hyst')
            if ~isempty(data2) % if we've recorded the quantum of conductance, can use that.
                if any(data2(j,:)<0)
                    [~,mi]= min(abs(data2(j,:)));
                    condQuant(i,j) = xvals(mi); % row is which scan, j is which gate.
                end
            else
                if any(data(j,:)<2e-9)
                    [~,mi]= min(abs(abs(data(j,:))-2e-9));
                    condQuant(i,j) = xvals(mi);
                end
            end
            if any(data(j,:)<1e-10)
                [~,mi]= min(abs(abs(data(j,:))-1e-10));
                zeroVal(i,j) = xvals(mi);
            end
            dispName = sprintf('%1.1f, %1.2f, %s', 1e3*rampRate(i),min(d{i}.scan.loops(1).rng),file{i}(end-7:end-4));
        elseif isopt(opts,'num')
            dispName = file{i}(end-7:end-4);
        else, dispName = '';
        end
                a=ga(mod(j-1,nrow*3)+1);
        plot(a,xvals,abs(data(j,:)),'DisplayName',dispName); hold(a,'on');
    end
end

minValList = unique(minVal);
colors = a.ColorOrder;
for j = 1:d{1}.scan.loops(2).npoints
    cQtmp = squeeze(condQuant(:,j,:)); % everything from a single gateset
    mV(j) = nanmean(cQtmp(:));         % mean condquant
    cQZtmp = squeeze(zeroVal(:,j,:));
    mVZ(j) = nanmean(cQZtmp(:));
end
out = struct('condQuant',condQuant,'zeroVal',zeroVal,'chans',{d{1}.scan.data.setchan}); 
for j = 1:d{i}.scan.loops(2).npoints
    if isfield(d{1}.scan,'data')
        xlab=d{1}.scan.data.setchan{j};
    else
        xlab = {'dummy'};
    end
    if ischar(xlab)
        xlab = {xlab};
    end
    a=ga(mod(j-1,nrow*3)+1);
    a.XLabel.String =xlab;
    if isopt(opts,'noppt')
        l = clickableLegend(a,'show');
        l.Box = 'off'; 
    else
        l=legend(a,'show');
    end
    l.FontSize = 8;
    l.Interpreter='none'; l.Location = 'northwest';
    
    a=gaa{1}(mod(j-1,nrow*3)+1); 
    for i = 1:length(d)
        ind(i) = find(minVal(i)==minValList);
        plot(a,1000*rampRate(i), squeeze(condQuant(i,j,:))-mV(j),'*','Color',colors(ind(i),:)); hold(a,'on');
    end
    a.Title.String = [xlab{1} ' ' num2str(mV(j),'%3.3f')];
    a.XLabel.Position(2) = a.YLim(1) - range(a.YLim)/10.5;
    a.XLabel.String =' Ramprate (mV /s)';
    a.XBaseline.Visible = 'on'; a.YBaseline.Visible = 'on';     
    a=gaa{1}(mod(j-1,nrow*3)+1+nrow*3);
    for i = 1:length(d)
        plot(a,1000*rampRate(i), squeeze(zeroVal(i,j,:))-mVZ(j),'*','Color',colors(ind(i),:)); hold(a,'on');
    end
    a.Title.String = [xlab{1} ' ' num2str(mVZ(j),'%3.3f')];
    a.XLabel.Position(2) = a.YLim(1) - range(a.YLim)/10.5;
    a.XLabel.String =' Ramprate (mV /s)';
    a.XBaseline.Visible = 'on'; a.YBaseline.Visible = 'on';     
    str = [str sprintf('%s : \t %3.3f   %3.3f   %3.3f  \n',xlab{1},mV(j),mVZ(j),mVZ(j)-mV(j))]; 
end

a=gaa{1}(6);
for i = 1:length(minValList)
    text(a,0.25,0.23+i*0.086,num2str(minValList(i)),'Color',colors(i,:))
end

formatFig(fignums,'qpc2'); 
if ~isopt(opts,'noppt') && ~isopt(opts,'autoword')
    ppt = guidata(pptplot);
    %fileStart = file{1}; fileEnd = file{end};
    %indEnd = strfind(fileStart,'\');    fileStart = fileStart(indEnd(end)+1:end);
    %indEnd = strfind(fileEnd,'\'); fileEnd = fileEnd(indEnd(end)+1:end);   
    %set(ppt.e_file,'String',file{1});
    %set(ppt.e_figures,'String',['[',sprintf('%d ',fignum),']']);
    %set(ppt.e_title,'String',[fileStart ' - ' fileEnd]);
    set(ppt.exported,'Value',0);
elseif isopt(opts,'autoword')
    pptinfo.body = str;
    scanfile=file{1};
    f = fignums;
    pptinfo.scanfile = scanfile;
    pptinfo.configch = d{1}.configch; pptinfo.scan = d{1}.scan; pptinfo.configvals = d{1}.configvals; 
    save2pptman(pptinfo,f)
end