function [file,out] = qpcPlot(opts,file,bigtitle)
% Plot a bunch of separate qpc files
% function file = qpcPlot(opts,file)
% opts: 
%   hyst: plot ramprate, zero val. 
%   num: plot file number in legend (default is date)
%   noise
%   autoword: automatically plot
%   noppt: no ppt gui
%   invis: don't show plot
%	cond: plot conductance.  
% bigtitle just puts a title over the top of the image. maybe stupid. 

global pptdata

if ~exist('opts','var')
    opts = '';
end
if ~exist('file','var') || isempty(file)
    [file,fpath]=get_files('sm_qpc*.mat');
    file0=file;
    if ~strcmp(fpath,pwd)         
        file = fullfile(fpath,file); 
    end    
else
    fpath = [];
    for i = 1 :length(file)
        foldInd = strfind(file{i},'\'); foldInd = [0 foldInd]; 
        file0{i} = file{i}(foldInd(end)+1:end);
    end
end % Pick files
if file{1}==0
    return 
end

minValList = [];
for i = 1:length(file) %Load data, get file info
    d{i} = load(file{i});
    setchan{i}=d{i}.scan.loops(1).setchan;
    minValList = [minValList min(d{i}.scan.loops(1).rng)];
    if ~isopt(opts,'num') %figure out what folder
        if isfield(pptdata,'qpcFolder') && (~isempty(strfind(file{i},pptdata.qpcFolder)) || ~isempty(strfind(fpath,pptdata.qpcFolder)))
            fileInd = strcmp({pptdata.qpcDir.name},file0{i});
            d{i}.name = pptdata.qpcDir(fileInd).date([1:6 12:17]);
        elseif (isfield(pptdata,'dataFolder') && ~isempty(strfind(file{i},pptdata.dataFolder))) || isempty(fpath) || (~isempty(strfind(fpath,pptdata.dataFolder)))
            fileInd = strcmp({pptdata.dir.name},file0{i});
            d{i}.name = pptdata.dir(fileInd).date([1:6 12:17]);
        elseif (isempty(strfind(file{i},'/')) && isempty(fpath)) || strcmp(fpath,pwd)
            if ~exist('dirInfo','var')
                dirInfo = dir;
            end
            fileInd = strcmp({dirInfo.name},file0{i});
            d{i}.name = dirInfo(fileInd).date([1:6 12:17]);
        end
    end
end
minValList = unique(minValList);
filesRem = 1:length(file);  % This is just the number of files. 
n = 1;
while length(filesRem) >= 1 % Find all unique sets of setchans, their matches, and keep removing files from list until list empty. 
    cellequal = @(x) isequal(setchan{filesRem(1)},x); % Function finds which setchans equal to given one. 
    inds=find(cellfun(cellequal,setchan));
    for i =1:length(inds)
        filesRem(filesRem==inds(i))=[]; % remove all the files that correspond to a match from the set of files. 
    end
    sameset{n} = inds;
    n=n+1;
end
str = ''; str2 = ''; h=[];
npl = length(sameset); % Find the number of figures have
nrow = ceil(npl/3); nfig = ceil(nrow/3); nrow = ceil(nrow/nfig);
fignum = 5:5+nfig-1; 
for i =1:nfig % Create figures 
    currFig = makeFigure(fignum(i));
    if isopt(opts,'invis')
        currFig.Visible = 'off';
    else
        currFig.Visible = 'on';
    end
    ga{i} = tight_subplot(nrow,3,[0.1 0.05],[0.07 0.05],0.03);
    f = gca; colors = f.ColorOrder;
    if exist('bigtitle','var') && ~isempty(bigtitle) && i == 1
        n = mtit(bigtitle); n.th.Position(2) = n.th.Position(2)+0.02;
    end
end
if isopt(opts,'hyst') % Create figures for plotting hyst info
    fignumH = fignum + 100;
    fignumH2 = fignum + 200;
    for i =1:nfig
        currFig = makeFigure(fignumH(i));
        if isopt(opts,'invis')        
        currFig.Visible = 'off';
        else 
            currFig.Visible = 'on'; 
        end
        ha{i} = tight_subplot(nrow,3,[0.1 0.05],[0.07 0.05],0.03); haNum(i) = 0; 
        currFig = makeFigure(fignumH2(i)); 
        if isopt(opts,'invis')        
            currFig.Visible = 'off';
        else 
            currFig.Visible = 'on'; 
        end
        ia{i} = tight_subplot(nrow,3,[0.1 0.05],[0.07 0.05],0.03);
    end
end

for i = 1:length(sameset) % Analyze, plot data. each samset index is set of setchans 
    condQuant{i} = nan(1,length(sameset{i})); zeroVal{i} = condQuant{i}; rampRate = condQuant{i}; minVal{i} = condQuant{i}; buff = condQuant{i}; % make them a set of nans. 
    xlab=d{sameset{i}(1)}.scan.loops(1).setchan; % get figure stuff ready 
    if ischar(xlab)
        xlab = {xlab};
    end
    for k=1:size(xlab,2)-1
        xlab{k}=[xlab{k} ', '];
    end
    xlab = cell2mat(xlab);    
    figind = ceil(i/nrow/3);
    a = ga{figind}(mod(i-1,nrow*3)+1); 
    for j = 1:length(sameset{i}) % each index j is diff scan with same setchans.
        c = sameset{i}(j); %index of scan in terms of original file list.
        buff(j) = isempty(d{c}.scan.loops(1).getchan) && ~isempty(d{c}.scan.loops(2).getchan); %check if data was buffered.
        xvals=linspace(d{c}.scan.loops(1).rng(1),d{c}.scan.loops(1).rng(2),d{c}.scan.loops(1).npoints);
        minVal{i}(j) = min(xvals);
        if isfield(d{c},'data')
            if isopt(opts,'cond') && length(d{c}.data)>1 % plot the conductance
                data = d{c}.data{2};
                data(data == -Inf) = NaN;
            else
                data = abs(d{c}.data{1});
            end
        else
            continue
        end
        if size(data,2) > 1 && size(data,1) > 1 % Skip file if data 2d.
            continue
        end
        if isopt(opts,'filt') && any(abs(data(:))>1e-1)
            data = 1e-8*data;
        end
        if length(d{c}.data) > 1 % if we've recorded can use that to find qoc.
            if any(d{c}.data{2}<0) % quantum of conductance data
                [~,mi]= min(abs(d{c}.data{2}));
                condQuant{i}(j) = xvals(mi); % row is which scan, j is which gate.
            end
        else
            if any(d{c}.data{1}<2e-9)
                [~,mi]= min(abs(abs(d{c}.data{1})-2e-9));
                condQuant{i}(j) = xvals(mi);
            end
        end
        if any(d{c}.data{1}<1e-10) % record voltage for zero current.
            [~,mi]= min(abs(abs(d{c}.data{1})-1e-10));
            zeroVal{i}(j) = xvals(mi);
        end
        if isopt(opts,'hyst') %
            rampRate(j) = diff(d{c}.scan.loops(1).rng)/d{c}.scan.loops(1).npoints/abs(d{c}.scan.loops(1).ramptime);
            dispName = sprintf('%1.0f, %s', 1e3*rampRate(j),d{c}.name(1:6));%file{c}(end-7:end-4));
        elseif isopt(opts,'num')
            dispName = file0{c}(end-7:end-4);
        else
            dispName = d{c}.name;
        end
        if length(d{c}.configch)==2 || ~isempty(strfind(file{c},'4k')) % for a while at least 4k station had 2 configchans. probably want to add in name filter.
            dispName = [dispName ' 4k'];
        end                      
        h(end+1) = plot(a,xvals,data,'DisplayName',dispName); hold(a,'on');
        if ~isempty(d{c}.configvals) && length(d{c}.configvals)>=18 && any(abs(d{c}.configvals(1:18))>1e-3) % print nonzero channels .
            inds = find(abs(d{c}.configvals(1:18))>=1e-3);
            str2 = [str2 sprintf('%s\t',file0{c}(4:end-4))];
            for m = 1:length(inds)
                str2 = [str2 sprintf('%s: %3.3f, ',d{c}.configch{inds(m)},d{c}.configvals(inds(m)))];
            end
            str2 = [str2 sprintf('\n')];
        end
    end
    a.XLabel.String =xlab;
    name = file0{c}; 
    if ~isempty(strfind(name,'\'))
        ind = strfind(name, '\'); 
        name = name(ind+1:end); 
    end
    a.Title.Interpreter = 'None';
    a.Title.String = name(8:end-4); 
    if isopt(opts,'noppt') 
        l = clickableLegend(a,'show'); 
    else
        l=legend(a,'show');         
    end
    l.Interpreter='none'; l.Location = 'northwest';
    l.FontSize = 7;
    mV(i) = nanmean(condQuant{i}); 
    mVZ(i) = nanmean(zeroVal{i});
    str = [str sprintf('%s : \t %3.3f   %3.3f   %3.3f  \n',xlab,mV(i),mVZ(i),mVZ(i)-mV(i))]; % This prints out all the zero data.           
    if isopt(opts,'hyst')  % Plot condQuant and zeroVal data against ramprate. 
        out.chans{i} = xlab; 
        if length(unique(rampRate))>1
            haNum(figind) = 1;
            a = ha{figind}(mod(i-1,nrow*3)+1);
            for j = 1:length(rampRate)
                ind(j) = find(minVal{i}(j)==minValList);
                if ~buff(j)
                    plot(a,1000*rampRate(j), condQuant{i}(j)-mV(i),'*','Color',colors(ind(j),:)); hold(a,'on');
                else
                    plot(a,1000*rampRate(j), condQuant{i}(j)-mV(i),'o','Color',colors(ind(j),:)); hold(a,'on');
                end
            end
            a.XBaseline.Visible = 'on'; a.YBaseline.Visible = 'on'; xlimVal = a.XLim;
            a.XLimMode = 'manual'; a.XLim = xlimVal; 
            a.Title.String = [xlab ': ' num2str(mV(i),'%3.3f')];
            a.XLabel.String =' Ramprate (mV /s)';
                        
            a =ia{figind}(mod(i-1,nrow*3)+1);
            for j = 1:length(rampRate)
                if ~buff(j)
                    plot(a,1000*rampRate(j), zeroVal{i}(j)-mVZ(i),'*','Color',colors(ind(j),:)); hold(a,'on');
                else
                    plot(a,1000*rampRate(j), zeroVal{i}(j)-mVZ(i),'o','Color',colors(ind(j),:)); hold(a,'on');
                end
            end
            a.XLabel.String =' Ramprate (mV /s)';
            a.Title.String = [xlab ': ' num2str(mVZ(i),'%3.3f')];
            a.XBaseline.Visible = 'on'; a.YBaseline.Visible = 'on'; xlimVal = a.XLim; 
            a.XLimMode = 'manual'; a.XLim = xlimVal; 
        end        
    end    
end

out.condQuant = condQuant; out.zeroVal = zeroVal;
if isopt(opts,'noise')
    figure(44); clf; hold(a,'on');
    for i =1:length(d)
        Fs = abs(1./ d{i}.scan.loops(1).ramptime);
        [psdData,freqs]=psd(d{i}.data{1},Fs);
        plot(freqs(3:end),sqrt(psdData(3:end)),'.-')
        pause
    end
    xlabel('Frequency (Hz)'); ylabel('Noise')
end
if isopt(opts,'hyst') && any(haNum)
    plotFigs = find(haNum); 
    a=ha{plotFigs(end)}(end); 
    for i = 1:length(minValList)
        text(a,0.25,0.23+i*0.086,num2str(minValList(i)),'Color',colors(i,:))
    end            
    fignum = [fignum fignumH(plotFigs) fignumH2(plotFigs)];
end
formatFig(fignum,'qpc');
if ~isopt(opts,'noppt') && ~isopt(opts,'autoword')
    ppt = guidata(pptplot3);
    fileStart = file0{1}; fileEnd = file0{end};    
    set(ppt.e_file,'String',file{1});
    set(ppt.e_figures,'String',['[',sprintf('%d ',fignum),']']);
    set(ppt.e_title,'String',[fileStart ' - ' fileEnd]);
    set(ppt.exported,'Value',0);
    set(ppt.e_body,'String',str)
elseif isopt(opts,'autoword')
    slideInfo.body = str;    slideInfo.body2 = str2; 
    slideInfo.scanfile=file{1};    
    slideInfo.configch = d{1}.configch; slideInfo.scan = d{1}.scan; slideInfo.configvals = d{1}.configvals;     
    save2pptauto2(slideInfo,fignum)
end
end