function [file,out] = plotChrg(opts,file,config)
% plot 2D non pulsed data.
% function file = plotChrgB(opts,file,cond)
% default is to plot both data and differentiated data, do line plot of all
% files when you click color plot.
% possible options:
%   sens: fits slope of line for line cut.
%   sens, both: both plots both data and diff data. sens plots diff data.
%   glitch: remove data more than a couple SD away from mean.
%   ydiff : differentiate data in y direction
%   nolow : don't do low pass filtering on diff plot.
%   flat: remove polynomial from data -- useful for sensing
%   two: assume there are two data channels (e.g. for charge sensing and CB)
%   square: make square plot
%   cond: plot conductance of data
%   log: plot log of data
%   hyst: clicking on last plot plots line cuts of all data sets. (default)
%   chron: plot files in chronological order
%   water: waterfall plot.
%   next: option for loading incrementing set of scans.
%   cbar: Use small colorbars, try to use a single color ball for everything.
%   zero : make all negative data a NaN. Useful for Ithaco offset. 
%   replot: still in process. want to only reload last file when refresh.
%   pick: use fancy filter for picking files.
%   invis: make all plots invisible. Makes it easier when ppting data.
%   autoplot: automatically plot to ppt. This makes figures invisible, doesn't show chanDisp.
%   noppt: don't bring up ppt gui
%   flip : analysis for shortlived attempt to flip scan dir after each line to speed up.
%   probably hysteresis makes this untenable.
%   comp: compare
%   stop: don't run both/hyst automatically.
%   noise: analyze noise in scan
%can press keys on figure to cause functions to run: (uses func keyPressed)
%   r: refresh file (useful when actively taking data)
%	m: measure distance between, slope of CB peaks
%	s: Run spyview on data
%	x: Load next set of scans in folder
%	l: Drag a rectangle on screen and replot just that data in new figure.
%	(useful for when color bar saturated).
%   j: Drag a rectangle on screen and replot just that data in same figure.
%	d: Click twice, gives distance between points.
%	p: Click once, give single point location
%   c: doublePt: Plot clicked point as dot on both diff and charge (figures 1 and 3)
%   f: Change button down function to plot line from figures 1 and 3 as double Y Plt.
%   e: hyst: change button down function to hyst.
% can give a conditions struct, 'cond' with fields
%	vout, i0 to calculate conductance
%   sens/chrg which are indices you want plotted for diff data, normal data.
%   glitch: # standard deviations to get rid of
%   filt: a filter pattern, such as sensL
%   sortConf / precis,
% other notes:
%     rotate

%% Configure options, load file, set up figures.
if ~exist('opts','var'), opts = ''; end
if ~isopt(opts,'stop'), opts = [opts ' both hyst']; end
if ~exist('config','var') || isempty(config), config = struct(); end
if isopt(opts,'sens') % If sens, button down function fits sensor slope.
    btnFcn = @btnFit;
else
    btnFcn = @btn;
end
if isopt(opts,'both'),  opts = [opts ' chrg sens']; end % give both differentiated and normal data.
if isopt(opts,'next') % allow one to sort through data.
    fileFolder = dir;
    dates = [fileFolder(:).datenum]; % you may want to eliminate . and .. first.
    [~,inds] = sort(dates);
    fileNames = {fileFolder(inds).name}; % Cell array of names in order by datenum.
    fileToke=regexp([file{:}], '.*?(\d{4}).mat','tokens'); fileToke = [fileToke{:}];
    fileNums = cell2mat(fileToke'); [~,sortInds]=sort(str2double(string(fileNums)));
    lastFile = find(strcmp(file{sortInds(end)},fileNames));
    file = fileNames(lastFile+1:lastFile + length(file));
end 
if (~exist('file','var') || isempty(file)) && ~isopt(opts, 'pick')% grab files: if chron, show in chronological order.
    if ~isfield(config,'filt')
        [file,fpath]=getFiles('sm*.mat');
    else
        [file,fpath] = getFiles(['sm_*' config.filt '*']);
    end
    if ~isopt(opts,'chron'),	file = fliplr(file);    end
elseif (~exist('file','var') || isempty(file)) && isopt(opts, 'pick')
    file = uipickfiles('FilterSpec','sm_*.mat');
end
if file{1}==0, return;  end
nfiles = length(file);
xLabOff = [32,16,14]; yLabOff = [23 11 7];
cBarPos = [0.03 0.015 0.01]; pos = [0.85 0.375 0.208]; %85/94
if ~isopt(opts,'water')
    if nfiles<=2
        nrow =1; ncol = nfiles;
        plotSpace = {nrow,ncol, [0.063, 0.12], [0.06 0.045], [0.06, 0.1]};
    elseif nfiles <= 4
        nrow = 2; ncol = 2;
        plotSpace = {nrow,ncol, [0.063, 0.12], [0.06 .045], [0.06, 0.1]};
    elseif nfiles<=6
        nrow = 2; ncol = 3;
        plotSpace = {nrow,ncol, [0.063, 0.073], [0.06 .045], [0.06, 0.1]};
    else
        nrow = 3; ncol = 3;
        plotSpace = {nrow,ncol, [0.063, 0.073], [0.06 0.045], [0.06, 0.1]};
    end
else
    if nfiles <= 2
        subPlotArg = [1,nfiles];
    elseif nfiles <= 4
        subPlotArg = [2,2];
    else
        subPlotArg = [3, ceil(nfiles/3)];
    end
end
if isopt(opts,'sens') && ~isopt(opts,'replot')
    fSens = makeFigure(3);
    if isopt(opts,'invis')
        fSens.Visible = 'off';
    else
        fSens.Visible = 'on';
        figure(3);
    end
    if ~isopt(opts,'autoplot'), fSens.Name = ['Diff ' file{1}(4:end-4)]; end
    ha = tight_subplot(plotSpace{:});
end
if isopt(opts,'chrg') && ~isopt(opts,'replot')
    fChrg = makeFigure(1);
    if isopt(opts,'invis')
        fChrg.Visible = 'off';
    else
        fChrg.Visible = 'on';
        figure(1);
    end
    if ~isopt(opts,'autoplot'),fChrg.Name = file{1}(4:end-4); end
    ga = tight_subplot(plotSpace{:});
end
if isopt(opts,'comp')
    fComp=figure(6); clf; fComp.Name = file{1}(4:end-4);
    ia = tight_subplot(nfiles,2,[0.1 0.05],0.05,[0.06 0.1]);
end
if isopt(opts,'water')
    fWater = makeFigure(7);
    if isopt(opts,'invis')
        fWater.Visible = 'off';
    else
        fWater.Visible = 'on';
    end
    if ~isopt(opts,'autoplot'),fWater.Name = file{1}(4:end-4); end
    ka = tight_subplot(subPlotArg(1),subPlotArg(2),[0.08 0.06],0.05,[0.06 0.1]);
end
%% Analyze data
if isopt(opts,'replot'), file = file(1); end
n=1;
str =[]; chgstr = [];
for i = 1:length(file) % Check that file exists, has non NaN data, and is 2D.
    if exist('fpath','var')
        try
            dtmp=load([fpath file{i}]);
        catch
            break
        end
    else
        dtmp=load([file{i}]);
    end
    fileInd = strfind(file{i},'\');
    if ~isempty(fileInd)
        fileName0{i} = file{i}(fileInd(end)+4:end-4);
    else
        fileName0{i} = file{i}(4:end-4);
    end
    if ~isfield(dtmp,'data') || isempty(dtmp.data) || all(isnan(dtmp.data{1}(:))) || length(dtmp.scan.loops(1).rng)==1 || ndims(dtmp.data{1})>2  %#ok<*ISMAT>
        str = [str fileName0{i} ' bad file ' newline];
    else
        fileList{n}=file{i};
        fileName{n} = fileName0{i};
        d(n)=dtmp;
        n=n+1;
    end
end
if isfield(config,'sortConf')
    if isfield(config,'chan')
        if isfield(config,'precis')
            [d,fileList] = sortConfig(d,fileList,config.sortConf,config.chan,config.precis);
        else
            [d,fileList] = sortConfig(d,fileList,config.sortConf,config.chan);
        end
    else
        [d,fileList] = sortConfig(d,fileList,config.sortConf);
    end
end
if ~exist('fileList','var')
    warning('No data in files')
    return
end
if length(fileList)> 9 % only use first 9 files.
    fileList = fileList(1:9);
end
for i=1:length(fileList)
    data = d(i).data{1}; scan = d(i).scan;
    if i == 1
        slideInfo.configch = d(i).configch; slideInfo.scan = d(i).scan; slideInfo.configvals = d(i).configvals;
        if ~isopt(opts,'autoplot') && ~isopt(opts,'invis')
            chanDisp(d(i).configch,d(i).configvals);
        end
    end
    if isfield(d(i).scan,'data') && isfield(d(i).scan.data,'flip') && d(i).scan.data.flip
        for j = 2:2:size(data,1)
            data(j,:) = fliplr(data(j,:));
        end
        rng1 = scanfn(d(i).scan,1);
        d(i).scan.loops(1).rng = rng1(2,:);
        rng2 = scanfn(d(i).scan,2);
        d(i).scan.loops(2).rng = rng2(2,:);
        d(i).scan.loops(2).setchan = d(i).scan.loops(2).setchan{2:end};
        d(i).scan.loops(1).setchan = d(i).scan.loops(1).setchan{2:end};
    end
    if isfield(d(i).scan,'data') && isfield(d(i).scan.data,'opts') && strcmp(d(i).scan.data.opts,'rotate')
        scantmp = d(i).scan;
        d(i).scan.loops(1).setchan = scantmp.loops(2).setchan;
        d(i).scan.loops(2).setchan = scantmp.loops(1).setchan;
        d(i).scan.loops(1).rng = scantmp.loops(2).rng;
        d(i).scan.loops(2).rng = scantmp.loops(1).rng;
        d(i).scan.loops(2).npoints = scantmp.loops(1).npoints;
        d(i).scan.loops(1).npoints = scantmp.loops(2).npoints;
        data = data';
        scan = d(i).scan;
    end
    if isopt(opts,'zero')
        if nanmean(data(:))<0
            data =-data;
        end
        data(sign(data)==-1)=NaN;
    end
    %% Print info on scan.
    pointSpacingX = diff(d(i).scan.loops(1).rng)/d(i).scan.loops(1).npoints;
    pointSpacingY = diff(d(i).scan.loops(2).rng)/d(i).scan.loops(2).npoints;
    ramprate = abs(pointSpacingX/d(i).scan.loops(1).ramptime);
    str = [str, sprintf('%s: ', fileName{i})];
    str =[str, sprintf('Spacing: X = %3.3f mV Y = %3.3f mV. Ramprate: %3.0f mV/s', abs(pointSpacingX(1))*1e3, abs(pointSpacingY(1))*1e3, ramprate(1)*1e3)];
    if isopt(opts,'buff')
        if d(i).scan.loops(1).ramptime < 0
            str = [str, ' Self-Ramping.'];
        end
        if ~isempty(strfind(d(i).scan.loops(2).getchan,'buf'))
            str = [str, ' Buffered'];
        end
    end
    if length(d(i).data)>1,	str = [str, '>1 chan'];     end
    if d(i).scan.loops(1).rng(1) > d(i).scan.loops(1).rng(2)
        str = [str, ' R to L. '];
    else
        str = [str, ' L to R. '];
    end
    [good,rng]=smramprate(d(i).scan); %decide later if you want to use range. Seems to depend on computer.
    if ~good
        %scan.loops(1).rng = rng;
        str = [str 'Fast. '];
    end
    str = sprintf('%s \n',str);
    xvals=scanRng(scan,1); yvals=scanRng(scan,2);
    cache(i).data = data; cache(i).xvals = xvals; cache(i).yvals = yvals; cache(i).file=fileList{i}; % set up cache for button down functions.
    xlab = makeLabel(d(i).scan.loops(1).setchan); ylab = makeLabel(d(i).scan.loops(2).setchan);
    if i > 1
        chgstrNew = changeConfigGen(d(i),oldconfig,oldconfigch,fileName{i});
        chgstr = [chgstr chgstrNew];
    end
    if ischar(d(i).scan.loops(2).getchan) || isempty(d(i).scan.loops(2).getchan)
        d(i).scan.loops(2).getchan = {d(i).scan.loops(2).getchan};
    end
    if ~isempty(d(i).configvals)
        oldconfig = d(i).configvals;
        oldconfigch = d(i).configch;
    else
        oldconfig = [];
        oldconfigch = [];
    end
    %% Plot everything
    if isopt(opts,'sens')
        if isopt(opts,'two') || isfield(config,'sens')
            if isfield(config,'sens')
                data = d(i).data{config.sens};
            else
                data = d(i).data{1};
            end
        end
        if ~isopt(opts,'ydiff')
            dataDiff=diff(data,[],2);
            fsigma=[0.4 0.1];
            if ~isopt(opts,'nolow')
                dataFilt = lowpass(fsigma,dataDiff);
            else
                dataFilt = dataDiff;
            end
            dataNorm = dataFilt./diff(xvals);
        else
            dataDiff=diff(data);
            dataNorm = dataDiff./diff(yvals)';
        end
        m=nanmean(dataNorm(:)); s=nanstd(dataNorm(:));
        if isopt(opts,'glitch')
            if ~isfield(config,'glitch') || isempty(config.glitch), 	config.glitch = 6;            end
            dataNorm(abs(dataNorm)-m>config.glitch*s)=NaN;
        end
        cacheSens(i).xvals = 1/2*(xvals(1:end-1)+xvals(2:end)); cacheSens(i).yvals = yvals; cacheSens(i).data = dataNorm; cacheSens(i).file =fileList{i}; %#ok<*AGROW>
        if isopt(opts,'ydiff')
            cacheSens(i).xvals = xvals; cacheSens(i).yvals = 1/2*(yvals(1:end-1)+yvals(2:end));
        end
        if ~isopt('opts','replot')
            a = ha(i);
        else
            fSens = figure(3);
            axesInds = find(isgraphics(fSens.Children,'axes'));
            a = axesInds(1);
        end
        if i == length(fileList), 	fSens.WindowKeyPressFcn=@(h,e) keyPressed(e,opts,fileList,i,config);        end
        g=imagesc(a,xvals,yvals,dataNorm,'ButtonDownFcn',@btn);  a.YDir = 'norm';
        a.Title.Interpreter = 'None'; a.Title.String = fileName{i};
        a.YLabel.Position(1) = a.XLim(1) - range(a.XLim)/yLabOff(ncol);
        a.XLabel.Position(2) = a.YLim(1) - range(a.YLim)/xLabOff(nrow);
        a.XLabel.String =xlab;         a.YLabel.String =ylab;
        if ~isopt(opts,'cbar')
            axPos = a.Position;
            colorbar(a,'Position',[axPos(1) + pos(ncol)+cBarPos(ncol), axPos(2), 0.015, axPos(4)], 'Location','manual');
            a.Position(3) = pos(ncol);
        end
        a.YDir = 'norm'; axis tight;
        if isopt(opts,'square'), 	axis equal;        end
        if isopt(opts,'hyst'), 	g.ButtonDownFcn = @(src,clk) hyst(src,clk,cacheSens);        end
    end
    if isopt(opts,'chrg')
        if isopt(opts,'two') || isfield(config,'chrg')
            if isfield(config,'chrg')
                data = d(i).data{config.chrg};
            else
                data = d(i).data{2};
            end
        end
        [data, newStr] = formatData(data,config,xvals,opts);
        a=ga(i);
        cache(i).data = data;
        h=imagesc(a,xvals,yvals,data,'ButtonDownFcn',btnFcn);
        if ~isopt(opts,'cbar')
            axPos = a.Position;
            colorbar(a,'Position',[axPos(1) + pos(ncol)+cBarPos(ncol), axPos(2), 0.015, axPos(4)], 'Location','manual');
            a.Position(3) = pos(ncol);
        end
        a.Title.Interpreter = 'None'; a.Title.String = fileName{i};
        a.YLabel.Position(1) = a.XLim(1) - range(a.XLim)/yLabOff(ncol);
        a.XLabel.Position(2) = a.YLim(1) - range(a.YLim)/xLabOff(nrow);
        a.XLabel.String =xlab; a.YLabel.String =ylab;
        a.YDir = 'norm';
        if isopt(opts,'hyst'), 	h.ButtonDownFcn = @(src,clk) hyst(src,clk,cache);  end
        fChrg.WindowKeyPressFcn=@(h,e) keyPressed(e,opts,fileList,i,config);
        axis tight;
        if isopt(opts,'square'), 	axis equal;        end
    end
    if isopt(opts,'water')
        a = ka(i);
        dataInds = find(~all(isnan(data')));
        data = data(dataInds,:);
        nLines = size(data,1);
        yleg = cellstr(num2str(yvals')); loc = 'northwest';
        if ~isempty(strfind(fileList{i},'oneGateL'))
            yleg = {'2a','2b','1a','1b','2a, 2b','1a, 1b'};
            yleg = yleg(dataInds);
            loc= 'southeast';
        elseif ~isempty(strfind(fileList{i},'oneGateR'))
            yleg = {'3a','3b','4a','4b','3a, 3b','4a, 4b'};
            yleg = yleg(dataInds);
        end
        for k = 1 : nLines
            plot(a,xvals,data(k,:),'DisplayName',yleg{k}); hold(a,'on');
        end
        a.Title.Interpreter = 'None'; a.Title.String = fileName{i};
        a.XLabel.String =xlab; a.YLabel.String =ylab;
        axis tight;
        l=legend(a,'show');
        l.Location = loc;
    end
    if isopt(opts,'comp')
        figure(6); a=ia(2*i-1);
        imagesc(a,xvals,yvals,dataNorm,'ButtonDownFcn',@btn);  a.YDir = 'norm';
        a=ia(2*i);
        if length(d(i).data)>1
            data = d(i).data{2};
        else
            warning('Only one data chan, can''t run comp')
            continue
        end
        h=imagesc(a,xvals,yvals,data,'ButtonDownFcn',@btn);  a.YDir = 'norm';
        cacheComp(1).data = dataNorm; cacheComp(1).xvals = xvals(1:end-1); cacheComp(1).yvals = yvals;
        cacheComp(2).data = data; cacheComp(2).xvals = xvals; cacheComp(2).yvals = yvals;
        cacheComp(1).fileList = fileList{i}; cacheComp(2).fileList = fileList{i};
        h.ButtonDownFcn = @(src,clk) doubleYPlt(src,clk,cacheComp);
    end
    %% Other analysis
    if isopt(opts,'noise')
        Fs = abs(1./ d(i).scan.loops(1).ramptime);
        psdData=NaN(size(data,1),floor(size(data,2)/2)+1);
        
        for j = 1:length(yvals)
            [sigf,freqs]=psd(data(j,:),Fs); %#ok<*DPSD>
            psdData(j,:)=sigf;
        end
        figure(44); clf;
        meanPsd = nanmean(psdData);
        plot(freqs(3:end),sqrt(meanPsd(3:end)),'.')
        xlabel('Frequency (Hz)'); ylabel('Noise')
        approxNoise = sum(meanPsd(3:end));
        fprintf('Noise rms level is %3.3f \n',approxNoise);
    end
    if isopt(opts,'maxSlp')
        midPt = round(scan.loops(2).npoints/2);
        midVal = yvals(midPt);
        [maxSens,ind]=max(dataNorm(midPt,:));
        maxVal = 1/2*(xvals(ind)+xvals(ind+1));
        out = struct('sens',maxSens,'yval',midVal,'xval',maxVal);
        return;
    end
end
%% Configure Plotting
if isopt(opts,'cbar') && isopt(opts,'chrg')
    betterCbar(fChrg,length(fileList),ncol)
end
if isopt(opts,'cbar') && isopt(opts,'sens')
    betterCbar(fSens,length(fileList),ncol)
end
if isopt(opts,'sens') && isopt(opts,'chrg')
    f = [1,3];
elseif isopt(opts,'sens')
    f = 3;
elseif isopt(opts,'chrg')
    f = 1;
end
if isopt(opts,'water')
    f = 7;
end
if ~isopt(opts,'noform')
    formatFig(f(1),'chrg');
    formatFig(f(2:end),'sens');
end
if ~isopt(opts, 'noppt') && ~isopt(opts,'autoppt') % PPT gui
    ppt = guidata(pptplot);
    set(ppt.e_file,'String',fileList{1});
    if isopt(opts,'sens') && isopt(opts,'chrg')
        set(ppt.e_figures,'String',['[',sprintf('%d %d',1,3),']']);
    elseif isopt(opts,'sens')
        set(ppt.e_figures,'String',['[',sprintf('%d ',3),']']);
    elseif isopt(opts,'chrg')
        set(ppt.e_figures,'String',['[',sprintf('%d ',1),']']);
    end
    if isopt(opts,'water')
        set(ppt.e_figures,'String',['[',sprintf('%d ',4),']']);
    end
    set(ppt.e_title,'String','');
    set(ppt.e_body,'String',[chgstr str]);
    set(ppt.exported,'Value',0);
elseif isopt(opts,'autoppt')
    slideInfo.body = str;
    slideInfo.body2 = chgstr;
    slideInfo.comments = '';    slideInfo.title = ''; scanfile=fileList{1};
    slideInfo.scanfile = scanfile;
    if isopt(opts,'difffirst')
        f = fliplr(f);
    end
    save2pptauto(slideInfo,f)
end
end