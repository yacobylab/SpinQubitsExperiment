function out = plotOne(opts,file,config)
% For plotting oneD pulsed scans -- Mostly autotune.
% function out = plotOne(opts,file,config)
if ~exist('opts','var'),    opts = ''; end
if ~exist('config','var') || isempty(config),   config = struct(); end
if isopt(opts,'typ')
    config.filt = {'*stp*';'*tl*';'*load*';'*lead*'};
end
if (~exist('file','var') || isempty(file)) % grab files: if chron, show in chronological order.
    if ~isfield(config,'filt')
        [file,fpath]=getFiles('sm*.mat');
    else        
        [file,fpath] = getFiles(config.filt);
    end
    if ~isopt(opts,'chron'), file = fliplr(file); end
end
if file{1}==0, return;  end

%%
n=1; chgstr = '';
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
        if n>1
            chgstrNew = changeConfigGen(d(i),oldconfig,oldconfigch,fileName{n});
            chgstr = [chgstr chgstrNew];
        end
        if ~isempty(d(n).configvals)
            oldconfig = d(n).configvals;
            oldconfigch = d(n).configch;
        else
            oldconfig = [];
            oldconfigch = [];
        end
        n=n+1;
    end
end
if ~exist('fileList','var')
    warning('No data in files')
    return
end

nPlots = ceil(length(fileList)/12);
axesList = [];
plotSpace = {4,3, [0.063, 0.073], [0.06 0.045], [0.06, 0.1]};
for i = 1:nPlots
    figure(12+i); clf;
    axesList = [axesList, tight_subplot(plotSpace{:})];
end
for i = 1:length(fileList)
    a = axesList(i);
    data = nanmean(d(i).data{1});
    if isopt(opts,'center')
        data=data-mean(data);
    end
    if ~isopt(opts,'nox') && isfield(d(i).scan.data.pulsegroups(1),'varpar')
        xvals = d(i).scan.data.pulsegroups(1).varpar(:,1)';
        plot(a,xvals,data);
    else
        plot(a,data);
    end
    a.Title.Interpreter = 'None'; a.Title.String = fileName{i};
    %a.XLabel.String =xlab;
end
fprintf(chgstr); 
end