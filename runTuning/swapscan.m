function swapscan(setname,opts)
%  swap scandata sides in/out
%    swapscan()             -- list available
%    swapscan('foo','save') -- save current scandata as foo
%    swapscan('bar')        -- save current scandata, load bar

global scandata;
nocopy={'sets'}; % these do not belong to sides, so don't copy
if ~exist('setname','var')
    fprintf('Available sets: (%s currently loaded)\n',scandata.name);
    fprintf('\t%s\n',scandata.sets.name);
else
    if ~exist('opts','var')
        opts='load';
    end
    if isempty(setname)
        setname=scandata.name;
    end
    i=find(strcmp({scandata.sets.name},setname));
    switch opts 
        case 'rename'
            scandata.name=setname;
            scandata.sets(i).name = setname; 
        case 'save'
            scandata.name=setname;
            f=fields(rmfield(scandata,nocopy));
            if(isempty(i))
                i=length(scandata.sets)+1;
            end
            for j=1:length(f)
                scandata.sets(i).(f{j})=scandata.(f{j});
            end
        case 'load'
            swapscan(scandata.name,'save');
            if isempty(i)
                error('Unable to load set "%s"\n',setname);                
            end
            % Merge set data into scandata so we don't nuke fields we neither copy nor understand..
            mergeset=scandata.sets(i);
            f = fieldnames(mergeset);
            for j=1:length(f)
                scandata.(f{j})= mergeset.(f{j});
            end
            fprintf('Loaded scandata set "%s"\n',scandata.name);
    end
end