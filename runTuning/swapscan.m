function swapscan(setname,opts)
%  swap scandata sides in/out
%    swapscan()             -- list available
%    swapscan('foo','save') -- save current scandata as foo
%    swapscan('bar')        -- save current scandata, load bar
% opts: 
%   save
%   rename
%   load
global scandata;
nocopy={'sets','file'}; % these do not belong to sides, so don't copy
if ~exist('setname','var')
    fprintf('Available sets: (%s currently loaded)\n',scandata.name);
    fprintf('\t%s\n',scandata.sets.name);
else
    if ~exist('opts','var'), opts='load'; end
    if isempty(setname), setname=scandata.name; end
    ind=find(strcmp({scandata.sets.name},setname)); % index of qubit
    switch opts 
        case 'rename'
            scandata.name=setname;
            scandata.sets(ind).name = setname; 
        case 'save' % Copy all current fields into sets(oldInd)
            scandata.name=setname;
            flds=fields(rmfield(scandata,nocopy));
            if(isempty(ind))
                ind=length(scandata.sets)+1;
            end
            for i=1:length(flds)
                scandata.sets(ind).(flds{i})=scandata.(flds{i});
            end
        case 'load'
            swapscan(scandata.name,'save'); 
            if isempty(ind)
                error('Unable to load set "%s", doesn''t exist yet. \n',setname);                
            end
            % Merge set data into scandata instead of replacing scandata,
            % so we don't accidentally delete any data. 
            mergeset=scandata.sets(ind);
            flds = fieldnames(mergeset);
            for i=1:length(flds)
                scandata.(flds{i})= mergeset.(flds{i});
            end
            fprintf('Loaded scandata set "%s"\n',scandata.name);
    end
end